# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize_scalar
import warnings

warnings.filterwarnings("ignore")

# ===========================================================================
# 1. TAKE-OFF FLIGHT ENVIRONMENT & AIRCRAFT CONFIG
# ===========================================================================
altitude = 0      # Sea level for TO density
num_engines = 4   # The 4 engines on the aircraft

# Drivetrain Efficiencies
eta_motor = 0.90
eta_esc = 0.99
eta_elec = eta_motor * eta_esc

atmo = asb.Atmosphere(altitude=altitude)
rho = atmo.density() 
speed_of_sound = atmo.speed_of_sound()

m_TO = 198
S = 41.5
s_TOG = 20     # Take-off ground distance [m]
C_L_TO = 1.0   # Take-off lift coefficient (assumed)
f_LW = 1.2

# ROSKAM Equations
v_TO = f_LW * np.sqrt(2 * 9.81 * m_TO / (C_L_TO * rho * S))  
P_TO_Total_Watts = ((m_TO * 9.81)**2) / (s_TOG * rho * S * C_L_TO)  

print("\n=======================================================")
print(f"Estimated Take-Off Speed : {v_TO:.2f} m/s")
print(f"Roskam Power Limit (Max) : {P_TO_Total_Watts/1000:.2f} kW")
print("=======================================================")

# ===========================================================================
# 2. PROPELLER GEOMETRY (From Sizing Results)
# ===========================================================================
Nb = 2            # Number of blades per propeller
D = 2.1357        # Propeller Diameter (m)
R_abs = D / 2.0   # Tip radius (m)
airfoil_name = "SD7037"

# Radial geometry and blade distributions
r = np.linspace(0.1, 1.0, 100)   # Normalized radius (r/R)
r_abs = r * R_abs                # Absolute radius in meters
beta_07 = 20.5                   # Blade pitch setting at 70% radius (degrees)

# Chord and Twist distributions
b = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D
beta_deg = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)

# ===========================================================================
# 3. SEA-LEVEL AERODYNAMICS (DEEP STALL)
# ===========================================================================
alphas_sweep = np.linspace(-30, 85, 300) 
airfoil = asb.Airfoil(airfoil_name)

print("Fetching viscous aerodynamic data for Sea Level...")
aero_data = airfoil.get_aero_from_neuralfoil(
    alpha=alphas_sweep,
    Re=1_000_000,   
    mach=0.05
)

cl_interp = interp1d(alphas_sweep, aero_data['CL'], kind='linear', bounds_error=False, fill_value="extrapolate")
cd_interp = interp1d(alphas_sweep, aero_data['CD'], kind='linear', bounds_error=False, fill_value="extrapolate")

# ===========================================================================
# 4. RPM SWEEP & BEMT SOLVER
# ===========================================================================
# Sweep from a low idle up to your target maximum 1800 RPM
rpms_sweep = np.linspace(100, 1800, 30)

total_thrusts = []
total_elec_powers_kw = []
efficiencies = []

print(f"Running Take-Off RPM Sweep (100 to 1800 RPM) at {v_TO:.2f} m/s...")

for rpm in rpms_sweep:
    n_rps = rpm / 60.0
    omega = 2 * np.pi * n_rps
    J = v_TO / (n_rps * D)

    phi_0_rad = np.arctan(v_TO / (omega * r_abs))

    dT_dr = np.zeros_like(r_abs)
    dM_dr = np.zeros_like(r_abs)

    for i in range(len(r_abs)):
        def equilibrium_equation(alpha_guess):
            cl_val = cl_interp(alpha_guess)
            phi_guess_rad = np.radians(beta_deg[i] - alpha_guess)
            lhs = cl_val * (Nb * b[i]) / (2 * np.pi * r_abs[i])
            rhs = 4 * np.sin(phi_guess_rad) * np.tan(phi_guess_rad - phi_0_rad[i])
            return lhs - rhs

        try:
            # Expanded search brackets for steep root stall angles
            true_alpha_deg = brentq(equilibrium_equation, -25, 85)
        except ValueError:
            # Fallback for extreme deep stall
            res = minimize_scalar(lambda a: abs(equilibrium_equation(a)), bounds=(-25, 85), method='bounded')
            true_alpha_deg = res.x

        try:
            true_phi_rad = np.radians(beta_deg[i] - true_alpha_deg)
            cl_local = cl_interp(true_alpha_deg)
            cd_local = cd_interp(true_alpha_deg)

            K = (cl_local * Nb * b[i] * np.cos(true_phi_rad)) / (8 * np.pi * r_abs[i] * (np.sin(true_phi_rad)**2))
            K_effective = np.clip(K, -0.5, 0.5)
            one_plus_a = 1.0 / (1.0 - K_effective)
            velocity_term = (one_plus_a**2) / (np.sin(true_phi_rad)**2)

            f_prandtl = (Nb / 2.0) * (R_abs - r_abs[i]) / (R_abs * np.sin(true_phi_rad))
            f_prandtl = np.clip(f_prandtl, 0.0, 100.0)
            F_local = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))

            dT_dr[i] = b[i] * velocity_term * (cl_local * np.cos(true_phi_rad) - cd_local * np.sin(true_phi_rad)) * F_local
            dM_dr[i] = b[i] * velocity_term * (cl_local * np.sin(true_phi_rad) + cd_local * np.cos(true_phi_rad)) * r_abs[i] * F_local

        except ZeroDivisionError:
            dT_dr[i] = 0.0
            dM_dr[i] = 0.0    

    # -----------------------------------------------------------------------
    # 5. INTEGRATION & 4-ENGINE SCALING
    # -----------------------------------------------------------------------
    F0 = 0.5 * rho * (v_TO**2) * Nb
    Thrust_Per_Prop = F0 * trapezoid(dT_dr, r_abs)
    Torque_Per_Prop = F0 * trapezoid(dM_dr, r_abs)
    Mech_Power_Per_Prop = Torque_Per_Prop * omega

    # Total Aircraft Metrics (x4 Engines)
    Total_Aircraft_Thrust = Thrust_Per_Prop * num_engines
    Total_Aircraft_Mech_Power = Mech_Power_Per_Prop * num_engines
    
    # Calculate True Electrical Power Draw
    Total_Aircraft_Elec_Power_kW = (Total_Aircraft_Mech_Power / eta_elec) / 1000.0

    # Dimensionless Coefficients
    C_T = Thrust_Per_Prop / (rho * (n_rps**2) * (D**4))
    C_P = Mech_Power_Per_Prop / (rho * (n_rps**3) * (D**5))

    efficiency = (J * C_T) / C_P if C_P > 0 else 0.0
    
    total_thrusts.append(Total_Aircraft_Thrust)
    total_elec_powers_kw.append(Total_Aircraft_Elec_Power_kW)
    efficiencies.append(efficiency * 100.0)

# ===========================================================================
# 6. PLOTTING THE SWEEP RESULTS
# ===========================================================================
fig, axs = plt.subplots(1, 3, figsize=(16, 5))

# Plot 1: Total Thrust vs RPM
axs[0].plot(rpms_sweep, total_thrusts, 'b-', linewidth=2.5)
axs[0].set_title(f'Total Aircraft Thrust (4 Engines)\nat {v_TO:.2f} m/s Runway Speed', fontsize=11)
axs[0].set_xlabel('Rotational Speed (RPM)')
axs[0].set_ylabel('Total Thrust (N)')
axs[0].grid(True, linestyle=':', alpha=0.6)

# Plot 2: Total Electrical Power vs RPM with Roskam Limit
axs[1].plot(rpms_sweep, total_elec_powers_kw, 'r-', linewidth=2.5, label='Required Elec Power')
axs[1].axhline(y=P_TO_Total_Watts/1000.0, color='k', linestyle='--', linewidth=2, label='Roskam Power Limit')
axs[1].set_title(f'Total Required Electrical Power\nat {v_TO:.2f} m/s Runway Speed', fontsize=11)
axs[1].set_xlabel('Rotational Speed (RPM)')
axs[1].set_ylabel('Total Power (kW)')
axs[1].grid(True, linestyle=':', alpha=0.6)
axs[1].legend(loc='upper left')

# Plot 3: Aerodynamic Efficiency vs RPM
axs[2].plot(rpms_sweep, efficiencies, 'g-', linewidth=2.5)
axs[2].set_title(f'Propeller Efficiency\nat {v_TO:.2f} m/s Runway Speed', fontsize=11)
axs[2].set_xlabel('Rotational Speed (RPM)')
axs[2].set_ylabel('Efficiency (%)')
axs[2].grid(True, linestyle=':', alpha=0.6)

plt.suptitle(f"Take-Off Performance Sweep | D = {D:.2f} m | Sea Level", fontsize=14, weight='bold')
plt.tight_layout()
plt.show()