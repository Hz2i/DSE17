# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import brentq

# SPEED + ALTITUDE
v_inf = 27.6  # Freestream velocity in m/s
altitude = 18288  # 60,000 ft operating altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
rho = atmo.density() # (kg/m^3)
speed_of_sound = atmo.speed_of_sound()
dynamic_viscosity = atmo.dynamic_viscosity() # Added for accurate Re calculation

# Propeller Constants
Nb = 2            # Number of blades
D = 2.0           # Propeller Diameter (m)
R_abs = D / 2.0   # Tip radius (m)

airfoil_names = ["SD7037"]

# Expanded bounds for aggressive maneuvers / deep stall
alphas_sweep = np.linspace(-35, 35, 250)

# Radial geometry: Stopping at 0.995 to prevent Prandtl tip singularity (Divide by zero)
r = np.linspace(0.3, 0.8, 100)  
r_abs = r * R_abs                
beta_07 = 20.5                   

# GEOMETRY FROM PAPER
b = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D
beta_deg = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)

# advance ratio sweep
lambda_J = np.linspace(0.3, 1.0, 30)

# storage for results per airfoil
results = {}

for name in airfoil_names:
    print(f"\nEvaluating Airfoil: {name}...")
    try:
        airfoil = asb.Airfoil(name)
    except Exception:
        print(f"Warning: failed to load airfoil {name}, skipping.")
        continue

    ct_curve = []
    cp_curve = []
    eta_curve = []

    for J in lambda_J:
        n_rps = v_inf / (D * J)
        omega = 2 * np.pi * n_rps
        phi_0_rad = np.arctan(v_inf / (omega * r_abs))

        # 1. DYNAMIC KINEMATICS (Fixing the hardcoded Re/Mach issue)
        v_local = np.sqrt(v_inf**2 + (omega * r_abs)**2)
        Re_dist = (rho * v_local * b) / dynamic_viscosity
        mach_dist = v_local / speed_of_sound

        dT_dr = np.zeros_like(r_abs)
        dM_dr = np.zeros_like(r_abs)
        
        bemt_failures = 0 # Track failures per J sweep

        for i in range(len(r_abs)):
            try:
                # Dynamically fetch aero data for THIS specific element's Re and Mach
                aero_data = airfoil.get_aero_from_neuralfoil(
                    alpha=alphas_sweep,
                    Re=Re_dist[i],
                    mach=mach_dist[i]
                )

                # 3. LINEAR EXTRAPOLATION (Fixing the dangerous cubic splines)
                cl_interp = interp1d(alphas_sweep, aero_data['CL'], kind='linear', bounds_error=False, fill_value="extrapolate")
                cd_interp = interp1d(alphas_sweep, aero_data['CD'], kind='linear', bounds_error=False, fill_value="extrapolate")

                # EQ 26 PAPER
                def equilibrium_equation(alpha_guess):
                    cl_val = cl_interp(alpha_guess)
                    lhs = cl_val * (Nb * b[i]) / (2 * np.pi * r_abs[i])
                    rhs = 4 * np.sin(np.radians(beta_deg[i] - alpha_guess)) * \
                              np.tan(np.radians(beta_deg[i] - alpha_guess) - phi_0_rad[i])
                    return lhs - rhs

                # Expanded search space to match linear extrapolation
                true_alpha_deg = brentq(equilibrium_equation, -30, 30)
                true_phi_rad = np.radians(beta_deg[i] - true_alpha_deg)

                cl_local = cl_interp(true_alpha_deg)
                cd_local = cd_interp(true_alpha_deg)

                # EQ 18 PAPER
                K = (cl_local * Nb * b[i] * np.cos(true_phi_rad)) / (8 * np.pi * r_abs[i] * (np.sin(true_phi_rad)**2))
                K = np.clip(K, -0.5, 0.5)
                one_plus_a = 1.0 / (1.0 - K)

                velocity_term = (one_plus_a**2) / (np.sin(true_phi_rad)**2)

                # 4. TIP LOSS REFINEMENT (Naturally protected by max r=0.995)
                f_prandtl = (Nb / 2.0) * (R_abs - r_abs[i]) / (R_abs * np.sin(true_phi_rad))
                F_local = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))

                dT_dr[i] = b[i] * velocity_term * (cl_local * np.cos(true_phi_rad) - cd_local * np.sin(true_phi_rad)) * F_local
                dM_dr[i] = b[i] * velocity_term * (cl_local * np.sin(true_phi_rad) + cd_local * np.cos(true_phi_rad)) * r_abs[i] * F_local

            except (ValueError, ZeroDivisionError):
                # 2. FAILURE TRACKING (Instead of silent fails)
                bemt_failures += 1
                dT_dr[i] = 0.0
                dM_dr[i] = 0.0    
        
        if bemt_failures > 5:
            print(f"  -> Warning: BEMT failed for {bemt_failures} elements at J={J:.2f}. Expect C_T/C_P dips.")

        # Integration
        F0 = 0.5 * rho * (v_inf**2) * Nb
        Total_Thrust = F0 * trapezoid(dT_dr, r_abs)
        Total_Torque = F0 * trapezoid(dM_dr, r_abs)
        
        C_T = Total_Thrust / (rho * (n_rps**2) * (D**4))
        C_P = (Total_Torque * omega) / (rho * (n_rps**3) * (D**5))
        
        eta = (J * C_T) / C_P if C_P > 0 else 0.0
        
        ct_curve.append(C_T)
        cp_curve.append(C_P)
        eta_curve.append(eta)

    results[name] = {
        'J': lambda_J,
        'CT': ct_curve,
        'CP': cp_curve,
        'ETA': eta_curve
    }


# Plotting: compare airfoils
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
colors = plt.cm.tab10.colors

for idx, name in enumerate(results.keys()):
    ax1.plot(results[name]['J'], results[name]['ETA'], '-', color=colors[idx % len(colors)], linewidth=2, label=f"{name}")
ax1.set_xlabel('Advance Ratio (J)', fontsize=11)
ax1.set_ylabel('Propeller Efficiency ($\\eta$)', fontsize=11)
ax1.set_title('Propeller Efficiency Curve Comparison', fontsize=12, weight='bold')
ax1.set_xlim(0.3, 1)
ax1.set_ylim(0.2, 1.0)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='lower left')

for idx, name in enumerate(results.keys()):
    ax2.plot(results[name]['J'], results[name]['CT'], 'o-', color=colors[idx % len(colors)], linewidth=1.5, label=f"{name} C_T", markersize=4)
    ax2.plot(results[name]['J'], results[name]['CP'], 's--', color=colors[idx % len(colors)], linewidth=1.5, label=f"{name} C_P", markersize=4)
ax2.set_xlabel('Advance Ratio (J)', fontsize=11)
ax2.set_ylabel('Coefficient Value', fontsize=11)
ax2.set_title('Thrust and Power Coefficients vs. Advance Ratio', fontsize=12, weight='bold')
ax2.set_xlim(0.3, 0.85)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper right', fontsize='small')

plt.tight_layout()
plt.show()

# print optimal advance ratio for sd7037
if 'SD7037' in results:
    optimal_idx = np.argmax(results['SD7037']['ETA'])
    print(f"\nOptimal Advance Ratio for SD7037: J = {results['SD7037']['J'][optimal_idx]:.3f} with Efficiency = {results['SD7037']['ETA'][optimal_idx]:.3f}")

# plot beta over radius
plt.figure(figsize=(8, 5))
plt.plot(r, beta_deg, 'b-', linewidth=2)
plt.xlabel('Normalized Radius (r/R)', fontsize=11)
plt.ylabel('Blade Pitch Angle (degrees)', fontsize=11)
plt.title('Blade Pitch Distribution Along the Span', fontsize=12, weight='bold')
plt.grid(True, linestyle=':', alpha=0.6)
plt.show()