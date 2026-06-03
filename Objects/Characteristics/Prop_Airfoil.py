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
mach_number = v_inf / speed_of_sound

# Propeller Constants
Nb = 2            # Number of blades
D = 2.0           # Propeller Diameter (m)
R_abs = D / 2.0   # Tip radius (m)

airfoil_names = ["S1223", "NACA4412", "E387", "SD7037", "FX63137"] 

alphas_sweep = np.linspace(-30, 30, 250)

# radial geometry and blade distributions
r = np.linspace(0.1, 1, 100)  # Normalized radius (r/R)
r_abs = r * R_abs                # Absolute radius in meters
beta_07 = 20.5                   # Blade pitch setting at 70% radius (degrees)

# GEOMETRY FROM PAPER
b = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D
beta_deg = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)  # Convert to degrees


# advance ratio sweep
lambda_J = np.linspace(0.3, 1.0, 30)

# storage for results per airfoil
results = {}

for name in airfoil_names:
    try:
        airfoil = asb.Airfoil(name)
        aero_data = airfoil.get_aero_from_neuralfoil(
            alpha=alphas_sweep,
            Re=100_000,
            mach=mach_number
        )

        cl_interp = interp1d(alphas_sweep, aero_data['CL'], kind='cubic', bounds_error=False, fill_value="extrapolate")
        cd_interp = interp1d(alphas_sweep, aero_data['CD'], kind='cubic', bounds_error=False, fill_value="extrapolate")

    except Exception:
        # if neuralfoil lookup fails, skip this airfoil
        print(f"Warning: failed to get aero data for {name}, skipping.")
        continue

    # Compute performance curves for this airfoil
    ct_curve = []
    cp_curve = []
    eta_curve = []

    for J in lambda_J:
        # Compute the unique rotational speed matching this specific advance ratio
        n_rps = v_inf / (D * J)
        omega = 2 * np.pi * n_rps
        # Compute uninduced geometric inflow angle distribution
        phi_0_rad = np.arctan(v_inf / (omega * r_abs))

        # Arrays to accumulate the continuous force integrals across the span
        dT_dr = np.zeros_like(r_abs)
        dM_dr = np.zeros_like(r_abs)

        for i in range(len(r_abs)):
            # EQ 26 PAPER
            def equilibrium_equation(alpha_guess):
                cl_val = cl_interp(alpha_guess)
                lhs = cl_val * (Nb * b[i]) / (2 * np.pi * r_abs[i])
                rhs = 4 * np.sin(np.radians(beta_deg[i] - alpha_guess)) * \
                          np.tan(np.radians(beta_deg[i] - alpha_guess) - phi_0_rad[i])
                return lhs - rhs

            try:
                # Solve for the true operating angle of attack
                true_alpha_deg = brentq(equilibrium_equation, -25, 25)
                true_phi_rad = np.radians(beta_deg[i] - true_alpha_deg)

                # Extract true viscous forces at this operating point
                cl_local = cl_interp(true_alpha_deg)
                cd_local = cd_interp(true_alpha_deg)

                # EQ 18 PAPER
                K = (cl_local * Nb * b[i] * np.cos(true_phi_rad)) / (8 * np.pi * r_abs[i] * (np.sin(true_phi_rad)**2))
                # potential localized mathematical boundaries
                K = np.clip(K, -0.5, 0.5)
                # The true axial induction multiplier
                one_plus_a = 1.0 / (1.0 - K)

                # vector projection
                velocity_term = (one_plus_a**2) / (np.sin(true_phi_rad)**2)


                # local tip loss factor for torque and thrust
                f_prandtl = (Nb / 2.0) * (R_abs - r_abs[i]) / (R_abs * np.sin(true_phi_rad))
                f_prandtl = np.clip(f_prandtl, 0.0, 100.0)
                F_local = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))


                dT_dr[i] = b[i] * velocity_term * (cl_local * np.cos(true_phi_rad) - cd_local * np.sin(true_phi_rad)) * F_local
                dM_dr[i] = b[i] * velocity_term * (cl_local * np.sin(true_phi_rad) + cd_local * np.cos(true_phi_rad)) * r_abs[i] * F_local

            except (ValueError, ZeroDivisionError):
                # if stall or numerical issue
                dT_dr[i] = 0.0
                dM_dr[i] = 0.0    
       
        # Integration
        F0 = 0.5 * rho * (v_inf**2) * Nb
        Total_Thrust = F0 * trapezoid(dT_dr, r_abs)
        Total_Torque = F0 * trapezoid(dM_dr, r_abs)
        
        # Non-Dimensional Performance Coefficients
        C_T = Total_Thrust / (rho * (n_rps**2) * (D**4))
        C_P = (Total_Torque * omega) / (rho * (n_rps**3) * (D**5))
        
        # efficiency
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

# Plot 1: Efficiency Curve for all airfoils
for idx, name in enumerate(results.keys()):
    ax1.plot(results[name]['J'], results[name]['ETA'], '-', color=colors[idx % len(colors)], linewidth=2, label=f"{name}")
# ax1.axvline(0.75, color='gray', linestyle='--', alpha=0.7, label='Design Point (J=0.75)')
ax1.set_xlabel('Advance Ratio (J)', fontsize=11)
ax1.set_ylabel('Propeller Efficiency ($\\eta$)', fontsize=11)
ax1.set_title('Propeller Efficiency Curve Comparison', fontsize=12, weight='bold')
ax1.set_xlim(0.3, 1)
ax1.set_ylim(0.2, 1.0)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='lower left')

# Plot 2: Dimensionless Coefficients Curve for all airfoils
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