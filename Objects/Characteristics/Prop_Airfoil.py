# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.interpolate import interp1d

v_inf = 30.0  # Freestream velocity in m/s (typical for HAPS)
altitude = 20000  # 20 km altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
reynolds_number = 200000
mach_number = v_inf/atmo.speed_of_sound()
Nb = 3 # Number of blades, needed for blade geometry calculations
# print(f"At {altitude} m: Re = {reynolds_number:,}, Mach = {mach_number:.2f}")
# rho = atmo.density()
# print(f"Air density at {altitude} m: {rho:.4f} kg/m³")

# Angle of Attack sweep range
alphas = np.linspace(-4, 14, 100) 

airfoil_names = ["S1223", "NACA4412", "E387", "SD7037", "FX63137"]
results = {}

print(f"Running aerodynamic sweeps at Re = {reynolds_number:,}...")

# ===========================================================================
# SOURCE: Performance Calculation and Design of Stratospheric Propeller 
# ===========================================================================
r = np.linspace(0.28, 0.9, 90)  # Normalized radial position along the blade
# r = 0.7  # Focus on 70% radius for blade geometry calculations
D = 2.5  # Assume a reference diameter of 1 meter for scaling
beta_07 = 20.5  # Initial blade angle at 70% radius (degrees)
rpm = 960  # Rotational speed in RPM
n = rpm / 60  # Convert RPM to revolutions per second

b = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D
beta = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)  # Convert to degrees
phi_0 = np.arctan2(v_inf, 2 * np.pi * n * (r * D/2)) * (180/np.pi)



# ===========================================================================
# Airfoil results
# ===========================================================================
for name in airfoil_names:
    airfoil = asb.Airfoil(name)
    
    # Pull viscous aerodynamics data directly using the integrated NeuralFoil engine
    aero_data = airfoil.get_aero_from_neuralfoil(
        alpha=alphas,
        Re=reynolds_number,
        mach=mach_number
    )
    
    # Calculate Cl/Cd ratio cleanly
    aero_data['Cl_Cd'] = aero_data['CL'] / aero_data['CD']
    
    results[name] = aero_data
    print(f"Successfully analyzed: {name}")


# # plot the blade geometry
# fig, ax = plt.subplots(figsize=(10, 5))
# # pitch and chord distribution
# ax.plot(r*D, b, label='Chord Length (b)', linewidth=2)
# ax.plot(r*D, beta, label='Blade Angle (beta)', linewidth=2)
# ax.plot(r*D, phi_0, label='Inflow Angle (phi_0)', linewidth=2)
# ax.set_xlabel('Radial Position (m)')
# ax.set_title('Stratospheric Propeller Blade Geometry')
# ax.grid(True, linestyle=':', alpha=0.6)
# plt.show()
# # ax.legend()


# 1. Create an interpolation function for your NeuralFoil Cl data
# so the solver can look up Cl for ANY angle of attack instantly.
cl_interp = interp1d(alphas, results['S1223']['CL'], kind='cubic', bounds_error=False, fill_value="extrapolate")

# # plot the Cl interpolation to verify it looks correct
# alpha_test = np.linspace(-4, 14, 200)
# plt.figure(figsize=(8, 5))
# plt.plot(alphas, results['S1223']['CL'], 'o', label='NeuralFoil Data (S1223)', markersize=5)
# plt.plot(alpha_test, cl_interp(alpha_test), '-', label='Cubic Interpolation', linewidth=2)
# plt.show()


# 2. Define the paper's Vortex Theory root-finding equation
def interference_eq(alpha_deg, r_local, b_local, beta_deg, phi0_deg, n_rpm):
    alpha_rad = np.radians(alpha_deg)
    beta_rad = np.radians(beta_deg)
    phi0_rad = np.radians(phi0_deg)
    
    # Get Cl from our NeuralFoil interpolation
    cl = cl_interp(alpha_deg) 

    # # Prevent division by zero near the hub
    # if np.abs(np.sin(phi0_rad)) < 1e-6:
    #     F = 1.0
    # else:
    #     # Eq 34 from the paper
    #     f_prandtl = (Nb / 2) * (D/2 - r_local * D/2) / (D/2 * np.sin(phi0_rad))
    #     f_prandtl = np.clip(f_prandtl, 0, 100) # Prevent math overflow at the exact tip
    #     # Eq 33 from the paper
    #     F = (2 / np.pi) * np.arccos(np.exp(-f_prandtl))
    
    # The paper's f(alpha) = 0 equation 
    # Note: absolute radius used here too! absolute_r = r_local * (D/2)
    left_side = cl * (Nb * b_local) / (2 * np.pi * n_rpm/60 * (r_local * D/2)) # PLEASE CHECK THIS FORMULA!!
    right_side = 4 * np.sin(beta_rad - alpha_rad) * np.tan(beta_rad - alpha_rad - phi0_rad)
    
    return left_side - right_side

# 3. Loop through your blade stations and solve!
actual_alphas = []
for i in range(len(r)):
    # Find the angle of attack (between -4 and 14 degrees) that makes the equation = 0
    solved_alpha = brentq(interference_eq, -4, 14, args=(r[i], b[i], beta[i], phi_0[i], rpm))
    actual_alphas.append(solved_alpha)

alpha_i = beta - phi_0 - actual_alphas

# test formula with plot
# formula = cl_interp(alphas) * (Nb * b) / (2 * np.pi * n) - 4 * np.sin(np.radians(beta) - np.radians(alphas)) * np.tan(np.radians(beta) - np.radians(alphas) - np.radians(phi_0))

# plot formula vs alpha for all values of r (which affects b, beta, and phi_0)
# plt.figure(figsize=(10, 6))
# for i in range(len(r)):
#     formula_i = cl_interp(alphas) * (Nb * b[i]) / (2 * np.pi * n * (r[i]*D/2)) - 4 * np.sin(np.radians(beta[i]) - np.radians(alphas)) * np.tan(np.radians(beta[i]) - np.radians(alphas) - np.radians(phi_0[i]))
#     plt.plot(alphas, formula_i, label=f'r={r[i]:.2f}', linewidth=2)
# plt.axhline(0, color='red', linestyle='--', label='Zero Line')
# plt.xlabel(r'Angle of Attack, $\alpha$ (degrees)')
# plt.title('Interference Equation Value vs. Angle of Attack for Different Radial Positions')
# plt.grid(True, linestyle=':', alpha=0.6)
# plt.legend()
# plt.show()

# plot actual alphas vs r
plt.figure(figsize=(8, 5))
plt.plot(r/D, actual_alphas, label='Solved Angle of Attack', linewidth=2)
plt.xlabel('Radial Position (m)')
plt.ylabel(r'Angle of Attack, $\alpha$ (degrees)')
plt.title('Solved Angle of Attack vs. Radial Position on Propeller Blade')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.show()

# plot alpha i vs r
plt.figure(figsize=(8, 5))
plt.plot(r/D, alpha_i, label=r'$\alpha_i$ (Induced Angle)', linewidth=2)
plt.xlabel('Radial Position (m)')
plt.ylabel(r'Induced Angle, $\alpha_i$ (degrees)')
plt.title('Induced Angle vs. Radial Position on Propeller Blade')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.show()


# plt.figure(figsize=(8, 5))
# plt.plot(alphas, formula, label='Interference Equation Value', linewidth=2)
# plt.axhline(0, color='red', linestyle='--', label='Zero Line')
# plt.xlabel(r'Angle of Attack, $\alpha$ (degrees)')
# plt.title('Interference Equation Value vs. Angle of Attack')
# plt.grid(True, linestyle=':', alpha=0.6)
# plt.legend()
# plt.show()




# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

# # Plot 1: Lift Curve (Alpha vs. Cl)
# for name, data in results.items():
#     ax1.plot(alphas, data['CL'], label=name, linewidth=2)
# ax1.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')  # Added 'r' for raw string
# ax1.set_ylabel(r'Lift Coefficient, $C_l$')
# ax1.set_title(r'Lift Curves ($\alpha$ vs. $C_l$)')      # Fixed: Added 'r' here
# ax1.grid(True, linestyle=':', alpha=0.6)
# ax1.legend()

# # Plot 2: Drag Polar (Cd vs. Cl)
# for name, data in results.items():
#     ax2.plot(data['CD'], data['CL'], label=name, linewidth=2)
# ax2.set_xlabel(r'Drag Coefficient, $C_d$')
# ax2.set_ylabel(r'Lift Coefficient, $C_l$')
# ax2.set_title(r'Drag Polars ($C_d$ vs. $C_l$)')        # Added 'r' for consistency
# ax2.grid(True, linestyle=':', alpha=0.6)
# ax2.legend()

# # Plot 3: Aerodynamic Efficiency (Alpha vs. Cl/Cd)
# for name, data in results.items():
#     ax3.plot(alphas, data['Cl_Cd'], label=name, linewidth=2.5)
# ax3.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')  # Added 'r' for raw string
# ax3.set_ylabel(r'Aerodynamic Efficiency, $C_l/C_d$')
# ax3.set_title(r'Section Efficiency ($\alpha$ vs. $C_l/C_d$)') # Fixed: Added 'r' here
# ax3.grid(True, linestyle=':', alpha=0.6)
# ax3.legend()

# plt.suptitle(f"HAPS Airfoil Performance & Efficiency Matrix ($Re={reynolds_number:,}$ at 20 km)", fontsize=14, weight='bold')
# plt.tight_layout()

# # Additional plot: S1223 across multiple Reynolds numbers
# re_list = [50000, 100000, 200000, 500000, 1000000]
# s1223 = asb.Airfoil("S1223")
# s1223_results_re = {}
# print(f"Running S1223 Reynolds sweep for Re = {[int(r) for r in re_list]}...")
# for Re in re_list:
#     aero = s1223.get_aero_from_neuralfoil(
#         alpha=alphas,
#         Re=Re,
#         mach=mach_number
#     )
#     aero['Cl_Cd'] = aero['CL'] / aero['CD']
#     s1223_results_re[Re] = aero

# fig2, (ax4, ax5) = plt.subplots(1, 2, figsize=(14, 5))

# # Lift curves for S1223 at different Re
# for Re, data in s1223_results_re.items():
#     ax4.plot(alphas, data['CL'], label=f"Re={int(Re):,}")
# ax4.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')
# ax4.set_ylabel(r'Lift Coefficient, $C_l$')
# ax4.set_title(r'S1223: Lift Curves at Different $Re$')
# ax4.grid(True, linestyle=':', alpha=0.6)
# ax4.legend()

# # Efficiency curves for S1223 at different Re
# for Re, data in s1223_results_re.items():
#     ax5.plot(alphas, data['Cl_Cd'], label=f"Re={int(Re):,}")
# ax5.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')
# ax5.set_ylabel(r'Aerodynamic Efficiency, $C_l/C_d$')
# ax5.set_title(r'S1223: Efficiency ($C_l/C_d$) at Different $Re$')
# ax5.grid(True, linestyle=':', alpha=0.6)
# ax5.legend()

# plt.suptitle(f"S1223 Performance Across Reynolds Numbers (mach={mach_number})", fontsize=13)
# plt.tight_layout()
# plt.show()