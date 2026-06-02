# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt


altitude = 20000  # 20 km altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
reynolds_number = 200000
mach_number = 25.0/atmo.speed_of_sound()
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
r = np.linspace(0, 1, 100)  # Normalized radial position along the blade
D = 1.0  # Assume a reference diameter of 1 meter for scaling
beta_07 = 20.0  # Initial blade angle at 70% radius (degrees)

b = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D
beta = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)  # Convert to degrees


# # plot the blade geometry
# fig, ax = plt.subplots(figsize=(10, 5))
# # pitch and chord distribution
# ax.plot(r*D, b, label='Chord Length (b)', linewidth=2)
# ax.plot(r*D, beta, label='Blade Angle (beta)', linewidth=2)
# ax.set_xlabel('Radial Position (m)')
# ax.set_title('Stratospheric Propeller Blade Geometry')
# ax.grid(True, linestyle=':', alpha=0.6)
# ax.legend()

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


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

# Plot 1: Lift Curve (Alpha vs. Cl)
for name, data in results.items():
    ax1.plot(alphas, data['CL'], label=name, linewidth=2)
ax1.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')  # Added 'r' for raw string
ax1.set_ylabel(r'Lift Coefficient, $C_l$')
ax1.set_title(r'Lift Curves ($\alpha$ vs. $C_l$)')      # Fixed: Added 'r' here
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend()

# Plot 2: Drag Polar (Cd vs. Cl)
for name, data in results.items():
    ax2.plot(data['CD'], data['CL'], label=name, linewidth=2)
ax2.set_xlabel(r'Drag Coefficient, $C_d$')
ax2.set_ylabel(r'Lift Coefficient, $C_l$')
ax2.set_title(r'Drag Polars ($C_d$ vs. $C_l$)')        # Added 'r' for consistency
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend()

# Plot 3: Aerodynamic Efficiency (Alpha vs. Cl/Cd)
for name, data in results.items():
    ax3.plot(alphas, data['Cl_Cd'], label=name, linewidth=2.5)
ax3.set_xlabel(r'Angle of Attack, $\alpha$ (degrees)')  # Added 'r' for raw string
ax3.set_ylabel(r'Aerodynamic Efficiency, $C_l/C_d$')
ax3.set_title(r'Section Efficiency ($\alpha$ vs. $C_l/C_d$)') # Fixed: Added 'r' here
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.legend()

plt.suptitle(f"HAPS Airfoil Performance & Efficiency Matrix ($Re={reynolds_number:,}$ at 20 km)", fontsize=14, weight='bold')
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