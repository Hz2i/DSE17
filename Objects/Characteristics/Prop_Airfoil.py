# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt


altitude = 20000  # 20 km altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
reynolds_number = 100000
mach_number = 0.15

# Angle of Attack sweep range
alphas = np.linspace(-4, 14, 100) 

airfoil_names = ["S1223", "NACA4412", "E387", "SD7037", "FX63137"]
results = {}

print(f"Running aerodynamic sweeps at Re = {reynolds_number:,}...")

for name in airfoil_names:
    formatted_name = "fx63137" if name == "FX63137" else name.lower()
    
    try:
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
    except Exception as e:
        print(f"Error loading or analyzing {name}: {e}")


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
plt.tight_layout()
plt.show()