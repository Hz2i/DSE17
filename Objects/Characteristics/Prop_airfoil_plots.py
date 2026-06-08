import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt

# Constant parameters
v_TAS = 27.0 # m/s, from preliminary sizing
reynolds_number = 100000 # approximate
altitude = 18288  # m, 60000 feet altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
mach_number = v_TAS/atmo.speed_of_sound()

# Blade geometry


# Airfoil analysis
airfoil_names = ["S1223", "NACA4412", "E387", "SD7037", "FX63137"]
alphas = np.linspace(-80, 80, 100)  # Angle of attack
results = {}
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
    print(f"  Max Cl/Cd: {np.max(aero_data['Cl_Cd']):.2f} at alpha = {alphas[np.argmax(aero_data['Cl_Cd'])]:.2f} degrees")

# Plot Cl/Cd vs. alpha for all airfoils
plt.figure(figsize=(10, 6))
plt.plot(alphas, results["S1223"]['Cl_Cd'], label="S1223")
plt.plot(alphas, results["NACA4412"]['Cl_Cd'], label="NACA4412")
plt.plot(alphas, results["E387"]['Cl_Cd'], label="E387")
plt.plot(alphas, results["SD7037"]['Cl_Cd'], label="SD7037")
plt.plot(alphas, results["FX63137"]['Cl_Cd'], label="FX63137")
plt.xlabel("Angle of Attack (degrees)")
plt.ylabel("Cl/Cd Ratio")
plt.title(f"Cl/Cd vs. Angle of Attack for Various Airfoils at Re={reynolds_number}")
plt.legend()
plt.grid()
plt.show()

# CL vs. alpha plot
plt.figure(figsize=(10, 6))
plt.plot(alphas, results["S1223"]['CL'], label="S1223")
plt.plot(alphas, results["NACA4412"]['CL'], label="NACA4412")
plt.plot(alphas, results["E387"]['CL'], label="E387")
plt.plot(alphas, results["SD7037"]['CL'], label="SD7037")
plt.plot(alphas, results["FX63137"]['CL'], label="FX63137")
plt.xlabel("Angle of Attack (degrees)")
plt.ylabel("Lift Coefficient (CL)")
plt.title(f"CL vs. Angle of Attack for Various Airfoils at Re={reynolds_number}")
plt.legend()
plt.grid()
plt.show()

# CD vs. alpha plot
plt.figure(figsize=(10, 6))
plt.plot(alphas, results["S1223"]['CD'], label="S1223")
plt.plot(alphas, results["NACA4412"]['CD'], label="NACA4412")
plt.plot(alphas, results["E387"]['CD'], label="E387")
plt.plot(alphas, results["SD7037"]['CD'], label="SD7037")
plt.plot(alphas, results["FX63137"]['CD'], label="FX63137")
plt.xlabel("Angle of Attack (degrees)")
plt.ylabel("Drag Coefficient (CD)")
plt.title(f"CD vs. Angle of Attack for Various Airfoils at Re={reynolds_number}")
plt.legend()
plt.grid()
plt.show()
