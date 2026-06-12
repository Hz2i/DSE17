import numpy as np
import matplotlib.pyplot as plt
import aerosandbox as asb
import csv

# ==========================================
# INPUT PARAMETERS
# ==========================================
D_ref = 1.7890472741146266*1000   # Reference Diameter in meters
beta_07 = 20.5   # Pitch angle at 75% station in degrees 

num_stations = 20
r_stations = np.linspace(0.0, 1.0, num_stations) 

# ==========================================
# AIRFOIL COORDINATES
# ==========================================
print("Fetching SD7037 coordinates...")
sd7037 = asb.Airfoil("sd7037").repanel(n_points_per_side=50)

x_base = sd7037.coordinates[:, 0] - 0.25 # Shift to quarter-chord
y_base = sd7037.coordinates[:, 1]


catia_rows = []
catia_rows.append(['StartLoft', '', ''])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

print(f"Generating {num_stations} 3D sections for D = {D_ref}m...")

for idx, r in enumerate(r_stations):
    b_ref = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D_ref
    beta_deg = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)
    beta_rad = np.radians(beta_deg)
    
    # Scale and Twist
    x_scaled = x_base * b_ref
    y_scaled = y_base * b_ref
    
    x_twisted = x_scaled * np.cos(-beta_rad) - y_scaled * np.sin(-beta_rad)
    y_twisted = x_scaled * np.sin(-beta_rad) + y_scaled * np.cos(-beta_rad)
    
    # Translate
    z_station = r * (D_ref / 2.0)
    z_array = np.full_like(x_twisted, z_station)
    
    # Plotting
    ax.plot(x_twisted, y_twisted, z_array, color='b')
    
    # --- CATIA MACRO FORMATTING ---
    catia_rows.append(['StartCurve', '', ''])
    for i in range(len(x_twisted)):
        catia_rows.append([x_twisted[i], y_twisted[i], z_array[i]])
    catia_rows.append(['EndCurve', '', ''])

# Close the macro commands
catia_rows.append(['EndLoft', '', ''])
catia_rows.append(['End', '', ''])

# CSV
output_filename = "HAPS_Propeller_Macro_Input.csv"

with open(output_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(catia_rows)

print(f"Success! Data exported to {output_filename}")
print("To use: Open this CSV, copy all cells, and paste into Columns A, B, and C of 'GSD_PointSplineLoftFromExcel.xls'")

# Format plot
ax.set_xlabel('Chordwise (X) [m]')
ax.set_ylabel('Thickness (Y) [m]')
ax.set_zlabel('Spanwise (Z) [m]')
ax.set_title('3D Propeller Blade Verification (SD7037)')
ax.set_box_aspect((1, 1, 3)) 
plt.show()