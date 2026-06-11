import matplotlib.pyplot as plt
import aerosandbox as asb
import aerosandbox.numpy as np
from objects_detailed.Characteristics.Airframe import *
from objects_detailed.Methods.CG_comp import *
from objects_detailed.Methods.SparGeometryParam import *

airframe = Airframe.airframe(S=55.8, A=20, qc_sweep=15.0/180*np.pi, taper=1.0, dihedral=0.0 , airfoil=asb.Airfoil("mh91"), display=False, init_polar=False)
airfoil_geometry = SparGeometryParam.AirfoilGeometry(Airframe = airframe)
plt.figure(figsize=(10, 5))
plt.plot(airfoil_geometry.xu, airfoil_geometry.yu, color='blue', label="Upper Surface")
plt.plot(airfoil_geometry.xl, airfoil_geometry.yl, color='red', label="Lower Surface")
plt.plot(airfoil_geometry.xu_off, airfoil_geometry.yu_off, color='blue', label="Upper Surface Offset", linestyle='--')
plt.plot(airfoil_geometry.xl_off, airfoil_geometry.yl_off, color='red', label="Lower Surface Offset", linestyle='--')
plt.axvline(x=airfoil_geometry.x_max_thickness, color='green', linestyle='--', label='Max Thickness Location')
plt.axvline(x=airfoil_geometry.x_min, color='green', linestyle=':', label='Available Width Boundaries')
plt.axvline(x=airfoil_geometry.x_max, color='green', linestyle=':', label='Available Width Boundaries')
plt.axhline(y=airfoil_geometry.y_upper, color='orange', linestyle='--', label='Straight Line Region')
plt.axhline(y=airfoil_geometry.y_lower, color='yellow', linestyle='--', label='Straight Line Region')
plt.scatter(airfoil_geometry.x_centroid, airfoil_geometry.y_centroid, color='black', s=100, label="Centroid Location, Conservative")
plt.fill_betweenx([airfoil_geometry.y_lower, airfoil_geometry.y_upper], airfoil_geometry.x_min, airfoil_geometry.x_max, color='gray', alpha=0.3, label='Available Region')
plt.title("Airfoil Cross-Section with Available Region")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.legend()
plt.axis('equal')
plt.grid()
plt.show()