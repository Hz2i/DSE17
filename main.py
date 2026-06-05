import numpy as np
import aerosandbox as asb
import matplotlib.pyplot as plt

from objects_detailed.Characteristics.Airframe import airframe, fuselage, nacelles
from objects_detailed.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from objects_detailed.Characteristics.PropulsionSystem import PropulsionSystem
from objects_detailed.Characteristics.ReferenceGeometries import *
from objects_detailed.Constants import Constants

# from objects_detailed.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.5
payload_apprx_frac = 0.2

MTOW_initial = 120.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
use_batt = True
energy_delta = 0.0
DoD = 0.7
night_time = 0.0

S = 36.0


# Flying wing planform:
fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
nac_geo = nacelles(nr_of_engines=4)
structure = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0,fus=fus_geo, nac=nac_geo, display=True)


K2 = structure.K2
print("Oswald efficiency: ", 1/(K2 * structure.AR * np.pi))
print("Max CL/CD:", structure.CL_CD_max)


# aero = structure.llt_analysis(series=True, alpha=np.linspace(-10.0, 20.0, 30))
#
# aero_single, llt_an = structure.llt_analysis(alpha=10.0)
# print("Panel coords:", llt_an.front_left_vertices)
#
# lift = aero["CL"]
# drag = aero["CD"]
#
# plt.plot(np.linspace(-10.0, 20.0, 30), lift)
# plt.show()
# plt.plot(lift, drag)
# plt.plot(lift, structure.CD0 + structure.K1 * lift + structure.K2 * lift**2, c='r')
# plt.show()


F, M, coords = structure.compute_force_distribution(alpha=10)

print("Force vectors:", F)
print("Vortex coordinates:", coords)

print("Span:", structure.b)


plt.scatter(coords[:,1], F[:,0])
plt.show()

plt.scatter(coords[:,1], F[:,1])
plt.show()

plt.scatter(coords[:,1], F[:,2])
plt.show()
