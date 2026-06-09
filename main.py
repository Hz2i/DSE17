import numpy as np
import aerosandbox as asb
import matplotlib.pyplot as plt

from objects_detailed.AircraftGeneral.Aircraft import Aircraft

from objects_detailed.Characteristics.Airframe import airframe, fuselage, nacelles
from objects_detailed.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from objects_detailed.Characteristics.ReferenceGeometries import *
from objects_detailed.Constants import Constants

# from objects_detailed.AircraftGeneral.Aircraft import Aircraft


pow_frac_prev = 0.5
payload_frac_prev = 0.1
struct_frac_prev = 0.35
gen_subsys_frac_prev = 0.05

OEM_frac = payload_frac_prev + gen_subsys_frac_prev + struct_frac_prev

MTOW_initial = 100.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
use_batt = True
energy_delta = 0.05
DoD = 0.7
night_time = 0.0

S = 36.0


# Flying wing planform:
fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
nac_geo = nacelles(nr_of_engines=4)
planform = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0,fus=fus_geo, nac=nac_geo, display=False, init_polar=True)
planform = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0,fus=fus_geo, nac=nac_geo, display=False, init_polar=True)


MTOW = MTOW_initial

# Compute initial error:
AHAPS = Aircraft(MTOW_guess=MTOW, OEM_frac=OEM_frac, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

planform.S = AHAPS.airframe.S
MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.internal_struct.total_structure_weight + AHAPS.compute_subsys_mass()

pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
payload_frac = AHAPS.payload.mass_payload / MTOW
struct_frac = AHAPS.internal_struct.total_structure_weight / MTOW
gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW

OEM_frac = payload_frac + gen_subsys_frac + struct_frac

error = (abs(pow_frac - pow_frac_prev)/pow_frac_prev + abs(payload_frac - payload_frac_prev)/payload_frac_prev + abs(struct_frac - struct_frac_prev)/struct_frac_prev + abs(gen_subsys_frac - gen_subsys_frac_prev)/gen_subsys_frac_prev + abs(MTOW-MTOW_current)/MTOW)/5.0

# error = abs(MTOW-MTOW_current)/MTOW

error_vec = np.ones(5) * error
monitoring_var = np.linalg.norm(error_vec)

iterations = 0
while monitoring_var > 5e-3 or iterations < 5:
    AHAPS = Aircraft(MTOW_guess=MTOW, OEM_frac=OEM_frac, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

    MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.internal_struct.total_structure_weight + AHAPS.Prop_mass + AHAPS.compute_subsys_mass()

    # print(f'subsystem masses: {AHAPS.compute_subsys_mass()} kg')
    # print(f'internal structure mass: {AHAPS.internal_struct.total_structure_weight} kg')
    # print(f'total mass spar {AHAPS.internal_struct.total_mass_spar} kg')
    # print(f'weight skin {AHAPS.internal_struct.Weight_skin} kg')
    # print(f'power storage mass: {AHAPS.pow_store.mass} kg')
    # print(f'solar mass: {AHAPS.solar.mass} kg')
    pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW_current
    payload_frac = AHAPS.payload.mass / MTOW_current
    struct_frac = AHAPS.internal_struct.total_structure_weight / MTOW_current
    gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW_current

    error = (abs(pow_frac - pow_frac_prev)/pow_frac_prev + abs(payload_frac - payload_frac_prev)/payload_frac_prev + abs(struct_frac - struct_frac_prev)/struct_frac_prev + abs(gen_subsys_frac - gen_subsys_frac_prev)/gen_subsys_frac_prev + abs(MTOW-MTOW_current)/MTOW)/5.0

    MTOW += (MTOW_current-MTOW) * 0.05
    planform.S = AHAPS.airframe.S

    pow_frac_prev = pow_frac
    payload_frac_prev = payload_frac
    struct_frac_prev = struct_frac
    gen_subsys_frac_prev = gen_subsys_frac

    OEM_frac = payload_frac + gen_subsys_frac + struct_frac

    iterations += 1
    error_vec = np.roll(error_vec, 1)
    error_vec[0] = error
    monitoring_var = np.linalg.norm(error_vec)

    print("Iteration:", iterations)
    print("Monitoring variable:", monitoring_var)
    print("Current error:", error)
    print('MTOW current:', MTOW_current)
    print("New MTOW estimate:", MTOW)
    print("Current power system mass fraction estimate:", pow_frac)
    print("Current structural mass fraction estimate:", struct_frac)
    print("Current payload mass fraction estimate:", payload_frac)
    print("Current general subsystem mass fraction estimate:", gen_subsys_frac)
    print("___________________________________")


print("Final MTOW:", AHAPS.MTOW)
print("Final power consumption:", AHAPS.Pow_req)
print("Final surface area:", AHAPS.airframe.S)
print("Final solar panel area:", AHAPS.solar.area)
print("Final energy storage system mass:", AHAPS.pow_store.mass)
print("Final energy storage system volume:", AHAPS.pow_store.volume)


# K2 = structure.K2
# print("Oswald efficiency: ", 1/(K2 * structure.AR * np.pi))
# print("Max CL/CD:", structure.CL_CD_max)


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

# structure.compute_load_distribution(alpha=10.0, TAS=35.0)


# print("Force vectors:", F)
# print("Vortex coordinates:", coords)

# print("Total normal force:", np.sum(F[:,2]))
# print("Span:", structure.b)


# plt.scatter(structure.vortex_coords[:,1], structure.dFx_dy_current)
# plt.show()
#
# plt.scatter(structure.vortex_coords[:,1], structure.dFy_dy_current)
# plt.show()
#
# plt.scatter(structure.vortex_coords[:,1], structure.dFz_dy_current)
# plt.show()
