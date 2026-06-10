import numpy as np
import aerosandbox as asb
import matplotlib.pyplot as plt

from objects_detailed.AircraftGeneral.Aircraft import Aircraft

from objects_detailed.Characteristics.Airframe import airframe, fuselage, nacelles
from objects_detailed.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from objects_detailed.Characteristics.ReferenceGeometries import *
from objects_detailed.Constants import Constants
from objects_detailed.Methods.LandingSkids import m_skid

# from objects_detailed.AircraftGeneral.Aircraft import Aircraft


pow_frac_prev = 0.5
payload_frac_prev = 0.1
struct_frac_prev = 0.35
gen_subsys_frac_prev = 0.05


MTOW_initial = 100.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
use_batt = True
energy_delta = 0.0
DoD = 0.8
night_time = 0.0

S = 36.0


# Flying wing planform:
fus_geo = fuselage(D=0.0, L1=0.0, L2=0.0, L3=0.0)
nac_geo = nacelles(nr_of_engines=0, pos=[])
planform = airframe(S=S, A=21.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0,fus=fus_geo, nac=nac_geo, display=False, init_polar=True)


MTOW = MTOW_initial

# Compute initial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

planform.S = AHAPS.airframe.S
MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.internal_struct.total_structure_weight + AHAPS.Prop_mass + AHAPS.compute_subsys_mass() + m_skid()

pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
payload_frac = AHAPS.payload.mass_payload / MTOW
struct_frac = AHAPS.internal_struct.total_structure_weight / MTOW
gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW

error = (abs(pow_frac - pow_frac_prev)/pow_frac_prev + abs(payload_frac - payload_frac_prev)/payload_frac_prev + abs(struct_frac - struct_frac_prev)/struct_frac_prev + abs(gen_subsys_frac - gen_subsys_frac_prev)/gen_subsys_frac_prev + abs(MTOW-MTOW_current)/MTOW)/5.0

# error = abs(MTOW-MTOW_current)/MTOW

error_vec = np.ones(5) * error
monitoring_var = np.linalg.norm(error_vec)

iterations = 0
while monitoring_var > 5e-3 or iterations < 5:
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

    MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.internal_struct.total_structure_weight + AHAPS.Prop_mass + AHAPS.compute_subsys_mass() + m_skid()


    print(f'Difference between guess and current MTOW: {MTOW_current - MTOW:.2f} kg')
    # print(f'subsystem masses: {AHAPS.compute_subsys_mass():.2f} kg')
    # print(f'internal structure mass: {AHAPS.internal_struct.total_structure_weight:.2f} kg')
    # print(f'total mass spar {AHAPS.internal_struct.total_mass_spar:.2f} kg')
    # print(f'weight skin {AHAPS.internal_struct.Weight_skin:.2f} kg')
    # print(f'power storage mass: {AHAPS.pow_store.mass:.2f} kg')
    # print(f'solar mass: {AHAPS.solar.mass:.2f} kg')
    pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW_current
    payload_frac = AHAPS.payload.mass / MTOW_current
    struct_frac = AHAPS.internal_struct.total_structure_weight / MTOW_current
    gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW_current

    error = (abs(pow_frac - pow_frac_prev)/pow_frac_prev + abs(payload_frac - payload_frac_prev)/payload_frac_prev + abs(struct_frac - struct_frac_prev)/struct_frac_prev + abs(gen_subsys_frac - gen_subsys_frac_prev)/gen_subsys_frac_prev + abs(MTOW-MTOW_current)/MTOW)/5.0

    MTOW += (MTOW_current-MTOW) * 1.0
    planform.S = AHAPS.airframe.S

    pow_frac_prev = pow_frac
    payload_frac_prev = payload_frac
    struct_frac_prev = struct_frac
    gen_subsys_frac_prev = gen_subsys_frac

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


print("Final MTOW:", MTOW_current)
print("Final power consumption:", AHAPS.Pow_req)
print("Final surface area:", AHAPS.airframe.S)
print("Final solar panel area:", AHAPS.solar.area)
print("Final energy storage system mass:", AHAPS.pow_store.mass)
print("Final energy storage system volume:", AHAPS.pow_store.volume)


K2 = AHAPS.airframe.K2
print("Oswald efficiency: ", 1/(K2 * AHAPS.airframe.AR * np.pi))
print("Max CL/CD:", AHAPS.airframe.CL_CD_max)


# save_bool = input("Save results? (Y/N)")
#
# if save_bool == "Y":
#     AHAPS_ID = input("Please input the ID of the preliminary design:")
#     FILE_ID = "outputs/final/" + AHAPS_ID + ".txt"
#     out_file = open(FILE_ID, "w")
#
#     print("GENERAL AIRCRAFT PARAMETERS:", file=out_file)
#     print(" - Final MTOW:", AHAPS.MTOW, file=out_file)
#     print(" - CL/CD at cruise:", AHAPS.CL_CD, file=out_file)
#     print(" - CD0:", AHAPS.CD0, file=out_file)
#     print(" - Final power consumption:", AHAPS.Pow_req, file=out_file)
#     print(" - Final surface area:", AHAPS.wing.S, file=out_file)
#     print(" - Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload, file=out_file)
#     print("___________________________________", file=out_file)
#     print("POWER SYSTEM PARAMETERS:", file=out_file)
#     print(" - Final solar panel area:", AHAPS.solar.area, file=out_file)
#     print(" - Final solar panel mass:", AHAPS.solar.mass, file=out_file)
#     print(" - Final solar panel power generation:", AHAPS.solar.daylight_power_req, file=out_file)
#     print(" - Final energy storage system mass:", AHAPS.pow_store.mass, file=out_file)
#     print(" - Final energy storage system volume:", AHAPS.pow_store.volume, file=out_file)
#     print(" - Final energy storage system capacity:", AHAPS.pow_store.mass*400*3600, file=out_file)
#     print("___________________________________", file=out_file)
#     print("PROPULSION SYSTEM PARAMETERS:", file=out_file)
#     print(" - Thrust required at cruise:", AHAPS.T_req, file=out_file)
#     print(" - TAS at cruise:", AHAPS.TAS_cruise, file = out_file)
#     print(" - Lambda Advance Ratio:", AHAPS.prop.lambda_adv, file=out_file)
#     print("___________________________________", file=out_file)
#     print("FUSELAGE PARAMETERS:", file=out_file)
#     print(" - Fuselage diameter:", AHAPS.fus.D, file=out_file)
#     print(" - Fuselage L1:", AHAPS.fus.L1, file=out_file)
#     print(" - Fuselage L2 (main body):", AHAPS.fus.L2, file=out_file)
#     print(" - Fuselage L3:", AHAPS.fus.L3, file=out_file)
#     print("___________________________________", file=out_file)
#     print("WING AND EMPENNAGE PARAMETERS:", file=out_file)
#     print(" - Wing surface area:", AHAPS.wing.S, file=out_file)
#     print(" - Wing aspect ratio:", AHAPS.wing.AR, file=out_file)
#     print(" - Wing full span:", AHAPS.wing.b, file=out_file)
#     print(" - Wing chord:", AHAPS.wing.root_chord, file=out_file)
#     print(" - Wing sweep:", AHAPS.wing.qc_sweep, file=out_file)
#     print(" - Wing dihedral:", AHAPS.wing.dihedral, file=out_file)
#     print(" - Wing taper:", AHAPS.wing.taper, file=out_file)
#     print(" - Sh/S:", AHAPS.Sh_S, file=out_file)
#     print(" - Sv/S:", AHAPS.Sv_S, file=out_file)
#     print(" - Vertical tail aspect ratio:", AHAPS.emp.AR_v, file=out_file)
#     print(" - Horizontal tail aspect ratio:", AHAPS.emp.AR_h, file=out_file)
#     print(" - Vertical tail Surface:", AHAPS.emp.Sv, file=out_file)
#     print(" - Horizontal tail Surface:", AHAPS.emp.Sh, file=out_file)
#     print(" - Vertical tail span:", AHAPS.emp.Sv/AHAPS.emp.root_chord_v, file=out_file)
#     print(" - Horizontal tail span:", AHAPS.emp.Sh/AHAPS.emp.root_chord_h, file=out_file)
#     print(" - Vertical tail qc sweep:", AHAPS.emp.qc_sweep_v, file=out_file)
#     print(" - Horizontal tail qc sweep:", AHAPS.emp.qc_sweep_h, file=out_file)
#     print(" - Vertical tail root chord:", AHAPS.emp.root_chord_v, file=out_file)
#     print(" - Horizontal tail root chord:", AHAPS.emp.root_chord_h, file=out_file)
#     print(" - Vertical tail taper:", AHAPS.emp.taper_v, file=out_file)
#     print(" - Horizontal tail taper:", AHAPS.emp.taper_h, file=out_file)
