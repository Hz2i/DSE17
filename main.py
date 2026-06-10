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
planform = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0, twist=-4.675, winglet_h=2.1, fus=fus_geo, nac=nac_geo, display=True, init_polar=True)

MTOW = MTOW_initial

# Compute initial error:
AHAPS = Aircraft(MTOW_guess=MTOW, m_skid=m_skid(), TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

planform.S = AHAPS.airframe.S
MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.airframe.m_total + AHAPS.Prop_mass + AHAPS.compute_subsys_mass()

pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
payload_frac = AHAPS.payload.mass_payload / MTOW
struct_frac = AHAPS.airframe.m_total / MTOW
gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW

error = (abs(pow_frac - pow_frac_prev)/pow_frac_prev + abs(payload_frac - payload_frac_prev)/payload_frac_prev + abs(struct_frac - struct_frac_prev)/struct_frac_prev + abs(gen_subsys_frac - gen_subsys_frac_prev)/gen_subsys_frac_prev + abs(MTOW-MTOW_current)/MTOW)/5.0

# error = abs(MTOW-MTOW_current)/MTOW

error_vec = np.ones(5) * error
monitoring_var = np.linalg.norm(error_vec)

iterations = 0
while monitoring_var > 5e-3 or iterations < 5:
    AHAPS = Aircraft(MTOW_guess=MTOW, m_skid=m_skid(), TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

    MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.airframe.m_total + AHAPS.Prop_mass + AHAPS.compute_subsys_mass()


    print(f'Difference between guess and current MTOW: {MTOW_current - MTOW:.2f} kg')
    # print(f'subsystem masses: {AHAPS.compute_subsys_mass():.2f} kg')
    # print(f'internal structure mass: {AHAPS.internal_struct.total_structure_weight:.2f} kg')
    # print(f'total mass spar {AHAPS.internal_struct.total_mass_spar:.2f} kg')
    # print(f'weight skin {AHAPS.internal_struct.Weight_skin:.2f} kg')
    # print(f'power storage mass: {AHAPS.pow_store.mass:.2f} kg')
    # print(f'solar mass: {AHAPS.solar.mass:.2f} kg')
    pow_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
    payload_frac = AHAPS.payload.mass_payload / MTOW
    struct_frac = (AHAPS.airframe.m_total) / MTOW
    gen_subsys_frac = AHAPS.compute_subsys_mass() / MTOW

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


save_bool = input("Save results? (Y/N)")

if save_bool == "Y":
    AHAPS_ID = input("Please input the ID of the final design:")
    FILE_ID = "outputs/final/" + AHAPS_ID

    weights_file_ID = FILE_ID + "_weights_par.txt"
    wing_file_ID = FILE_ID + "_wing_par.txt"
    powers_file_ID = FILE_ID + "_powers_par.txt"
    cruise_file_ID = FILE_ID + "_cruise_aero_par.txt"
    int_struct_file_ID = FILE_ID + "_int_struct_par.txt"

    out_file = open(weights_file_ID, "w")

    print(" - Final MTOW:", AHAPS.MTOW, file=out_file)
    print(" - Energy storage system mass:", AHAPS.pow_store.mass, file=out_file)
    print(" - Energy generation system mass:", AHAPS.solar.mass, file=out_file)
    print(" - Total structural mass:", AHAPS.airframe.m_total, file=out_file)
    print(" - Total skid mass:", AHAPS.m_skid, file=out_file)
    print(" - Total winglet mass (both included):", AHAPS.mass_winglet, file=out_file)
    print(" - Propulsion system mass:", AHAPS.Prop_mass, file=out_file)
    print(" - Total general subsystem mass:", AHAPS.compute_subsys_mass(), file=out_file)

    out_file = open(wing_file_ID, "w")

    print(" - Final wing surface area:", AHAPS.airframe.S, file=out_file)
    print(" - Solar panel area:", AHAPS.solar.area, file=out_file)
    print(" - Wing airfoil:", AHAPS.airframe.foil, file=out_file)
    print(" - Wing aspect ratio:", AHAPS.airframe.AR, file=out_file)
    print(" - Total wing span:", AHAPS.airframe.b, file=out_file)
    print(" - Wing root chord:", AHAPS.airframe.c_r, file=out_file)
    print(" - Wing tip chord:", AHAPS.airframe.c_t, file=out_file)
    print(" - Wing taper:", AHAPS.airframe.taper, file=out_file)
    print(" - Wing QC sweep:", AHAPS.airframe.qc_sweep, file=out_file)
    print(" - Wing LE sweep:", AHAPS.airframe.le_sweep, file=out_file)
    print(" - Wing dihedral:", AHAPS.airframe.dihedral, file=out_file)
    print(" - Wing twist:", AHAPS.airframe.twist, file=out_file)
    print(" - Wing internal structure parameters:", AHAPS.internal_struct.optimized_geometry, file=out_file)

    out_file = open(powers_file_ID, "w")

    print(" - Final total power consumption:", AHAPS.Pow_req, file=out_file)
    print(" - Propulsive power at cruise:", AHAPS.Pow_motor, file=out_file)
    print(" - Propeller diameter:", AHAPS.prop.D, file=out_file)
    print(" - Heating power at cruise:", AHAPS.Pow_heat, file=out_file)
    print(" - Power storage system mass:", AHAPS.pow_store.mass, file=out_file)
    print(" - Power storage system volume:", AHAPS.pow_store.volume, file=out_file)
    print(" - Power generation system mass:", AHAPS.solar.mass, file=out_file)
    print(" - Solar panel area:", AHAPS.solar.area, file=out_file)
    print(" - Total subsystem power:", AHAPS.compute_subsys_pow(), file=out_file)

    out_file = open(cruise_file_ID, "w")

    print(" - TAS at cruise:", AHAPS.TAS_cruise, file=out_file)
    print(" - Thrust at cruise:", AHAPS.T_req, file=out_file)
    print(" - CD0:", AHAPS.airframe.CD0, file=out_file)
    print(" - K1 (drag polar coefficient):", AHAPS.airframe.K1, file=out_file)
    print(" - K2 (drag polar coefficient):", AHAPS.airframe.K2, file=out_file)
    print(" - Cruise CL:", AHAPS.CL_cruise, file=out_file)
    print(" - Cruise CD:", AHAPS.CD_cruise, file=out_file)
    print(" - Lift over drag at cruise:", AHAPS.CL_CD, file=out_file)
    print(" - Cruise AoA:", AHAPS.alpha, file=out_file)
    print(" - Maximum CL:", AHAPS.airframe.CL_max, file=out_file)
    print(" - Lift gradient:", AHAPS.airframe.CL_alpha, file=out_file)
    print(" - CL at 0 AoA:", AHAPS.airframe.CL0, file=out_file)
