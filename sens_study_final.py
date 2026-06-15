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

from objects_detailed.Characteristics.Components_Materials import battery, solar_panel


def size_aircraft(added_mass=0.0, added_pow=0.0, batt_en_rho=500.0, solar_eff=0.3, MTOW_in=120.0, S_in=36.0, A=20.0, TAS_initial=25.0, qc_sweep=15.0*np.pi/180.0, taper=1.0, twist=-4.675, tol=5e-3, mon_iter=5, gamma = 0.0, h_cruise = 18500.0, lat = 30.0, day_margin = 0, use_batt = True, energy_delta = 0.0, DoD = 0.8, night_time = 0.0):
    pow_frac_prev = 0.5
    payload_frac_prev = 0.1
    struct_frac_prev = 0.35
    gen_subsys_frac_prev = 0.05

    battery_mod = battery(E_m=batt_en_rho)
    solar_mod = solar_panel(eff=solar_eff*0.97**2*0.95)

    fus_geo = fuselage(D=0.0, L1=0.0, L2=0.0, L3=0.0)
    nac_geo = nacelles(nr_of_engines=0, pos=[])
    planform = airframe(S=S_in, A=A, qc_sweep=qc_sweep, taper=taper, dihedral=0.0*np.pi/180.0, twist=twist, winglet_h=2.1, fus=fus_geo, nac=nac_geo, display=False, init_polar=True)

    MTOW = MTOW_in

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
    while monitoring_var > tol or iterations < mon_iter:
        AHAPS = Aircraft(MTOW_guess=MTOW, m_skid=m_skid(), TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)

        MTOW_current = AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass + AHAPS.airframe.m_total + AHAPS.Prop_mass + AHAPS.compute_subsys_mass() + AHAPS.parasite_mass


        print("___________________________________")
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

    return AHAPS.MTOW, AHAPS.airframe.S, AHAPS.airframe.b, AHAPS.Pow_req


default_values = {
    "MTOW": 319.48,
    "Power": 3745.5569,
    "A":20.0,
    "batt_en_rho":500.0,
    "solar_eff":0.3,
    "taper":1.0,
    "twist":0.675
    }

modified_values = {
    "added_mass": 0.0,
    "added_pow": 0.0,
    "A":20.0,
    "batt_en_rho":500.0,
    "solar_eff":0.3,
    "taper":1.0,
    "twist":0.675
    }

delta = 0.05

AHAPS_ID = input("Please input the ID of the results file:")
FILE_ID = "outputs/sensitivity/" + AHAPS_ID + ".txt"
out_file = open(FILE_ID, "w")

for param in modified_values.keys():
    default = modified_values[param]

    if ~(param == "added_mass") and ~(param == "added_pow"):
        modified_values[param] = (1.0 + delta) * modified_values[param]
    elif param == "added_mass":
        modified_values[param] = (1.0) * delta * default_values["MTOW"]
    else:
        modified_values[param] = (1.0) * delta * default_values["Power"]

    MTOW_new, S_new, b_new, pow_new = size_aircraft(added_mass=modified_values["added_mass"], added_pow=modified_values["added_pow"], A=modified_values["A"], batt_en_rho=modified_values["batt_en_rho"], solar_eff=modified_values["solar_eff"], taper=modified_values["taper"], twist=modified_values["twist"])

    print("===================================", file=out_file)
    print("+0.05 factor change in:", param, file=out_file)

    print("MTOW effect:", MTOW_new - 319.48, file=out_file)
    print("Surface area effect:", S_new - 62.9389, file=out_file)
    print("Total span effect:", b_new - 35.4793, file=out_file)
    print("Total power consumption effect:", pow_new - 3745.5569, file=out_file)

    modified_values[param] = default


    if ~(param == "added_mass") and ~(param == "added_pow"):
        modified_values[param] = (1.0 - delta) * modified_values[param]
    elif param == "added_mass":
        modified_values[param] = (-1.0) * delta * default_values["MTOW"]
    else:
        modified_values[param] = (-1.0) * delta * default_values["Power"]

    MTOW_new, S_new, b_new, pow_new = size_aircraft(added_mass=modified_values["added_mass"], added_pow=modified_values["added_pow"], A=modified_values["A"], batt_en_rho=modified_values["batt_en_rho"], solar_eff=modified_values["solar_eff"], taper=modified_values["taper"], twist=modified_values["twist"])

    print("___________________________________", file=out_file)
    print("-0.05 factor change in:", param, file=out_file)

    print("MTOW effect:", MTOW_new - 319.48, file=out_file)
    print("Surface area effect:", S_new - 62.9389, file=out_file)
    print("Total span effect:", b_new - 35.4793, file=out_file)
    print("Total power consumption effect:", pow_new - 3745.5569, file=out_file)

    modified_values[param] = default
