import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage, nacelles
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants
from Objects.Performance.ScissorPlot import ScissorPlot

from Objects.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.5    # From the NASA paper (mass fraction of the power system): 0.30 for the fuel cells, 0.5 or 0.6 for batteries
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
Sh_S = 0.15
Sv_S = 0.1

# Choose planform type (uncomment the required one):

# Traditional wing planform:
# wing_geo = wing(S=S,A=25.0, qc_sweep=0.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0)
# fus_geo = fuselage()
# emp_geo = empennage(S_h = S*Sh_S, S_v = S*Sv_S)
# nac_geo = nacelles(nr_of_engines=4)

# Flying wing planform:
wing_geo = wing(S=S,A=25.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0, airfoil=airfoil_e334())
fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
emp_geo = empennage(S_h = 0.0, S_v = 0.0)
nac_geo = nacelles(nr_of_engines=4)

# General subsystem parameters may be changed using the following (commented) code block; Sensible defaults should already be implemented

# gen_comp = ComputerSystem()
# gen_comms = CommunicationSystem()
# gen_flightCon = FlightConditionsSystem()
# gen_ctrl = ControlSystem()
# gen_payload = PayloadSystem()


MTOW = MTOW_initial

# Compute initial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, Sh_S = Sh_S, Sv_S = Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt, energy_delta=energy_delta)

powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
payload_frac = AHAPS.payload.mass_payload / MTOW
error = abs(powM_frac - powM_frac_target)/powM_frac_target + abs(payload_frac - payload_apprx_frac)/payload_apprx_frac

# powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
# error = abs(powM_frac - powM_frac_target)/powM_frac_target

error_vec = np.ones(10) * error
monitoring_var = np.linalg.norm(abs(error_vec - np.mean(error_vec)))

iterations = 0
while monitoring_var > 5e-3 or iterations < 10:
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD,Sh_S = Sh_S, Sv_S = Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt, energy_delta=energy_delta)

    powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
    payload_frac = AHAPS.payload.mass_payload / MTOW
    error = abs(powM_frac - powM_frac_target)/powM_frac_target + abs(payload_frac - payload_apprx_frac)/payload_apprx_frac

    # powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
    # error = abs(powM_frac - powM_frac_target)/powM_frac_target

    MTOW = (AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass_payload)/(powM_frac_target + payload_apprx_frac)

    # MTOW = (AHAPS.pow_store.mass + AHAPS.solar.mass)/(powM_frac_target)

    iterations += 1
    error_vec = np.roll(error_vec, 1)
    error_vec[0] = error
    monitoring_var = np.linalg.norm(abs(error_vec - np.mean(error_vec)))

    print("Iteration:", iterations)
    print("Monitoring variable:", monitoring_var)
    print("Current error:", error)
    print("Current MTOW estimate:", MTOW)
    print("Current power system mass fraction estimate:", powM_frac)
    print("___________________________________")


print("Final MTOW:", AHAPS.MTOW)
print("Final power consumption:", AHAPS.Pow_req)
print("Final surface area:", AHAPS.wing.S)
print("Final solar panel area:", AHAPS.solar.area)
print("Final energy storage system mass:", AHAPS.pow_store.mass)
print("Final energy storage system volume:", AHAPS.pow_store.volume)
print("Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload)
print("Lambda Advance Ratio:", AHAPS.prop.lambda_adv)


save_bool = input("Save results? (Y/N)")

if save_bool == "Y":
    AHAPS_ID = input("Please input the ID of the preliminary design:")
    FILE_ID = "outputs/prelim_concept_" + AHAPS_ID + ".txt"
    out_file = open(FILE_ID, "w")

    print("GENERAL AIRCRAFT PARAMETERS:", file=out_file)
    print(" - Final MTOW:", AHAPS.MTOW, file=out_file)
    print(" - CL/CD at cruise:", AHAPS.CL_CD, file=out_file)
    print(" - CD0:", AHAPS.CD0, file=out_file)
    print(" - Final power consumption:", AHAPS.Pow_req, file=out_file)
    print(" - Final surface area:", AHAPS.wing.S, file=out_file)
    print(" - Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload, file=out_file)
    print("___________________________________", file=out_file)
    print("POWER SYSTEM PARAMETERS:", file=out_file)
    print(" - Final solar panel area:", AHAPS.solar.area, file=out_file)
    print(" - Final solar panel mass:", AHAPS.solar.mass, file=out_file)
    print(" - Final energy storage system mass:", AHAPS.pow_store.mass, file=out_file)
    print(" - Final energy storage system volume:", AHAPS.pow_store.volume, file=out_file)
    print("___________________________________", file=out_file)
    print("PROPULSION SYSTEM PARAMETERS:", file=out_file)
    print(" - Thrust required at cruise:", AHAPS.T_req, file=out_file)
    print(" - TAS at cruise:", AHAPS.TAS_cruise, file = out_file)
    print(" - Lambda Advance Ratio:", AHAPS.prop.lambda_adv, file=out_file)
    print("___________________________________", file=out_file)
    print("WING AND EMPENNAGE PARAMETERS:", file=out_file)
    print(" - Wing surface area:", AHAPS.wing.S, file=out_file)
    print(" - Wing aspect ratio:", AHAPS.wing.AR, file=out_file)
    print(" - Wing full span:", AHAPS.wing.b, file=out_file)
    print(" - Wing chord:", AHAPS.wing.root_chord, file=out_file)
    print(" - Wing sweep:", AHAPS.wing.qc_sweep, file=out_file)
    print(" - Wing dihedral:", AHAPS.wing.dihedral, file=out_file)
    print(" - Wing taper:", AHAPS.wing.taper, file=out_file)
    print(" - Sh/S:", AHAPS.Sh_S, file=out_file)
    print(" - Sv/S:", AHAPS.Sv_S, file=out_file)
    print(" - Vertical tail aspect ratio:", AHAPS.emp.AR_v, file=out_file)
    print(" - Horizontal tail aspect ratio:", AHAPS.emp.AR_h, file=out_file)
    print(" - Vertical tail qc sweep:", AHAPS.emp.qc_sweep_v, file=out_file)
    print(" - Horizontal tail qc sweep:", AHAPS.emp.qc_sweep_h, file=out_file)
    print(" - Vertical tail root chord:", AHAPS.emp.root_chord_v, file=out_file)
    print(" - Horizontal tail root chord:", AHAPS.emp.root_chord_h, file=out_file)
    print(" - Vertical tail taper:", AHAPS.emp.taper_v, file=out_file)
    print(" - Horizontal tail taper:", AHAPS.emp.taper_h, file=out_file)




scissorplot = ScissorPlot(wing=AHAPS.wing,empennage=AHAPS.emp,fuselage=AHAPS.fus)
scissorplot.compute_required_coefs()
print("Sh/S sufficient: ", Sh_S > scissorplot.minimum_Sh_S(x_cg_min=0.2,x_cg_max=0.4))


scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)
