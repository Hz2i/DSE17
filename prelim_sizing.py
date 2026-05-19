import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage, nacelles
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants
from Objects.Performance.ScissorPlot import ScissorPlot

from Objects.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.25    # From the NASA paper (mass fraction of the power system)
MTOW_initial = 120.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
use_batt = False
energy_delta = 0.0
DoD = 0.7
night_time = 0.0

S = 36.0
Sh_S = 0.15
Sv_S = 0.1

wing_geo = wing(S=S,A=25.0, qc_sweep=0.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0)
fus_geo = fuselage()
emp_geo = empennage(S_h = S*Sh_S, S_v = S*Sv_S)
nac_geo = nacelles(nr_of_engines=4)

# General subsystem parameters may be changed using the following (commented) code block; Sensible defaults should already be implemented

# gen_comp = ComputerSystem()
# gen_comms = CommunicationSystem()
# gen_flightCon = FlightConditionsSystem()
# gen_ctrl = ControlSystem()
# gen_payload = PayloadSystem()


MTOW = MTOW_initial

# Compute initial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, Sh_S = Sh_S, Sv_S = Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt)
powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
error = abs(powM_frac - powM_frac_target)/powM_frac_target


iterations = 0
while error > 1e-3:
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD,Sh_S = Sh_S, Sv_S = Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt)
    powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
    error = abs(powM_frac - powM_frac_target)/powM_frac_target

    MTOW = (AHAPS.pow_store.mass + AHAPS.solar.mass)/powM_frac_target

    iterations += 1

    print("Iteration:", iterations)
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
print("Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload )
print("Lambda Advance Ratio:", AHAPS.prop.lambda_adv)

scissorplot = ScissorPlot(wing=AHAPS.wing,empennage=AHAPS.emp,fuselage=AHAPS.fus)
scissorplot.compute_required_coefs()
print("Sh/S sufficient: ", Sh_S > scissorplot.minimum_Sh_S(x_cg_min=0.2,x_cg_max=0.4))


scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)


AHAPS_ID = "1"
FILE_ID = "outputs/prelim_concept_" + AHAPS_ID + ".txt"
out_file = open(FILE_ID, "w")

print("Final MTOW:", AHAPS.MTOW, file=out_file)
print("Final power consumption:", AHAPS.Pow_req, file=out_file)
print("Final surface area:", AHAPS.wing.S, file=out_file)
print("Final solar panel area:", AHAPS.solar.area, file=out_file)
print("Final energy storage system mass:", AHAPS.pow_store.mass, file=out_file)
print("Final energy storage system volume:", AHAPS.pow_store.volume, file=out_file)
print("Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload, file=out_file)
print("Lambda Advance Ratio:", AHAPS.prop.lambda_adv, file=out_file)
