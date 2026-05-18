import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage, nacelles
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants

from Objects.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.6     # From the NASA paper (mass fraction of the power system)
MTOW_initial = 120.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
DoD = 0.8
night_time = 0.0

wing_geo = wing(A=30, qc_sweep=0.0*np.pi/180, taper=1.0, dihedral=5.0*np.pi/180.0)
fus_geo = fuselage()
emp_geo = empennage()
nac_geo = nacelles(nr_of_engines=2)

# General subsystem parameters may be changed using the following (commented) code block; Sensible defaults should already be implemented

# gen_comp = ComputerSystem()
# gen_comms = CommunicationSystem()
# gen_flightCon = FlightConditionsSystem()
# gen_ctrl = ControlSystem()
# gen_payload = PayloadSystem()


MTOW = MTOW_initial

# Compute intial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo)
powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass)/MTOW
error = abs(powM_frac - powM_frac_target)/powM_frac_target


iterations = 0
while error > 1e-3:
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo)
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
print("Final surface area:", AHAPS.wing.S, AHAPS)
print("Final solar panel area:", AHAPS.solar.area)
print("Final battery mass:", AHAPS.pow_store.mass)
print("Final battery volume:", AHAPS.pow_store.volume)
print("Final remaining mass (MTOW - Power System Mass - Payload Mass):", AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload )
