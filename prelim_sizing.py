import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants

from Objects.AircraftGeneral.Aircraft import Aircraft


dryM_frac_target = 0.3
MTOW_initial = 120.0
TAS_initial = 20.0
gamma = 0.0
h_cruise = 18000.0
lat = 40.0
day_margin = 30
DoD = 0.8
night_time = 0.0

wing_geo = wing(A=26.5, qc_sweep=0.0, taper=1.0, dihedral=5.0*np.pi/180.0)
fus_geo = fuselage()
emp_geo = empennage()

# General subsystem parameters may be changed using the following (commented) code block; Sensible defaults should already be implemented

# gen_comp = ComputerSystem()
# gen_comms = CommunicationSystem()
# gen_flightCon = FlightConditionsSystem()
# gen_ctrl = ControlSystem()
# gen_payload = PayloadSystem()


MTOW = MTOW_initial
dM = 0.1
gradient_damping = 0.0001
iterations = 0

# Compute intial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo)

dryM_frac = (MTOW - AHAPS.pow_store.mass - AHAPS.solar.mass)/(MTOW)

error = abs(dryM_frac - dryM_frac_target)/dryM_frac_target

# # Consider implementing error monitoring code (stop iterations if convergence stagnates):
# error_vec = np.ones(5)
# error_vec *= error
# error = np.linalg.norm(error_vec)
# # Note: current implementation is non-functional: vector exists, but its variance is not monitored



# while error > 1e-3 and iterations < 1e3:
while error > 1e-3:
    # dM = 0.01*MTOW
    AHAPS = Aircraft(MTOW_guess=(MTOW+dM), TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo)

    dryM_frac = (MTOW + dM - AHAPS.pow_store.mass - AHAPS.solar.mass)/(MTOW + dM)
    error_current = abs(dryM_frac - dryM_frac_target)/dryM_frac_target

    # Note: dry mass fraction is the current optimization target
    grad = (error_current - error)/dM

    MTOW -= gradient_damping * grad

    iterations += 1
    error = error_current
    print("Iteration:", iterations)
    print("Current error:", error_current)
    print("Grad:", grad)
    print("Current MTOW estimate:", MTOW)
    print("Current dry mass fraction estimate:", dryM_frac)
    print("___________________________________")

print("Final MTOW:", MTOW)
