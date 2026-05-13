import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants

from Objects.AircraftGeneral.Aircraft import Aircraft


dryM_frac_target = 0.3
MTOW_initial = 120.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18000.0
lat = 50.0
day_margin = 0
DoD = 0.8
night_time = 0.0

wing_geo = wing(A=26.5, qc_sweep=0.0, taper=0.45, dihedral=5.0*np.pi/180.0)
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
gradient_damping = 0.05
iteration = 0

# Compute intial error:
AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo)

dryM_frac = (MTOW - AHAPS.pow_store.mass - AHAPS.solar.mass)/(MTOW)

error = ((dryM_frac - dryM_frac_target)**2/dryM_frac_target)**0.5

# Consider implementing error monitoring code (stop iterations if convergence stagnates):
error_vec = np.ones(5)
error_vec *= error
error = np.linalg.norm(error_vec)
# Note: current implementation is non-functional: vector exists, but its variance is not monitored

while error > 1e-3 and iterations < 1e3:
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD, wing=wing_geo, fus=fus_geo, emp=emp_geo)

    dryM_frac = (MTOW + dM - AHAPS.pow_store.mass - AHAPS.solar.mass)/(MTOW + dM)
    error_current = ((dryM_frac - dryM_frac_target)**2/dryM_frac_target)**0.5

    # Note: dry mass fraction is the current optimization target
    grad = (error_current - error_vec[0])/dM
    error_vec = np.roll(error_vec, 1)
    error_vec[0] = error_current
    error = np.linalg.norm(error_vec)

    MTOW -= gradient_damping * error_current

    iteration += 1
    print("Iteration:", iteration)
    print("Current error:", error_current)
    print("Current MTOW estimate:", MTOW)
    print("___________________________________")

print("Final MTOW:", MTOW)
