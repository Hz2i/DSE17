import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PowerSystem_sizing import power_storage, power_generation, solar_incidence
from Objects.Characteristics.PropulsionSystem import prop_eff_height
from Objects.Characteristics.ReferenceGeometries import foil, fuselage, empennage
from Objects.Constants import Constants

from Objects.AircraftGeneral.Aircraft_new import Aircraft


OEM_frac = 0.3
MTOW_initial = 120.0
TAS_initial = 25.0
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


error = 1.0e5
error_vec = np.ones(5)
error_vec *= error

AHAPS = Aircraft(MTOW_guess=MTOW_initial, TAS=TAS_initial, wing=wing_geo, fus=fus_geo, emp=emp_geo)
AHAPS.compute_motor_pow(h=h_cruise)
AHAPS.compute_total_pow()

pow_storage_sys = power_storage(AHAPS.Pow_req, latitude=lat, days_from_solstice=day_margin, DOD=DoD)
pow_gen = power_generation(AHAPS.Pow_req, latitude=lat, days_from_solstice=day_margin)



# solar_conditions = solar_incidence(latitude=lat, days_from_solstice=day_margin)
