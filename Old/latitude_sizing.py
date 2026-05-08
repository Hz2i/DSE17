import numpy as np
from Old.objects_prelim.Old_Power_Stuff.materials_components import battery, solar_panel, fuel_cell
from objects_prelim.aero_parameters import wing_par
from Old.objects_prelim.Old_Power_Stuff.power_sizing import power_storage, power_generation, solar_incidence, power_required


# Initialise known parameters:

MTOW_guess = 120 # kg
LD = 40 # L/D ratio
prop_eff = 0.8 # [-]
V_cruise = 25 # m/s
payload = 100 # W
payload_peak = 150 # W
payload_frac = 0.1 # [-]
margin = 300 # W

DoD = 0.8 # depth of discharge

latitude_arr = np.linspace(0.0, 60.0, 61)
days_arr = np.linspace(-27.0, 27.0, 55)

batt = battery()
sol_arr = solar_panel()
prelim_wing = wing_par(A=25, liftDrag=40)


# Perform computations:

power_req = power_required(mass, LD, prop_eff, V_cruise, payload, payload_peak, payload_frac, margin)
