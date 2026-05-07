import numpy as np
from objects_prelim.materials_components import battery, solar_panel, fuel_cell
from objects_prelim.aero_parameters import wing_par
from objects_prelim.power_sizing import power_storage, power_generation, solar_incidence


MTOW_guess = 90.0       # Initialise an MTOW guess
latitude_arr = np.linspace(0.0, 60.0, 61)
days_arr = np.linspace(-27.0, 27.0, 55)

batt = battery()
sol_arr = solar_panel()

prelim_wing = wing_par(A=25, liftDrag=40)
