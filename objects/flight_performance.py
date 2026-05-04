import numpy as np
from objects.class1_est import struct_groups                # Import weight fractions and other such preliminary values
from objects.aero_parameters import wing_par, drag_par      # Import aerodynamic parameters


class mission_profile:
    def __init__(self):                     # Initialise mission constants, according to assumptions

    def climb_phase(self):                  # Compute climb phase parameters

    def climb_reqs(self):                   # Obtain climb requirements

    def matching_diagram(self):             # Combine all requirements into diagram

    def choose_loading_thrust(self):        # Obtain lifting structure loading and thrust-to-weight | Don't forget to implement margin selection (bias towards larger or smaller areas for iteration)



class atmosphere_model:
    def __init__(self):                     # Foreseen to be a helper script; Initialise with known constants and compute atmospheric values/functions
