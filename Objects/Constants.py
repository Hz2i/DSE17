import numpy as np

class Constants:
    def __init__(self):                                     # Initialise with proper values
        self.g = 9.81                                      # Gravitational acceleration (m/s^2)
        self.R = 287.05                                    # Specific gas constant for air (J/kg*K)
        self.gamma = 1.4                                   # Ratio of specific heats for air