import numpy as np

class Constants:
    def __init__(self):                                     # Initialise with proper values
        self.g = 9.81                                      # Gravitational acceleration (m/s^2)
        self.R = 287.05                                    # Specific gas constant for air (J/kg*K)
        self.gamma = 1.4                                   # Ratio of specific heats for air
        self.solar_constant = 1378                          # W/m^2

        self.container_inner_length = 5.867                        # m, standard ISO 20-foot container internal length
        self.container_inner_width = 2.330                         # m, standard ISO 20-foot container internal width
        self.container_inner_height = 2.350                        # m, standard ISO 20-foot container internal height
        self.container_door_width = 2.286                           # m, standard ISO 20-foot container door dimensions
        self.container_door_height = 2.261                          # m, standard ISO 20-foot container door dimensions
        self.container_mass_capacity = 28280                        # kg, standard ISO 20-foot container maximum net load mass
        self.container_packing_efficiency = 0.6                     # Assumed packing efficiency allowing for easy packing and unpacking
        self.container_packet_mass_factor = 1.1                     # Factor to account for packaging materials and safety margins in mass calculations