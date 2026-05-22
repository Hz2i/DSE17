import numpy as np

class ClimbEnveloppe:
    def __init__(self, lat_range=None, alt=None, climbtime=None):                                 # Initialise with required values (add inputs after self, as per necessity)
        self.lat_range = lat_range  # Latitude range in degrees
        self.alt = alt  # Altitude in meters (converted from feet)
        self.climbtime = climbtime*3600  # Climb time in seconds (if provided)
    def compute_climb_enveloppe(self):
        return None
    
    def plot_climb_enveloppe(self):
        return None

if __name__ == "__main__":
    lat_range = [0,5,10,15,20,25,30]  # Latitude range in degrees
    alt = 60000*0.3048  # Altitude in meters (converted from feet)
    climbtime = 1  # Climb time in hours