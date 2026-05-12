import numpy as np

class ComputerSystem:
    def __init__(self):                                     # Initialise with proper values
        pass
    def compute_computer_system(self):                          # Compute all relevant characteristics of the computer system
        return None


class CommunicationSystem:
    def __init__(self):                                     # Initialise with proper values
            pass
    def compute_communication_system(self):                          # Compute all relevant characteristics of the communication system
        return None


class FlightConditionsSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = None
        self.x_pos = None
    def compute_flight_conditions_system(self):                          # Compute all relevant characteristics of the flight conditions system
        return None


class PayloadSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = 20.0 #kg
        self.volume = 0.17 #m^3
        self.power = 100 #W
        self.x_pos = None
    def compute_payload_system(self):                          # Compute all relevant characteristics of the payload system
        return None


class ControlSystem:
    def __init__(self):                                     # Initialise with proper values
        self.power = 160 #W (a little high, but based on report)
    def compute_control_system(self):                          # Compute all relevant characteristics of the control system
        return None
