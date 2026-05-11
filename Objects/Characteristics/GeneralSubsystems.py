import numpy as np

class ComputerSystem:
    def __init__(self, required_flops):                                     # Initialise with proper values
        self.computing_power_required = required_flops # FLOPS, order of 10^9 - 10^10 probably
        self.computing_efficiency = 2 * 10**9  # FLOPS/W, Raspberry Pi 5
        self.mass_density = 2.5 * 10**(-12)  # kg/FLOPS, Raspberry Pi 5
        self.volume_density = 95.2 * 10**(-6)  # m^3/FLOPS, Raspberry Pi 5
     
    def compute_computer_system(self):                          # Compute all relevant characteristics of the computer system
        self.comp_electrical_power_required = self.computing_power_required / self.computing_efficiency # W
        self.comp_mass = self.computing_power_required * self.mass_density               
        self.comp_volume = self.computing_power_required * self.volume_density     
        self.comp_x_pos = 0.0  # Needs implementation
        return self.comp_electrical_power_required, self.comp_mass, self.comp_volume, self.comp_x_pos                                                   

class CommunicationSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = self.compute_communication_system()
        self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
    def compute_communication_system(self):                          # Compute all relevant characteristics of the communication system
        return None

class FlightConditionsSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = 0.0                                     # Currently initialised with 0; Class 2 estimation methods required!
        self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
    def compute_flight_conditions_system(self):                          # Compute all relevant characteristics of the flight conditions system
        return None

class PayloadSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = self.compute_payload_system()
        self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
    def compute_payload_system(self):                          # Compute all relevant characteristics of the payload system
        return None

class ControlSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = self.compute_control_system()
        self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
    def compute_control_system(self):                          # Compute all relevant characteristics of the control system
        return None
