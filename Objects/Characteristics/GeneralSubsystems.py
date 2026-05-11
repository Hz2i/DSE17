import numpy as np

class ComputerSystem:
    def __init__(self, required_flops, computing_efficiency = 2 * 10**9, mass_density = 2.5 * 10**(-12), volume_density = 95.2 * 10**(-6)):                                     # Initialise with proper values
        self.computing_power_required = required_flops # FLOPS, order of 10^9 - 10^10 probably
        self.computing_efficiency = computing_efficiency  # FLOPS/W, Raspberry Pi 5
        self.mass_density = mass_density  # kg/FLOPS, Raspberry Pi 5
        self.volume_density = volume_density  # m^3/FLOPS, Raspberry Pi 5
     
    def compute_computer_system(self):                          # Compute all relevant characteristics of the computer system
        self.comp_electrical_power_required = self.computing_power_required / self.computing_efficiency # W
        self.comp_mass = self.computing_power_required * self.mass_density               
        self.comp_volume = self.computing_power_required * self.volume_density     
        self.comp_x_pos = 0.0  # Needs implementation
        return self.comp_electrical_power_required, self.comp_mass, self.comp_volume, self.comp_x_pos                                                   

class CommunicationSystem:
    def __init__(self, required_bitrate=2 * 10**6, required_distance=400 * 10**3, link_budget=10**(-12), transmit_power):                                     # Initialise with proper values
        self.bitrate_required = required_bitrate        # bps, payload needs 1 MBit/s, so 2 Mbit/s is a conservative estimate
        self.distance = required_distance                    # m, 400 km horizontally, approx. 20 km vertically
        self.link_budget = link_budget                 # 
        self.transmit_power = transmit_power        # W
        self.transmit_efficiency = 0.3              # W/W = -
        self.comms_mass_density = 0.030             # kg/W, typical for HAPS
        self.comms_volume_density =                 # m^3/W
        self.comms_x_pos = 0.0  # Needs implementation
        self.adsb_power = 1 #W, ping 20Si
        self.adsb_mass = 0.020 #kg, ping 20Si
        self.adsb_volume = 21.25 * 10**(-6) #m^3, ping 20Si
        self.elt_power = 5 #W Artex ELT 345
        self.elt_mass = 0.908 #kg Artex ELT 345
        self.elt_volume = 1.173842 * 10**(-3) #m^3 Artex ELT 345
        self.ssr_power = 18 #W, Garmin GTX 345 maximum
        self.ssr_mass = 1.04 #kg, Garmin GTX 345
        self.ssr_volume = 0.01709952 #m^3, Garmin GTX 345

    def compute_communication_system(self):                          # Compute all relevant characteristics of the communication system
        self.comms_mass = self.adsb_mass + self.elt_mass + self.ssr_mass + self.comms_mass_density * self.transmit_power
        self.comms_volume = self.adsb_volume + self.elt_volume + self.ssr_volume + self.comms_volume_density * self.transmit_power
        self.comms_x_pos = 0.0  # Needs implementation
        self.comms_electrical_power_required = self.adsb_power + self.elt_power + self.ssr_power + self.transmit_power * self.transmit_efficiency
        return self.comms_mass

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
        self.power_required = 0.0

    def compute_control_system(self):                          # Compute all relevant characteristics of the control system
        return None
