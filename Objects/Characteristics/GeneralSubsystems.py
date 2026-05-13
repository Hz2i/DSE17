import numpy as np

class ComputerSystem:
    def __init__(self, required_flops=1e9, computing_efficiency = 2 * 10**9, mass_density = 2.5 * 10**(-12), volume_density = 95.2 * 10**(-6)):                                     # Initialise with proper values
        """required_flops: conservative estimate of the required computing power, in FLOPS
            computing_efficiency: from Raspberry Pi 5, in FLOPS/W
            mass_density: from Raspberry Pi 5, in kg/FLOPS
            volume_density: from Raspberry Pi 5, in m^3/FLOPS
            """
        self.computing_power_required = required_flops      # FLOPS, order of 10^9 - 10^10 probably
        self.computing_efficiency = computing_efficiency    # FLOPS/W, Raspberry Pi 5
        self.mass_density = mass_density                    # kg/FLOPS, Raspberry Pi 5
        self.volume_density = volume_density                # m^3/FLOPS, Raspberry Pi 5

        self.compute_computer_system()
     
    def compute_computer_system(self):
        """Compute all relevant characteristics of the computer system, 
        by simple multiplication of the required computing power with the efficiency, mass density and volume density"""
        self.comp_electrical_power_required = self.computing_power_required / self.computing_efficiency # W
        self.comp_mass = self.computing_power_required * self.mass_density               
        self.comp_volume = self.computing_power_required * self.volume_density     
        self.comp_x_pos = 0.0  # Needs implementation
        # return self.comp_electrical_power_required, self.comp_mass, self.comp_volume, self.comp_x_pos

class CommunicationSystem:
    def __init__(self, required_bitrate=2*10**6, required_distance=400*10**3, link_budget=10**(-12), transmit_power=10, 
                 transmit_efficiency=0.3, comms_mass_density=0.030, comms_volume_density=0.0005, comms_x_pos=0.0, 
                 adsb_power=1, adsb_mass=0.020, adsb_volume=21.25 * 10**(-6), 
                 elt_power=5, elt_mass=0.908, elt_volume=1.173842 * 10**(-3), 
                 ssr_power=18, ssr_mass=1.04, ssr_volume=0.01709952):                                     # Initialise with proper values
        """required_bitrate: conservative estimate of the required bitrate for the communication system, in bps
            required_distance: conservative estimate of the required communication distance, in m
            link_budget: from other tools, in W
            transmit_power: from other tools, in W
            transmit_efficiency: typical for HAPS, in W/W
            comms_mass_density: typical for HAPS, in kg/W
            comms_volume_density: from Catalani et al. 2020, in m^3/W
            comms_x_pos: needs implementation, in m
            adsb_power: from Ping 20Si, in W
            adsb_mass: from Ping 20Si, in kg
            adsb_volume: from Ping 20Si, in m^3
            elt_power: from Artex ELT 345, in W
            elt_mass: from Artex ELT 345, in kg
            elt_volume: from Artex ELT 345, in m^3
            ssr_power: from Garmin GTX 345, in W
            ssr_mass: from Garmin GTX 345, in kg
            ssr_volume: from Garmin GTX 345, in m^3
            """
        
        self.bitrate_required = required_bitrate                # bps, payload needs 1 MBit/s, so 2 Mbit/s is a conservative estimate
        self.distance = required_distance                       # m, 400 km horizontally, approx. 20 km vertically
        self.link_budget = link_budget                          # 
        self.transmit_power = transmit_power                    # W
        self.transmit_efficiency = transmit_efficiency          # W/W = -
        self.comms_mass_density = comms_mass_density            # kg/W, typical for HAPS
        self.comms_volume_density = comms_volume_density        # m^3/W Catalani et al. 2020
        self.comms_x_pos = comms_x_pos                          # Needs implementation
        self.adsb_power = adsb_power                            #W, ping 20Si
        self.adsb_mass = adsb_mass                              #kg, ping 20Si
        self.adsb_volume = adsb_volume                          #m^3, ping 20Si
        self.elt_power = elt_power                              #W Artex ELT 345
        self.elt_mass = elt_mass                                #kg Artex ELT 345
        self.elt_volume = elt_volume                            #m^3 Artex ELT 345
        self.ssr_power = ssr_power                              #W, Garmin GTX 345 maximum
        self.ssr_mass = ssr_mass                                #kg, Garmin GTX 345
        self.ssr_volume = ssr_volume                            #m^3, Garmin GTX 345

        self.compute_communication_system()

    def compute_communication_system(self):                
        """Compute all relevant characteristics of the communication system, 
        by simple summation of the characteristics of the individual components, 
        plus the contribution from the transmit power""" 
        self.comms_mass = self.adsb_mass + self.elt_mass + self.ssr_mass + self.comms_mass_density * self.transmit_power
        self.comms_volume = self.adsb_volume + self.elt_volume + self.ssr_volume + self.comms_volume_density * self.transmit_power
        self.comms_x_pos = self.comms_x_pos  # Needs implementation
        self.comms_electrical_power_required = self.adsb_power + self.elt_power + self.ssr_power + self.transmit_power * self.transmit_efficiency
        # return self.comms_mass

class FlightConditionsSystem:
    def __init__(self):                                     # Initialise with proper values
        self.mass = 0.0                                     # Currently initialised with 0; Class 2 estimation methods required!
        self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
        self.IMU_mass = 0.055 #kg, STIM320
        self.IMU_volume = 3.5 * 10**(-7) #m^3, STIM320
        self.IMU_power_required = 1.5 #W, STIM320
        self.GNSS_mass = 0.0068 #kg, mosaic-x5
        self.GNSS_volume = 0.031 * 0.031 * 0.004 #m^3, mosaic-x5
        self.GNSS_power_required = 0.6 #W, mosaic-x5
        self.pitot_mass = 0.160 #kg, heated pitot tube
        self.pitot_volume = 0.19 * 0.019 * 0.019 #m^3, heated pitot tube
        self.pitot_power_required = 30 #W, heated pitot tube heating estimate

        self.compute_flight_conditions_system()

    def compute_flight_conditions_system(self):
        self.FCS_mass = self.mass + self.IMU_mass + self.GNSS_mass + self.pitot_mass
        self.FCS_volume = self.volume + self.IMU_volume + self.GNSS_volume + self.pitot_volume
        self.FCS_power_required = self.power_required + self.IMU_power_required + self.GNSS_power_required + self.pitot_power_required
        # return self.FCS_mass, self.FCS_volume, self.FCS_power_required

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
