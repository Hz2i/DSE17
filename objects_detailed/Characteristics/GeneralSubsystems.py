import numpy as np

class ComputerSystem:
    def __init__(self, cable_mass_density = 0.064008, cable_diameter = 0.0007239, cable_length = 45.0):                                     # Initialise with proper values
        """required_flops: conservative estimate of the required computing power, in FLOPS
            computing_efficiency: from Raspberry Pi 5, in FLOPS/W
            mass_density: from Raspberry Pi 5, in kg/FLOPS
            volume_density: from Raspberry Pi 5, in m^3/FLOPS
            cable_mass_density: from Prysmian cat5.e SF/UTP, in kg/m
            cable_diameter: from Prysmian cat5.e SF/UTP, in m
            connected_subsystems: number of subsystems connected to the computer system
            cable_length: assumed average distance between computer system and connected subsystems, in m
            """
        #required_flops=1e9, computing_efficiency = 2 * 10**9, mass_density = 2.5 * 10**(-12), volume_density = 95.2 * 10**(-6), 
        #self.computing_power_required = required_flops      # FLOPS, order of 10^9 - 10^10 probably
        #self.computing_efficiency = computing_efficiency    # FLOPS/W, Raspberry Pi 5
        #self.mass_density = mass_density                    # kg/FLOPS, Raspberry Pi 5
        #self.volume_density = volume_density                # m^3/FLOPS, Raspberry Pi 5
        #connected_subsystems = 7
        self.redundancy = 1.0
        self.cable_mass_density = cable_mass_density        # kg/m, Prysmian cat5.e SF/UTP
        self.cable_diameter = cable_diameter                # m, Prysmian cat5.e SF/UTP
        self.comp_mass = 0.05 # raspberry pi 5
        self.comp_volume = 0.000099528
        self.comp_power = 12
        #self.connected_subsystems = connected_subsystems    # Number of subsystems connected to the computer system
        self.cable_length = cable_length                    # m, assumed average distance between computer system and connected subsystems
        self.compute_computer_system()
     
    def compute_computer_system(self):
        """Compute all relevant characteristics of the computer system, 
        by simple multiplication of the required computing power with the efficiency, mass density and volume density"""
        #self.comp_electrical_power_required = self.computing_power_required / self.computing_efficiency # W
        #self.comp_mass = self.computing_power_required * self.mass_density               
        #self.comp_volume = self.computing_power_required * self.volume_density
        self.cable_mass = self.cable_mass_density * self.cable_length #* self.connected_subsystems # kg, mass of the cables connecting the computer system to the subsystems     
        self.cable_volume = np.pi * (self.cable_diameter/2)**2 * self.cable_length #* self.connected_subsystems # m^3, volume of the cables connecting the computer system to the subsystems
        self.mass = (self.comp_mass + self.cable_mass)*self.redundancy
        self.volume = (self.comp_volume + self.cable_volume)*self.redundancy
        self.power = self.comp_power
        self.comp_x_pos = 0.0  # Needs implementation
        # return self.comp_electrical_power_required, self.comp_mass, self.comp_volume, self.comp_x_pos

class CommunicationSystem:
    def __init__(self): 
        ''', required_bitrate=2*10**6, required_distance=400*10**3, link_budget=10**(-12), transmit_power=10, 
                 transmit_efficiency=0.3, comms_mass_density=0.030, comms_volume_density=0.0005, comms_x_pos=0.0, 
                 adsb_power=1, adsb_mass=0.020, adsb_volume=21.25 * 10**(-6), 
                 elt_power=5, elt_mass=0.908, elt_volume=1.173842 * 10**(-3), 
                 ssr_power=18, ssr_mass=1.04, ssr_volume=0.01709952'''
                                            # Initialise with proper values
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
        '''
        self.bitrate_required = required_bitrate                # bps, payload needs 1 MBit/s, so 2 Mbit/s is a conservative estimate
        self.distance = required_distance                       # m, 400 km horizontally, approx. 20 km vertically
        self.link_budget = link_budget                          # 
        self.transmit_power = transmit_power                    # W
        self.transmit_efficiency = transmit_efficiency          # W/W = -
        self.comms_mass_density = comms_mass_density            # kg/W, typical for HAPS
        self.comms_volume_density = comms_volume_density        # m^3/W Catalani et al. 2020
        self.comms_x_pos = comms_x_pos 
        '''                         # Needs implementation
        self.redundancy = 1.0 # factor everything is multiplied by if redundancy is added                           # Initialise with proper values
        self.antenna_power =0.49                                     # 1621 Chelton
        self.antenna_mass=0.005579247
        self.antenna_volume =0.0005
        self.adsb_power = 20                            #W, ping 20Si
        self.adsb_mass = 0.02                             #kg, ping 20Si
        self.adsb_volume = 0.0000215                          #m^3, ping 20Si
        self.elt_power = 5                              #W Artex ELT 345
        self.elt_mass = 0.90719                                #kg Artex ELT 345
        self.elt_volume = 0.001173842                            #m^3 Artex ELT 345
        self.ssr_power = 18                              #W, Garmin GTX 345 maximum
        self.ssr_mass = 1.04                                #kg, Garmin GTX 345
        self.ssr_volume = 0.01709952                            #m^3, Garmin GTX 345

        self.compute_communication_system()

    def compute_communication_system(self):                
        """Compute all relevant characteristics of the communication system, 
        by simple summation of the characteristics of the individual components, 
        plus the contribution from the transmit power""" 
        self.mass = (self.adsb_mass + self.elt_mass + self.ssr_mass + self.antenna_mass)*self.redundancy #+ self.comms_mass_density * self.transmit_power)
        self.volume = (self.adsb_volume + self.elt_volume + self.ssr_volume +self.antenna_volume)*self.redundancy#+ self.comms_volume_density * self.transmit_power
        self.comms_x_pos = 0.0  # Needs implementation
        self.power= self.adsb_power + self.elt_power + self.ssr_power + self.antenna_power#
        #self.transmit_power * self.transmit_efficiency
        # return self.comms_mass

class FlightConditionsSystem:
    def __init__(self):          
        self.redundancy = 1.0 # factor everything is multiplied by if redundancy is added                           # Initialise with proper values
        self.twopitottubes = 2.0
        #self.mass = 0.0                                     # Currently initialised with 0; Class 2 estimation methods required!
        #self.volume = 0.0
        self.x_pos = 0.0
        self.power_required = 0.0
        self.IMU_mass = 0.006 #kg, pixhawk, includes imu, magnetometer, barometer
        self.IMU_volume = 4.63e-5 #m^3, 
        self.IMU_power_required = 7.5 #W, 
        self.GNSS_mass = 0.0068 #kg, mosaic-x5
        self.GNSS_volume = 3.84e-6 #m^3, mosaic-x5
        self.GNSS_power_required = 1.1 #W, mosaic-x5
        self.pitot_mass = 0.455 #kg, heated pitot tube
        self.pitot_volume = 0.00728 #m^3, heated pitot tube
        self.pitot_power_required = 68.6 #W, heated pitot tube heating estimate


        self.compute_flight_conditions_system()

    def compute_flight_conditions_system(self):
        self.mass = (self.IMU_mass + self.GNSS_mass + self.twopitottubes* self.pitot_mass)*self.redundancy
        self.volume = (self.IMU_volume + self.GNSS_volume + self.twopitottubes*self.pitot_volume)*self.redundancy
        self.power = self.power_required + self.IMU_power_required + self.GNSS_power_required + self.pitot_power_required # dont multiply power with redundancy, as you only send power to one
        # return self.FCS_mass, self.FCS_volume, self.FCS_power_required

class PayloadSystem: # finish cables here
    def __init__(self):                                     # Initialise with proper values
        self.volume = 0.01 #m^3 common aerospace values
        self.x_pos = 0.0
        self.mass_payload = 20.0 #kg
        self.power_required = 100.0
        self.mass_mounting = 0.05 * self.mass_payload #kg, 5% of payload mass

        '''
        self.mass_connector = 0.050 #kg, MIL-DTL 38999 Series III
        self.volume_connector = 3.12e-5
        
        self.mass_cables = 0.020 #kg/m AWG16
        self.resistance_cable = 0.013 #Ohm/m AWG16
        self.current_cable = 35 #A, AWG16 maximum current
        self.power_loss = self.resistance_cable * self.current_cable**2 #W, power loss in the cables
        self.mass_mounting = 0.05 * self.mass_payload #kg, 5% of payload mass
        #self.mass = self.compute_payload_system()
        '''

        self.compute_payload_system()

    def compute_payload_system(self):     
        self.mass = self.mass_payload+self.mass_mounting#self.mass_connector + self.mass_cables + 
        self.volume = self.volume
        self.power = self.power_required
        # return self.PS_mass, self.PS_volume, self.PS_power

class ControlSystem:
    def __init__(self):                                     # Initialise with proper values
        self.power_required = 0.0
        self.actuator_mass = 0.065 #kg, HS-7955TG servo
        self.actuator_volume = 0.040 * 0.020 * 0.037 #m^3, HS-7955TG servo
        self.actuator_current = 1.5 #A, HS-7955TG servo estimated average current
        self.actuator_power = 6 * self.actuator_current #W, 6V * 1.5A, HS-7955TG servo estimated average power
        self.mass_cables = 0.020 #kg/m AWG16
        self.diameter_cable = 0.003 #m, estimated diameter of the cables
        self.volume_cable = np.pi * (self.diameter_cable/2)**2 * 1 #m^3, estimated volume of 1 m of cable
        self.resistance_cable = 0.013 #Ohm/m AWG16
        self.power_loss_cable = self.resistance_cable * self.actuator_current**2 #W, power loss in the cables
        self.length_pushrod = 0.7 #m, estimated length of the pushrod from the actuator to the control surface
        self.diameter_pushrod = 0.005 #m, estimated diameter of the pushrod
        self.volume_pushrod = np.pi * (self.diameter_pushrod/2)**2 * self.length_pushrod #m^3, volume of the pushrod
        self.aluminum_density = 2700 #kg/m^3, density of aluminum
        self.mass_pushrod = self.volume_pushrod * self.aluminum_density #kg, mass of the pushrod
        self.joints_mass_ratio = 0.5 #Estimated mass ratio of the joints compared to the pushrod mass

        self.compute_control_system()


    def compute_control_system(self): #Assumed a total of 4 actuators (ailerons, elevator, rudder), so all characteristics are multiplied by 4
        self.mass = (self.actuator_mass + self.mass_cables + self.mass_pushrod + self.joints_mass_ratio * self.mass_pushrod) * 4
        self.volume = (self.actuator_volume + self.volume_pushrod + self.volume_cable) * 4
        self.power = (self.actuator_power + self.power_loss_cable) * 4
        # return self.CS_mass, self.CS_volume, self.CS_power_required
