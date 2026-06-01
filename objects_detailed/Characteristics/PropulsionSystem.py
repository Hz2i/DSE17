from pathlib import Path

import numpy as np
import aerosandbox as asb
import ambiance as am
import matplotlib.pyplot as plt
import pandas as pd
#Get Data Points First


class PropulsionSystem:
    def __init__(self, plotdata=False, T=0.0, velocity=0.0, alt=0.0, rpm=0.0, torque=0.0, motor_temp=0.0, propeller_diameter=2.5):                                     # Initialise with proper values
        self.mass = None  # kg, estimated mass of the propulsion system
        self.alt = alt  # Altitude in meters
        self.velocity = velocity  # m/s, cruise airspeed
        self.T = T
        
        #Motor characteristics
        self.motor_temp = motor_temp
        
        # self.q = 0  # Dynamic pressure, to be calculated based on altitude and velocity
        # Propeller characteristics for a HAPS-scale low-speed propeller
        self.propeller_diameter = None  # meters
        self.propeller_area = None
        self.velocity = velocity  # m/s, cruise airspeed
        
        #Efficiencies
        self.rpm_out = rpm  # Initialize RPM output (assuming direct drive, no gearbox)
        self.motor_eff = self.calc_motor_eff(motor_temp, rpm, torque, plotdata)
        self.gearbox_eff = self.calc_gearbox_eff()
        self.propeller_eff = self.calc_propeller_eff()
        self.overall_eff = self.calc_overall_eff()

        #Power Required
        self.power_required = self.Calc_Power_Req()

    # Call AeroSandbox for propulsive calculations:
    
    def calc_motor_eff(self, motor_temp, rpm, torque, plotdata):                          # Compute all relevant characteristics of the propulsion system
        pass
    
    def calc_gearbox_eff(self):
        gearbox_eff = 0.95
        return gearbox_eff
    
    def calc_propeller_eff(self):
        pass
    
    def calc_overall_eff(self):
        if self.motor_eff is not None and self.gearbox_eff is not None and self.propeller_eff is not None:
            return self.motor_eff * self.gearbox_eff * self.propeller_eff
        else:
            return None
        
    def Calc_Power_Req(self):
        return float(self.T * self.velocity / self.overall_eff)
