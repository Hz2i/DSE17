import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Objects.Characteristics.Components_Materials import solar_panel, battery, fuel_cell


class power_storage:
    def __init__(self, power_req, latitude, days_from_solstice, DOD, batteries_used=True, energy_delta=0.0, f_c_efficiency=0.75, battery_input=battery(), fuel_cell_input=fuel_cell()):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude, days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.DOD = DOD # depth of discharge
        self.energy_delta = energy_delta # fraction of energy not recovered during daytime
        self.batteries_used = batteries_used
        self.bat = battery_input
        self.f_c = fuel_cell_input
        self.f_c_efficiency = f_c_efficiency

    def compute_weight_volume(self):        # Compute weight and volume of power storage components
        if self.batteries_used:
            energy_req = self.power_req*(86400-self.daylight_time) * (1.0+28.0*self.energy_delta) /self.DOD
            self.mass = energy_req/self.bat.massEnergy                  # Mass of Energy Storage System
            self.volume = energy_req/self.bat.volumeEnergy              # Volume of Energy Storage System
        else:
            energy_req = self.power_req*(86400-self.daylight_time) * (1.0+28.0*self.energy_delta )/self.f_c_efficiency
            self.mass = energy_req/self.f_c.massEnergy                  # Mass of Energy Storage System
            self.volume = energy_req/self.f_c.volumeEnergy              # Volume of Energy Storage System

class power_generation:
    def __init__(self, power_req, latitude, days_from_solstice, solar_panel_input=solar_panel(), energy_delta=0.0):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude,days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.avg_incidence = daylight_analysis.avg_incidence
        self.solar = solar_panel_input
        self.energy_delta = energy_delta # fraction of energy not recovered during daytime

    def compute_weight_surface(self):       # Compute weight and surface area required for power generation components
        avg_incidence = np.pi/2 - self.avg_incidence # converting to incidence angle wrt normal
        power_per_m = np.minimum(self.solar.powLimS,self.solar.efficiency*1378.0*np.cos(avg_incidence))
        self.daylight_power_req = (self.power_req*(86400-self.daylight_time) * self.energy_delta + self.power_req * self.daylight_time)/ self.daylight_time # watts
        self.area = self.daylight_power_req/power_per_m # m^2
        self.mass = self.area*self.solar.surfRho # kg

class solar_incidence:
    def __init__(self, latitude, days_from_solstice = 0):                     # Initialise known/computed values and constants
        self.axial_tilt = 23.43585/180*np.pi      # Solar Inclination Winter Solstice in Radians
        self.lat = latitude/180*np.pi             # Latitude in Radians
        self.day = days_from_solstice   # Days away from winter solstice


    def daylight_cycle(self):               # Compute incidence and daylight hours based on given parameters
        self.eq_inclination = self.axial_tilt * np.sin(2*np.pi/365*(274+self.day)) #radians
        self.max_incidence = np.pi/2 + self.eq_inclination - self.lat #max incidence angle in radians
        self.avg_incidence = self.max_incidence*2/np.pi         # radians
        self.daylight_time = 3600*2/15*180/np.pi*np.arccos(-np.tan(self.lat)*np.tan(self.eq_inclination))   # in seconds
        self.night_time = 86400-self.daylight_time                                                                              # in seconds


def power_required(self, mass=120, LD=40, prop_eff=0.8,V_cruise=25, payload=100,payload_peak=150, payload_frac=0.1,margin=300):
    output = mass*9.81/LD*V_cruise/prop_eff + payload + payload_peak*payload_frac + margin
    return output

# class power_required:
#     def __init__(self, mass=120, LD=40, prop_eff=0.8,V_cruise=25, payload=100,payload_peak=150, payload_frac=0.1,margin=300):
#         self.propulsive_power = mass*9.81/LD*V_cruise/prop_eff
#         self.payload_power = payload + payload_peak*payload_frac
#
#         self.T_int =
#         self.T_ext =
#
#
#     def thermal_power(self):

