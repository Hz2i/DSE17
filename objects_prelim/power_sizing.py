import numpy as np
from materials_components import solar_panel, battery, fuel_cell


class power_storage:
    def __init__(self,power_req,latitude, days_from_solstice,DOD):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude,days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.DOD = DOD
        bat = battery()
        self.massEnergy = bat.massEnergy                     # Mass energy density [J/kg]
        self.volumeEnergy = bat.volumeEnergy             # Volumetric energy density [J/m^3]
        self.massRho = bat.massRho    # Mass density [kg/m^3]

    def compute_weight_volume(self):        # Compute weight and volume of power storage components
        energy_req = self.power_req*(86400-self.daylight_time)/self.DOD
        self.mass = energy_req/self.massEnergy                  # Mass of Energy Storage System
        self.volume = energy_req/self.volumeEnergy              # Volume of Energy Storage System

class power_generation:
    def __init__(self,power_req,latitude, days_from_solstice):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude,days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.avg_incidence = daylight_analysis.avg_incidence
        solar = solar_panel()
        self.efficiency = solar.efficiency      # Efficiency
        self.surfRho = solar.surfRho            # Surface mass density [kg/m^2]
        self.powLimS = solar.powLimS            # Power limit per m^2 [W/m^2]
        self.powLimM = solar.powLimM            # Power limit per kg [W/kg]

    def compute_weight_surface(self):       # Compute weight and surface area required for power generation components
        avg_incidence = np.pi/2 - self.avg_incidence
        power_per_m = np.maximum(250.0,self.efficiency*1378.0*np.cos(avg_incidence))
        self.area = self.power_req/power_per_m
        self.mass = self.area*self.surfRho

class solar_incidence:
    def __init__(self, latitude, days_from_solstice = 0):                     # Initialise known/computed values and constants
        self.axial_tilt = 23.43585/180*np.pi      # Solar Inclination Winter Solstice in Degrees
        self.lat = latitude/180*np.pi             # Latitude
        self.day = days_from_solstice   # Days away from winter solstice


    def daylight_cycle(self):               # Compute energy generated during daylight hours based on given parameters
        self.eq_inclination = self.axial_tilt * np.sin(2*np.pi/365*(274+self.day))
        self.max_incidence = np.pi/2 + self.eq_inclination - self.lat
        self.avg_incidence = self.max_incidence*2/np.pi         
        self.daylight_time = 3600*2/15*180/np.pi*np.arccos(-np.tan(self.lat)*np.tan(self.eq_inclination))   # in seconds
        self.night_time = 86400-self.daylight_time                                                                              # in seconds

power_req = 1400 # watts
latitude = 60 # deg
days_from_solstice = 0 # days
DoD = 0.8

solar_properties = solar_incidence(latitude,days_from_solstice)
solar_properties.daylight_cycle()

energy_storage = power_storage(power_req,latitude,days_from_solstice,DoD)
energy_storage.compute_weight_volume()

energy_generation = power_generation(power_req,latitude,days_from_solstice)
energy_generation.compute_weight_surface()

print(f'\n At a latitude of {latitude} degrees: \n' + 
      f'Max Solar Incidence Angle: {solar_properties.max_incidence/np.pi*180:.2f} deg' +
      f'Avg Solar Incidence Angle: {solar_properties.avg_incidence/np.pi*180:.2f} deg' +
      f'Daylight Time: {}')
