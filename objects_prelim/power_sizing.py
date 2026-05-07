import numpy as np
from materials_components import solar_panel, battery, fuel_cell


class power_storage:
    def __init__(self, power_req, latitude, days_from_solstice, DOD, battery_input=battery()):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude, days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.DOD = DOD # depth of discharge
        self.bat = battery_input

    def compute_weight_volume(self):        # Compute weight and volume of power storage components
        energy_req = self.power_req*(86400-self.daylight_time)/self.DOD
        self.mass = energy_req/self.bat.massEnergy                  # Mass of Energy Storage System
        self.volume = energy_req/self.bat.volumeEnergy              # Volume of Energy Storage System

class power_generation:
    def __init__(self, power_req, latitude, days_from_solstice, solar_panel_input=solar_panel()):                     # Initialise known/computed values and constants
        daylight_analysis = solar_incidence(latitude,days_from_solstice)
        daylight_analysis.daylight_cycle()
        self.power_req = power_req
        self.daylight_time = daylight_analysis.daylight_time
        self.avg_incidence = daylight_analysis.avg_incidence
        self.solar = solar_panel_input

    def compute_weight_surface(self):       # Compute weight and surface area required for power generation components
        avg_incidence = np.pi/2 - self.avg_incidence # converting to incidence angle wrt normal
        power_per_m = np.minimum(self.solar.powLimS,self.solar.efficiency*1378.0*np.cos(avg_incidence))
        self.daylight_power_req = self.power_req * 86400 / self.daylight_time # watts
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

class power_required:
    def __init__(self, mass=120, LD=40, prop_eff=0.8,V_cruise=25, payload=100,payload_peak=150, payload_frac=0.1,margin=300):
        self.output = mass*9.81/LD*V_cruise/prop_eff + payload + payload_peak*payload_frac + margin


mass = 120 # kg
LD = 40 # L/D ratio
prop_eff = 0.8 # [-]
V_cruise = 25 # m/s
payload = 100 # W
payload_peak = 150 # W
payload_frac = 0.1 # [-]
margin = 300 # W

bat = battery()
solar_c = solar_panel()

power_req = power_required(mass,LD,prop_eff,V_cruise,payload,payload_peak,payload_frac,margin)
latitude = 50 # deg
days_from_solstice = 0 # days
DoD = 0.8 # depth of discharge

solar_properties = solar_incidence(latitude,days_from_solstice)
solar_properties.daylight_cycle()

energy_storage = power_storage(power_req.output,latitude,days_from_solstice,DoD,fuel_cell())
energy_storage.compute_weight_volume()

energy_generation = power_generation(power_req.output,latitude,days_from_solstice,solar_c)
energy_generation.compute_weight_surface()

print(f'\nAt a latitude of {latitude} degrees: \n' + 
      f'Max Solar Incidence Angle: {solar_properties.max_incidence/np.pi*180:.2f} deg \n' +
      f'Avg Solar Incidence Angle: {solar_properties.avg_incidence/np.pi*180:.2f} deg \n' +
      f'Daylight Time: {solar_properties.daylight_time/3600:.2f} hours \n' +
      f'Power Required: {power_req.output:.2f} W \n' +
      f'Battery Mass: {energy_storage.mass:.2f} kg \n' +
      f'Battery Volume: {energy_storage.volume:.3f} m^3 \n' +
      f'Solar Panel Size: {energy_generation.area:.2f} m^2 \n' +
      f'Solar Panel Mass: {energy_generation.mass:.2f} kg \n'+
      f'Remaining Mass: {mass - 20 - energy_generation.mass - energy_storage.mass:.2f} kg \n')
