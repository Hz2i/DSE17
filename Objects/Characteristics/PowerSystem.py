import numpy as np

class solar_incidence:
    def __init__(self, latitude=50, days_from_solstice = 0):                     # Initialise known/computed values and constants
        self.axial_tilt = 23.43585/180*np.pi      # Solar Inclination Winter Solstice in Radians
        self.lat = latitude/180*np.pi             # Latitude in Radians
        self.day = days_from_solstice   # Days away from winter solstice


    def daylight_cycle(self):               # Compute incidence and daylight hours based on given parameters
        self.eq_inclination = self.axial_tilt * np.sin(2*np.pi/365*(274+self.day)) #radians
        self.max_incidence = np.pi/2 + self.eq_inclination - self.lat #max incidence angle in radians
        self.avg_incidence = self.max_incidence*2/np.pi         # radians
        self.daylight_time = 3600*2/15*180/np.pi*np.arccos(-np.tan(self.lat)*np.tan(self.eq_inclination))   # in seconds
        self.night_time = 86400-self.daylight_time

class solar:
    def __init__(self,lat,A,x=None,days_from_solstice=0,powLimS=250,solar_efficiency=0.31,rhoS=0.25):
        self.solar_constant = 1378 # W/m^2
        self.powLimS = powLimS # W/m^2
        self.solar_efficiency = solar_efficiency # [-]
        self.rhoS = rhoS # kg/m^2

        self.daylight_analysis = solar_incidence(lat,days_from_solstice)
        self.daylight_analysis.daylight_cycle()

        ## Outputs ##
        self.A = A # m^2
        self.M = self.A*self.rhoS # kg
        self.x = x # m, CG position, to be decided!!! (based on wing position)
        self.Pow = self.P(A,h,lat,days_from_solstice) # W, overall average power output during daylight hours


    def P(A,h,lat,days_from_solstice):
        self.avg_incidence = np.pi/2 - self.daylight_analysis.avg_incidence # converting to incidence angle wrt normal
        power = A*np.minimum(self.powLimS,self.solar_efficiency*self.solar_constant*np.cos(self.avg_incidence))

        # How does efficiency change with height???

        return power


class power_storage:
    def __init__(self, E, x=None, nr_of_cycles=0,DoD=80, E_m=450, E_vol=1150):
        self.massEnergy = E_m*3600.0                        # Mass energy density [J/kg]
        self.volumeEnergy = E_vol*3600.0/0.001              # Volumetric energy density [J/m^3]
        self.massRho = self.volumeEnergy/self.massEnergy    # Mass density [kg/m^3]

        self.DoD = DoD                                      # Depth of Discharge [-]
        self.E = E                                          # Design Energy Capacity [J]
        self.V = self.E/self.volumeEnergy                   # Volume [m^3]
        self.M = self.E/self.massEnergy                     # Mass [kg]
        self.x = x                                          # m, CG position, to be decided!!! (based on wing position)

        def reduced_capacity(E,nr_of_cycles):
            # SA03 battery from Amprius
            # 100% DOD to 80% SOH: 300 cycles
            # 70% DOD to 90% SOH: 700 cycles
            # Assumption: discharge between -20 and 55 deg C; charge between 10 and 55 deg C
            return E * (1-nr_of_cycles/400*0.2)

        self.current_capacity = reduced_capacity(self.E,nr_of_cycles) # Current Energy Capacity [J], after nr_of_cycles
