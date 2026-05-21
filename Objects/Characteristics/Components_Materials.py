import numpy as np

'''
class material:
    def __init__(self):                 # Provide method to define materials with certain parameters
'''


class solar_panel:
    def __init__(self, eff=0.20, rhoS=0.665, powS=200.0, powM=1000.0):                 # Provide method to define solar panel types with certain parameters
        self.efficiency = eff           # Efficiency
        self.surfRho = rhoS             # Surface mass density [kg/m^2]
        self.powLimS = powS             # Power limit per m^2 [W/m^2]
        self.powLimM = powM             # Power limit per kg [W/kg]


class battery:
    def __init__(self, E_m=400.0, E_vol=872.0):                 # Provide method to define battery types with certain parameters
        self.massEnergy = E_m*3600.0                        # Mass energy density [J/kg]
        self.volumeEnergy = E_vol*3600.0/0.001              # Volumetric energy density [J/m^3]
        self.massRho = self.volumeEnergy/self.massEnergy    # Mass density [kg/m^3]
        self.cycle_limit_nr=400                             # Cycle limit
        self.cycle_limit_degradation=0.1                    # Cycle degradation at limit (fraction)


class fuel_cell:
    def __init__(self, E_m=750.0, E_vol=270.0):                 # Provide method to define fuel cell types with certain parameters
        self.massEnergy = E_m*3600.0                        # Mass energy density [J/kg]
        self.volumeEnergy = E_vol*3600.0/0.001              # Volumetric energy density [J/m^3]
        self.massRho = self.volumeEnergy/self.massEnergy    # Mass density [kg/m^3]
        self.cycle_limit_nr=400                             # Cycle limit
        self.cycle_limit_degradation=0.0                    # Cycle degradation at limit (fraction)
