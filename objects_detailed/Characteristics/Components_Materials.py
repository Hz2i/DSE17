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
        self.cycle_limit_degradation=0.2                    # Cycle degradation at limit (fraction)


class fuel_cell:
    def __init__(self, E_m=750.0, E_vol=270.0):                 # Provide method to define fuel cell types with certain parameters
        self.massEnergy = E_m*3600.0                        # Mass energy density [J/kg]
        self.volumeEnergy = E_vol*3600.0/0.001              # Volumetric energy density [J/m^3]
        self.massRho = self.volumeEnergy/self.massEnergy    # Mass density [kg/m^3]

class CFRP:
    def __init__(self):                 # Provide method to define CFRP with certain parameters
        self.E = 230e9                     # Young's modulus [Pa]
        self.sigma = 3400e6                  
        self.rho = 1720   #https://www.researchgate.net/figure/Properties-of-carbon-fiber-reinforced-polymer-CFRP_tbl3_342938721               # Density [kg/m^3]
        self.ply_thickness = 0.000335
        self.poisson = 0.297#https://www.sciencedirect.com/science/article/pii/S0142941824000321
        self.G = self.E/(2*(1+self.poisson)) # Shear modulus [Pa]

class PET:
    def __init__(self, E=2.8e9, rho=1380.0):                 # Provide method to define PET with certain parameters
        self.E = 2.95e9                      # Young's modulus [Pa]
        self.shear = 55e6 #max shear stress PEThttps://www.ensinger-pc.com/injection-molding-materials/our-plastic-stock-shapes/pet-thermoplastic-polyester/
        self.G = 1.375e9#https://designerdata.nl/materials/plastics/thermo-plastics/polyethylene-terephthalate
        self.rho = 1200

class GLARE:
    def __init__(self):                 # Provide method to define GLARE with certain parameters
        self.sigma=403.66e6#https://www.academia.edu/104528943/ON_FABRICATION_AND_TESTING_OF_GLARE
        self.rho=2400#https://www.researchgate.net/figure/Properties-of-Glare-components_tbl1_343586328