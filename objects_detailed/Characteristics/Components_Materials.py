import numpy as np

'''
class material:
    def __init__(self):                 # Provide method to define materials with certain parameters
'''


class solar_panel:
    def __init__(self, eff=0.30*0.97**2*0.95, rhoS=0.400, powS=270.0, powM=500.0):                 # Provide method to define solar panel types with certain parameters
        self.efficiency = eff           # Efficiency
        self.surfRho = rhoS             # Surface mass density [kg/m^2]
        self.powLimS = powS             # Power limit per m^2 [W/m^2]
        self.powLimM = powM             # Power limit per kg [W/kg]


class battery:
    def __init__(self, E_m=500.0*0.961, E_vol=872.0*0.961):                 # Provide method to define battery types with certain parameters
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
        self.min_thickness = 0.0005

class Mylar:#change to mylar
    def __init__(self, E=2.8e9, rho=1380.0):                 # Provide method to define PET with certain parameters
        self.E = 9.80665*510e6                      # Young's modulus [Pa]
        self.shear = 9.80665*15e6 #max shear stress PET
        self.poisson = 0.38
        self.G = self.E/(2*(1+self.poisson)) # Shear modulus [Pa]
        self.rho = 1390#https://stenbacka.fi/wp-content/uploads/sites/3/2016/07/mylar_a_fysikaaliset_ominaisuudet.pdf
        self.min_thickness = 0.0001 #Still Determine

class Titanium:
    def __init__(self):                 # Provide method to define GLARE with certain parameters
        self.sigma=880.0e6#https://www.academia.edu/104528943/ON_FABRICATION_AND_TESTING_OF_GLARE
        self.rho=4430#https://www.researchgate.net/figure/Properties-of-Glare-components_tbl1_343586328
        self.min_thickness=0.0005
        self.E = 55.5e9
        self.shear = 181.4e6 # https://scispace.com/pdf/plasticity-correction-factors-for-buckling-of-flat-4y7cnffw7y.pdf
        #mylar + ressst subsystems in code

class Silicone_Rubber:
    def __init__(self):
        self.E = 0.05e9                     # Young's modulus [Pa]
        self.poisson = 0.33
        self.rho = 1100  # Density [kg/m^3]

class Aluminum7075: #https://www.jeelix.com/7075-t6-aluminum-material-properties/
    def __init__(self):
        self.rho = 2810   # Density [kg/m^3]

class PA6: #https://designerdata.nl/materials/plastics/thermo-plastics/polyamide-6
    def __init__(self):
        self.rho = 1150   # Density [kg/m^3]
        self.E =4.6e9
        self.sigma = 66.5e6
