import numpy as np
from aero_computations import wing_par


class airframe:
    def __init__(self, wing, empennage, fuselage):
        self.wing = wing
        self.empennage = empennage
        self.fuselage = fuselage


class wing:
    def __init__(self, wing_par=wing_par(), mass=0.0, vol_int=0.0, vol_tot=0.0):
        self.wing_par=wing_par
        self.m = mass
        self.vol_int = vol_int
        self.vol_tot = vol_tot


class empennage:            #Implement similarly


class fuselage:             #Implement similarly
