import numpy as np


class airfoil_e334: #FLYING WING
    def __init__(self):                         # Study literature and add foil parameters, if applicable
        self.clmax = 1.4
        self.max_thickness = 0.119
        self.thickness_pos = 0.303
        self.max_camber = 0.04
        self.camber_pos = 0.254


class airfoil_e387: #TRADITIONAL
    def __init__(self):                         # Study literature and add foil parameters, if applicable
        self.clmax = 1.3
        self.max_thickness = 0.091
        self.thickness_pos = 0.311
        self.max_camber = 0.032
        self.camber_pos = 0.448

class airfoil_NACA0012: #Empennage
    def __init__(self):
        self.clmax = 1.0
        self.max_thickness = 0.12
        self.thickness_pos = 0.30
        self.max_camber = 0
        self.camber_pos = 0

    

         
