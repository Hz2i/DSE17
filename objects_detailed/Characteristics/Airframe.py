import numpy as np
import aerosandbox as asb
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from objects_detailed.Characteristics.ReferenceGeometries import *


class airframe:
    def __init__(self, S=36.0, A=25.0, qc_sweep=0.0, taper=1.0, dihedral=0.0 , airfoil=airfoil_e334(), fus = fuselage(), nac = nacelles()):
        self.foil = airfoil
        self.AR = A
        self.qc_sweep = qc_sweep
        self.taper = taper
        self.dihedral = dihedral
        self.S = S                
        self.CL_grad = None          # Currently initialised with None; Add method to compute
        self.CL_max = None           # Currently initialised with None; Add method to compute
        self.e = None
        self.x_ac = None

        self.m_wing = None                # Currently initialised with None; Class 2 estimation methods required!
        self.v_int_wing = None            # Currently initialised with None; Stress calculations and internal design required!
        self.v_tot_wing = None            # Currently initialised with None

        self.fuselage = fus             # Fuselage geometry and parameters
        self.nacelles = nac             # Fuselage Geometry and parameters

    def vlm_analysis(self):         # Call to AeroSandbox for VLM implementation
        pass


class fuselage:
    def __init__(self,D=0.3,L1=0.5,L2=3,L3=1.5):
        self.D = D
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3

        self.m = None           # Currently initialised with None; Class 2 estimation methods required!
        self.v_total = None     # Currently initialised with None



class nacelles:
    def __init__(self, nr_of_engines = 4, diameter = 0.5, length = 1.0, pos = [0.3, 0.6]):
        self.nr_of_engines = nr_of_engines
        self.positions = pos
        self.D = diameter
        self.L = length

        self.m = None
