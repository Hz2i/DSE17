import numpy as np
from Objects.Characteristics.ReferenceGeometries import foil, fuselage


class wing:
    def __init__(self, A=25, qc_sweep=0.0, taper=0.0, dihedral=0.0 , airfoil=foil()):
        self.foil = airfoil
        self.AR = A
        self.qc_sweep = qc_sweep
        self.taper = taper
        self.dihedral = dihedral
        self.S = 0.0                # Currently initialised with 0; Power sizing must be performed first to obtain it
        self.CL_grad = 0.0          # Currently initialised with 0; Add method to compute
        self.CL_max = 0.0           # Currently initialised with 0; Add method to compute
        self.CL_CD = 0.0            # Currently initialised with 0

        self.m = 0.0                # Currently initialised with 0; Class 2 estimation methods required!
        self.v_int = 0.0            # Currently initialised with 0; Stress calculations and internal design required!
        self.v_tot = 0.0            # Currently initialised with 0

    def compute_CL_grad(self):      # IMPLEMENTATION REQUIRED

    def compute_CL_max(self):       # IMPLEMENTATION REQUIRED

    def compute_oswald_eff(self):   # IMPLEMENTATION REQUIRED


class fuselage:
    def __init__(self):
        self.cd_profile = 0.0   # Compute with data from reference geometry
        self.m = 0.0            # Currently initialised with 0; Class 2 estimation methods required!
        self.v_total = 0.0      # Currently initialised with 0


class empennage:
    def __init__(self):
        self.Sh = 0.0           # Currently initialised with 0
        self.m = 0.0            # Currently initialised with 0
