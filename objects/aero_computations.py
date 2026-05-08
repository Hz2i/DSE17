import numpy as np
from objects.reference_geometries import foil, fuselage, empennage          # Import geometrical parameters


class wing_par:
    def __init__(self, A=25, qc_sweep=0.0, taper=0.0, self.dihedral=0.0 , airfoil=foil()):                  # With initial values, compute and initialise required coefficients
        self.AR = A
        self.qc_sweep = qc_sweep    # Quarter Chord sweep [rad]
        self.taper = taper
        self.dihedral = dihedral
        self.airfoil = airfoil

    def compute_CL_grad(self):

    def compute_CL_max(self):


class drag_par:
    def __init__(self):                                             # With initial values, compute and initialise required coefficients

    def compute_CD_profile(self):

    def drag_polar_coeffs(self):
