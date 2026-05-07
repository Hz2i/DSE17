import numpy as np
from objects.reference_geometries import foil, fuselage, empennage          # Import geometrical parameters


class wing_par:
    def __init__(self, A, liftDrag):                                             # With initial values, compute and initialise required coefficients
        self.AR = A
        self.L_D = liftDrag

    def update_coefficients(self, A, liftDrag):                                  # Implement similar method that updates existing values using new guesses
        self.AR = A
        self.L_D = liftDrag


class drag_par:
    def __init__(self):                                             # With initial values, compute and initialise required coefficients

    def update_coefficients(self):                                  # Implement similar method that updates existing values using new guesses
