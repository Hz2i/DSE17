from math import *
import numpy as np
import control.matlab as control

from Objects.Performance.Control import Mass_moments


class Coeff_Values(Mass_moments):

    def __init__(self):
        # Constant values concerning atmosphere and gravity
        self.rho0 = 1.2250  # air density at sea level [kg/m^3]
        self.lmda = -0.0065  # temperature gradient in ISA [K/m]
        self.Temp0 = 288.15  # temperature at sea level in ISA [K]
        self.R = 287.05  # specific gas constant [m^2/sec^2K]
        self.g = 9.81  # gravity constant [m/sec^2]

        # Air density [kg/m^3]
        self.rho = ---
        self.W = self.m * self.g  # aircraft weight [N]

        # Constant values concerning aircraft inertia
        self.muc = self.m / (self.rho * self.S * self.c)
        self.mub = self.m / (self.rho * self.S * self.b)
        self.KX2 = 0.019
        self.KZ2 = 0.042
        self.KXZ = 0.002
        self.KY2 = 1.25 * 1.114

        # Aerodynamic constants
        # self.Cmac = 0
        # self.CNwa = self.CLa
        # self.CNha = 2 * pi * self.Ah / (self.Ah + 2)
        # self.depsda = 4 / (self.A + 2)

        # Lift and drag coefficient
        self.CL = 2 * self.W / (self.rho * self.V0 ** 2 * self.S)
        self.CD = self.CD0 + (self.CLa * self.alpha0) ** 2 / (pi * self.A * self.e)

        # Stability derivatives
        self.CX0 = self.W * sin(self.th0) / (0.5 * self.rho * self.V0 ** 2 * self.S) # CHECK
        self.CXu = -0.0
        self.CXa = +0.9868407
        self.CXadot = +0.0
        self.CXq = 1.1470932
        self.CXde = ---
        self.CXdt = ---

        self.CZ0 = -self.W * cos(self.th0) / (0.5 * self.rho * self.V0 ** 2 * self.S) # CHECK
        self.CZu = -0.0
        self.CZa = -5.1064415
        self.CZadot = -0.0
        self.CZq = -5.1677895
        self.CZde = ---
        self.CZdt = ---

        self.Cm0 = ---
        self.Cmu = +0.0
        self.Cmadot = +0.0
        self.Cmq = -5.8458962
        self.CmTc = -0.0

        self.CYb = -0.1497744
        self.CYbdot = 0.
        self.CYp = -0.1608418
        self.CYr = +0.0366229
        self.CYda = ---
        self.CYdr = ---

        self.Clb = -0.0845726
        self.Clp = -0.7594234
        self.Clr = +0.0204966
        self.Clda = -0.23088 # next todo
        self.Cldr = +0.03440

        self.Cnb = +0.1348
        self.Cnbdot = 0
        self.Cnp = -0.0602
        self.Cnr = -0.2061
        self.Cnda = -0.0120
        self.Cndr = -0.0939

class Dynamic_Analysis(param = Coeff_Values):

    def __init__(self, param):

        # Symmetric Matrices
        self.A_s = [[param.V0 * param.CXu / (2 * param.c * param.muc), param.V0 * param.CXa / (2 * param.c * param.muc),
                     param.V0 * param.CZ0 / (2 * param.c * param.muc), 0],
                    [param.V0 * param.CZu / (param.c * (2 * param.muc - param.CZadot)),
                     param.V0 * param.CZa / (param.c * (2 * param.muc - param.CZadot)),
                     -param.V0 * param.CX0 / (param.c * (2 * param.muc - param.CZadot)),
                     param.V0 * (2 * param.muc + param.CZq) / (param.c * (2 * param.muc - param.CZadot))],
                    [0, 0, 0, param.V0 / param.c],
                    [param.V0 / param.c * (param.Cmu + param.CZu * (param.Cmadot / (2 * param.muc - param.CZadot))) / (
                            2 * param.muc * param.KY2),
                     param.V0 / param.c * (param.Cma + param.CZa * (param.Cmadot / (2 * param.muc - param.CZadot))) / (
                             2 * param.muc * param.KY2),
                     -param.V0 / param.c * (param.CX0 * (param.Cmadot / (2 * param.muc - param.CZadot))) / (
                             2 * param.muc * param.KY2),
                     param.V0 / param.c * (param.Cmq + param.Cmadot * (
                             (2 * param.muc + param.CZq) / (2 * param.muc - param.CZadot))) / (
                             2 * param.muc * param.KY2)]]

        B_s = [
            [param.V0 * param.CXde / (2 * param.c * param.muc)],
            [param.V0 * param.CZde / (param.c * (2 * param.muc - param.CZadot))],
            [0],
            [param.V0 / param.c * (param.Cmde + param.CZde * param.Cmadot / (2 * param.muc - param.CZadot)) / (
                    2 * param.muc * param.KY2)]]  # Assuming cmdt = 0

        C_s = [[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]]

        D_s = np.zeros((4, 1))

        # Asymmetric Matrices
        self.A_a = [[param.V0 * param.CYb / (2 * param.b * param.mub), param.V0 * param.CL / (2 * param.b * param.mub),
                     param.V0 * param.CYp / (2 * param.b * param.mub),
                     param.V0 * (param.CYr - 4 * param.mub) / (2 * param.b * param.mub)],
                    [0, 0, 2 * param.V0 / param.b, 0],
                    [param.V0 * (param.Clb * param.KZ2 + param.Cnb * param.KXZ) / (
                            4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)), 0,
                     param.V0 * (param.Clp * param.KZ2 + param.Cnp * param.KXZ) / (
                             4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)),
                     param.V0 * (param.Clr * param.KZ2 + param.Cnr * param.KXZ) / (
                             4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2))],
                    [param.V0 * (param.Clb * param.KXZ + param.Cnb * param.KX2) / (
                            4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)), 0,
                     param.V0 * (param.Clp * param.KXZ + param.Cnp * param.KX2) / (
                             4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)),
                     param.V0 * (param.Clr * param.KXZ + param.Cnr * param.KX2) / (
                             4 * param.b * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2))]]

        B_a = [[0, (param.V0 * param.CYb) / (param.b * 2 * param.mub)],
               [0, 0],
               [(param.V0 * (param.Clda * param.KZ2 + param.Cnda * param.KXZ)) / (
                       param.b * 4 * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)),
                (param.V0 * (param.Cldr * param.KZ2 + param.Cndr * param.KXZ)) / (
                        param.b * 4 * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2))],
               [(param.V0 * (param.Clda * param.KXZ + param.Cnda * param.KX2)) / (
                       param.b * 4 * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2)),
                (param.V0 * (param.Cldr * param.KXZ + param.Cndr * param.KX2)) / (
                        param.b * 4 * param.mub * (param.KX2 * param.KZ2 - param.KXZ ** 2))]]

        C_a = [[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]]

        D_a = np.zeros((4, 2))

        self.sys_s = control.ss(self.A_s, B_s, C_s, D_s)
        self.sys_a = control.ss(self.A_a, B_a, C_a, D_a)