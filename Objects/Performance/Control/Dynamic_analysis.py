from math import *

import numpy as np
import control
import matplotlib
import aerosandbox as asb
from aerosandbox import Atmosphere

from Old.objects_prelim.Aircraft import aircraft

matplotlib.use('TkAgg')   # or 'Qt5Agg' if you have PyQt5 installed
import matplotlib.pyplot as plt

from Objects.Performance.Control.Mass_moments import Mass_moments
from Objects.Performance.Control.Control_trial import FlyingWingWithWingletsAeroBuildup

class Coeff_Values(Mass_moments, FlyingWingWithWingletsAeroBuildup):

    def __init__(self):
        mm = Mass_moments()
        cs = FlyingWingWithWingletsAeroBuildup()

        (
            CXu, CXa, CXq, CXde,
            CZu, CZa, CZq, CZde,
            Cm, Cmu, Cma, Cmq, Cmde,
            CYb, CYp, CYr, CYda, CYdr,
            Clb, Clp, Clr, Clda, Cldr,
            Cnb, Cnp, Cnr, Cnda, Cndr
        ) = cs.Coefficients()

        inputs = cs.DynamicAnalysisInputs()

        V0 = inputs["V0"]
        b = inputs["b"]
        c = inputs["c"]
        S = inputs["S"]
        rho = inputs["rho"]

        # print(self.scalar(CXu))

        super().__init__()
        self.m = mm.m_total
        self.V0 = V0
        self.b = b
        self.c = c
        self.S = S

        # Atmosphere
        self.rho0 = 1.225
        self.lmda = -0.0065
        self.Temp0 = 288.15
        self.R = 287.05
        self.g = 9.81
        self.rho = rho
        self.W = self.m * self.g

        # # Aerodynamic constants — fill in your actual values
        # self.CD0 = 0.01  # zero-lift drag coefficient
        # self.CLa = 5.1064415  # lift curve slope [1/rad]
        # self.alpha0 = 0.0  # trim angle of attack [rad]
        # self.A = self.b ** 2 / self.S  # aspect ratio (or set manually: 20)
        # self.e = 0.85  # Oswald efficiency factor
        self.th0 = 0.0  # trim pitch angle [rad]
        self.Cma = self.scalar(Cma)  # pitch moment curve slope (placeholder)

        # Constant values concerning aircraft inertia
        self.muc = self.m / (self.rho * self.S * self.c)
        self.mub = self.m / (self.rho * self.S * self.b)
        self.KX2 = self.k_x_nd ** 2
        self.KZ2 = self.k_z_nd ** 2
        self.KXZ = self.k_xz_nd
        self.KY2 = self.k_y_nd ** 2

        # Aerodynamic constants
        # self.Cmac = 0
        # self.CNwa = self.CLa
        # self.CNha = 2 * pi * self.Ah / (self.Ah + 2)
        # self.depsda = 4 / (self.A + 2)

        # # Lift and drag coefficient
        self.CL = 0.828
        self.CD = 0.0196

        # Stability derivatives
        self.CX0 = self.W * sin(self.th0) / (0.5 * self.rho * self.V0 ** 2 * self.S) # CHECK
        self.CXu = self.scalar(CXu) # todo check sign
        self.CXa = self.scalar(CXa)
        self.CXadot = 0.0
        self.CXq = self.scalar(CXq)
        self.CXde = self.scalar(CXde)
        self.CXdt = 0.0

        self.CZ0 = - self.W * cos(self.th0) / (0.5 * self.rho * self.V0 ** 2 * self.S)
        self.CZu = self.scalar(CZu) # todo check sign
        self.CZa = self.scalar(CZa)
        self.CZadot = 0.0
        self.CZq = self.scalar(CZq)
        self.CZde = self.scalar(CZde)
        self.CZdt = 0.0

        self.Cm0 = 0.0  # or self.scalar(Cm)
        self.Cmu = self.scalar(Cmu)
        self.Cmadot = 0.0
        self.Cmq = self.scalar(Cmq)
        self.CmTc = 0.0
        self.Cmde = self.scalar(Cmde)

        self.CYb = self.scalar(CYb)
        self.CYbdot = 0.0
        self.CYp = self.scalar(CYp)
        self.CYr = self.scalar(CYr)
        self.CYda = self.scalar(CYda)
        self.CYdr = self.scalar(CYdr)

        self.Clb = self.scalar(Clb)
        self.Clp = self.scalar(Clp)
        self.Clr = self.scalar(Clr)
        self.Clda = self.scalar(Clda)
        self.Cldr = self.scalar(Cldr)

        self.Cnb = self.scalar(Cnb)
        self.Cnbdot = 0.0
        self.Cnp = self.scalar(Cnp)
        self.Cnr = self.scalar(Cnr)
        self.Cnda = self.scalar(Cnda)
        self.Cndr = self.scalar(Cndr)

    def scalar(self, x):
        return float(np.asarray(x).squeeze().item())


class Dynamic_Analysis:
    def __init__(self, param=None):
        if param is None:
            param = Coeff_Values()

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

class Plotting(Dynamic_Analysis):
    def __init__(self, param=None):
        if param is None:
            param = Coeff_Values()
        super().__init__(param)
        self.param = param

    def plot_phugoid(self, u_perturb=1.0, t_end=200.0, dt=0.01):
        x0 = [u_perturb, 0.0, 0.0, 0.0]
        t = np.arange(0, t_end, dt)
        U = np.zeros((1, len(t)))
        t, y = control.forced_response(self.sys_s, T=t, U=U, X0=x0)

        # Check for divergence
        if np.any(np.abs(y) > 1e6):
            print("WARNING: phugoid response is diverging — check symmetric poles")

        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        axes[0].plot(t, y[0])
        axes[0].set_ylabel('Δu [m/s]')
        axes[1].plot(t, np.degrees(y[2]))
        axes[1].set_ylabel('θ [deg]')
        axes[2].plot(t, np.degrees(y[1]))
        axes[2].set_ylabel('α [deg]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Phugoid — Δu₀ = {u_perturb} m/s')
        plt.tight_layout()
        plt.show()

    def plot_short_period(self, alpha_perturb_deg=2.0, t_end=10.0, dt=0.001):
        x0 = [0.0, np.radians(alpha_perturb_deg), 0.0, 0.0]
        t = np.arange(0, t_end, dt)
        U = np.zeros((1, len(t)))
        print(U)
        t, y = control.forced_response(self.sys_s, T=t, U=U, X0=x0)

        fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
        axes[0].plot(t, np.degrees(y[1]))
        axes[0].set_ylabel('α [deg]')
        axes[1].plot(t, np.degrees(y[3]))
        axes[1].set_ylabel('q [deg/s]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Short-period — Δα₀ = {alpha_perturb_deg}°')
        plt.tight_layout()
        plt.show()

    def plot_dutch_roll(self, beta_perturb_deg=5.0, t_end=30.0, dt=0.01):
        x0 = [np.radians(beta_perturb_deg), 0.0, 0.0, 0.0]
        t = np.arange(0, t_end, dt)
        U = np.zeros((2, len(t)))  # 2 inputs: aileron + rudder
        t, y = control.forced_response(self.sys_a, T=t, U=U, X0=x0)  # sys_a

        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        axes[0].plot(t, np.degrees(y[0]))
        axes[0].set_ylabel('β [deg]')
        axes[1].plot(t, np.degrees(y[1]))
        axes[1].set_ylabel('φ [deg]')
        axes[2].plot(t, np.degrees(y[3]))
        axes[2].set_ylabel('r [deg/s]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Dutch roll — Δβ₀ = {beta_perturb_deg}°')
        plt.tight_layout()
        plt.show()

    def plot_spiral(self, phi_perturb_deg=5.0, dt=0.01):
        x0 = [0.0, np.radians(phi_perturb_deg), 0.0, 0.0]

        t_long = np.arange(0, 120.0, dt)
        U_long = np.zeros((2, len(t_long)))
        t_l, y_l = control.forced_response(self.sys_a, T=t_long, U=U_long, X0=x0)

        fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axes[0].plot(t_l, np.degrees(y_l[1]))
        axes[0].set_ylabel('φ [deg]')
        axes[1].plot(t_l, np.degrees(y_l[3]))
        axes[1].set_ylabel('r [deg/s]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Spiral mode — Δφ₀ = {phi_perturb_deg}°')
        plt.tight_layout()
        plt.show()

    def plot_poles(self):
        fig, ax = plt.subplots(figsize=(7, 5))
        for sys, label, color in [
            (self.sys_s, 'symmetric', 'steelblue'),
            (self.sys_a, 'asymmetric', 'tomato'),
        ]:
            poles = control.poles(sys)
            ax.scatter(poles.real, poles.imag,
                       marker='x', s=100, linewidths=2,
                       color=color, label=label)
            for p in poles:
                ax.annotate(f'  {p:.3f}', (p.real, p.imag), fontsize=7)
        ax.axvline(0, color='k', lw=0.8, ls='--')
        ax.axhline(0, color='k', lw=0.8, ls='--')
        ax.set_xlabel('Real')
        ax.set_ylabel('Imaginary')
        ax.set_title('Poles')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    def plot_aileron_response(self, aileron_deflect_deg=20.0, t_end=20.0, dt=0.001):
        x0 = [0.0, np.radians(45.0), 0.0, 0.0]
        t = np.arange(0, t_end, dt)
        U = np.zeros((2, int(t_end / dt)))
        for i in range(int(t_end / dt)):
                U[0, i] = np.radians(aileron_deflect_deg)
        print(U)
        t, y = control.forced_response(self.sys_a, T=t, U=U, X0=x0)

        fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
        axes[0].plot(t, np.degrees(y[1]))
        axes[0].set_ylabel('thi [deg]')
        axes[1].plot(t, np.degrees(y[2]))
        axes[1].set_ylabel('p [deg/s]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Aileron response — Δa = {aileron_deflect_deg}°')
        plt.tight_layout()
        plt.show()

    def plot_elevator_response(self, elevator_deflect_deg=-10.0, t_end=5.0, dt=0.001):
        x0 = [0.0, 0.0, 0.0, 0.0]
        t = np.arange(0, t_end, dt)
        U = np.zeros((1, int(t_end / dt)))
        for i in range(int(t_end / dt)):
                U[0, i] = np.radians(elevator_deflect_deg)
        print(U)
        t, y = control.forced_response(self.sys_s, T=t, U=U, X0=x0)

        fig, axes = plt.subplots(4, 1, figsize=(10, 5), sharex=True)
        axes[0].plot(t, np.degrees(y[2]))
        axes[0].set_ylabel('theta [deg]')
        axes[1].plot(t, np.degrees(y[3]))
        axes[1].set_ylabel('q [deg/s]')
        axes[2].plot(t, np.degrees(y[1]))
        axes[2].set_ylabel('alpha [deg/s]')
        axes[3].plot(t, np.degrees(y[0]))
        axes[3].set_ylabel('u [deg/s]')
        for ax in axes:
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='k', lw=0.6)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle(f'Elevator response — Δe = {elevator_deflect_deg}°')
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    p = Coeff_Values()
    plot = Plotting(p)

    print("Symmetric poles:", control.poles(plot.sys_s))
    print(f"CXu:  {p.CXu}")
    print(f"CZu:  {p.CZu}")
    print(f"CX0:  {p.CX0}")
    print(f"CZ0:  {p.CZ0}")
    print(f"CL:   {p.CL}")
    print(f"CD:   {p.CD}")
    print(f"KY2:  {p.KY2}")
    print(f"muc:  {p.muc}")
    print(f"V0:   {p.V0}")
    print(f"rho:  {p.rho}")

    plot.plot_poles()
    plot.plot_phugoid(u_perturb=1.0, t_end=200)
    plot.plot_short_period(alpha_perturb_deg=2.0, t_end=10)
    plot.plot_dutch_roll(beta_perturb_deg=5.0, t_end=30)
    plot.plot_spiral(phi_perturb_deg=5.0)
    plot.plot_aileron_response()
    plot.plot_elevator_response()
