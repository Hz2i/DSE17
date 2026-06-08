import numpy as np
import aerosandbox as asb
import ambiance as am
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from objects_detailed.Characteristics.Airframe import airframe
from objects_detailed.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from objects_detailed.Characteristics.PowerSystem_sizing import power_storage, power_generation
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from objects_detailed.Constants import Constants

from objects_detailed.Methods.StructuralAnalysis import bending_stress_lift, bending_stress_drag, torsional_stress
from objects_detailed.Methods.SparGeometryParam import SparGeometryOptimization, AirfoilGeometry, optimize_variables, determined_geometry


# Note for Stefan: REFACTOR CODE TO FIT NEW FLOW DIAGRAM
# It can be assumed that the general logic of this class (at least pertaining to airframe sizing) shall not significantly change.

class Aircraft:
    def __init__(self, MTOW_guess=200.0, TAS=25.0, h=18500.0, gamma=0.0, lat=30.0, day_margin=0, DoD=0.8, airframe=airframe(), comp=ComputerSystem(), comms=CommunicationSystem(), flight_con=FlightConditionsSystem(), payload=PayloadSystem(), ctrls=ControlSystem(), use_batt=True, energy_delta=0.0):
        self.MTOW = MTOW_guess
        self.const = Constants()

        self.airframe = airframe
        self.comp = comp
        self.comms = comms
        self.flight_con = flight_con
        self.payload = payload
        self.ctrls = ctrls

        self.pow_store = None
        self.use_batt = use_batt
        self.solar = None
        self.energy_delta = energy_delta
        self.prop = None

        self.TAS = TAS
        self.TAS_cruise = TAS
        self.h = h
        self.lat = lat
        self.day_margin = day_margin
        self.DoD = DoD
        self.gamma = gamma

        self.e = None
        self.CL_CD = None
        self.CD0 = None
        self.CL_opt = None

        self.size_wing()
        self.size_structure()

    def update_flight_conditions(self, TAS_new=20.0, h_new=18000.0, gamma_new=0.0):
        self.TAS = TAS_new
        self.h = h_new
        self.gamma = gamma_new
        self.T_req = self.MTOW*self.const.g/self.wing.CL_CD + self.MTOW*self.const.g * np.sin(np.radians(self.gamma))
        self.prop = PropulsionSystem(T=self.T_req, velocity=self.TAS, alt=self.h)


    def size_wing(self):
        err = 1.0e3
        iterations = 0
        damping = 0.25
        surface_check = True

        # Ensure the first iteration of area sizing undershoots:
        self.airframe.S = self.airframe.S / 5.0

        # Ensure current Cl and CD values are in scope:
        CL_current = None
        CD_current = None

        simulation_required = True
        S_simulated = self.airframe.S
        final_sim = False

        while surface_check or final_sim:
            self.airframe.define_geometry()
            if simulation_required or final_sim:
                self.airframe.compute_polar(alt=self.h, TAS=self.TAS, res=5)
                S_simulated = self.airframe.S
                simulation_required = False

            self.CL_opt = (self.airframe.K1 + (self.airframe.K1**2 + 12*self.airframe.K2*self.airframe.CD0)**0.5)/(2 * self.airframe.K2)
            TAS_opt = (self.MTOW*self.const.g / (0.5 * am.Atmosphere(self.h).density[0] * self.airframe.S * self.CL_opt))**0.5

            if TAS_opt > self.TAS:
                CL_current = self.CL_opt
                CD_current = self.airframe.CD0 + self.airframe.K1 * CL_current + self.airframe.K2 * CL_current**2
                self.CL_CD = CL_current/CD_current
                self.TAS_cruise = TAS_opt
            else:
                CL_current = self.MTOW*self.const.g / (0.5 * am.Atmosphere(self.h).density[0] * self.TAS**2 * self.airframe.S)
                CD_current = self.airframe.CD0 + self.airframe.K1 * CL_current + self.airframe.K2 * CL_current**2
                self.CL_CD = CL_current/CD_current
                self.TAS_cruise = self.TAS

            if CL_current > 0.8* self.airframe.CL_max:
                CL_current = 0.8* self.airframe.CL_max
                CD_current = self.airframe.CD0 + self.airframe.K1 * CL_current + self.airframe.K2 * CL_current**2
                self.CL_CD = CL_current/CD_current
                self.TAS_cruise = (self.MTOW*self.const.g / (0.5 * am.Atmosphere(self.h).density[0] * self.airframe.S * CL_current))**0.5

            self.alpha = (CL_current - self.airframe.CL0)/self.airframe.CL_alpha

            self.T_req = (self.MTOW*self.const.g/self.CL_CD + self.MTOW*self.const.g * np.sin(np.radians(self.gamma)))/self.airframe.nacelles.nr_of_engines
            self.prop = PropulsionSystem(T=self.T_req, velocity=self.TAS_cruise, alt=self.h, rpm=1000.0, torque=4.0, motor_temp=-40.0,propeller_diameter=1.5)

            self.Pow_motor = self.prop.power_required * self.airframe.nacelles.nr_of_engines
            self.Pow_req = self.compute_subsys_pow() + self.Pow_motor

            self.pow_store = power_storage(self.Pow_req, latitude=self.lat, days_from_solstice=self.day_margin, DOD=self.DoD, batteries_used=self.use_batt, energy_delta=self.energy_delta)
            self.pow_store.compute_weight_volume()
            self.solar = power_generation(self.Pow_req, latitude=self.lat, days_from_solstice=self.day_margin, energy_delta=self.energy_delta)
            self.solar.compute_weight_surface()

            if self.solar.area < self.airframe.S/1.025:
                surface_check = False
                final_sim = ~final_sim
            else:
                surface_check = True
                final_sim = False
                self.airframe.S += np.maximum(0.01, damping * (1.025 * self.solar.area - self.airframe.S))

                if (self.airframe.S - S_simulated)/S_simulated > 0.25:
                    simulation_required = True

                iterations += 1
                print("Inner iteration:", iterations)
                print("Wing surface:", self.airframe.S)
                print("Surface difference:", self.solar.area - self.airframe.S)

        print("Optimal CL:", self.CL_opt)
        print("CL:", CL_current)
        print("CD0:", self.airframe.CD0)
        print("CL/CD:", self.CL_CD)


    def size_structure(self):
        self.airframe.compute_load_distribution(alpha=self.alpha, TAS=self.TAS_cruise, alt=self.h, res=20)

        I_lift_spar, I_lift_connection = bending_stress_lift(airframe=self.airframe, ult_safety_factor=1.5)
        I_drag_spar, I_drag_connection = bending_stress_drag(airframe=self.airframe, ult_safety_factor=1.5)
        t_skin = torsional_stress(airframe=self.airframe, ult_safety_factor=1.5)

        print("I_xx_spar_req:", I_lift_spar)
        print("I_yy_spar_req:", I_drag_spar)
        print("I_xx_sleeve_req:", I_lift_connection)
        print("I_yy_sleeve_req:", I_drag_connection)

        airfoil_geometry = AirfoilGeometry(self.airframe, t_skin_airfoil=t_skin, Available_width=0.1, plot=False)

        self.internal_struct = SparGeometryOptimization(
            I_xx_spar_req=I_lift_spar,
            I_yy_spar_req=I_drag_spar,
            I_xx_sleeve_req=I_lift_connection,
            I_yy_sleeve_req=I_drag_connection,
            n_sections=6,
            min_eccentricity_factor=1.1,
            airframe=self.airframe,
            airfoil_geometry=airfoil_geometry,
            Plot=False
        )


    def compute_subsys_pow(self):
        return self.comp.power + self.comms.power + self.flight_con.power+ self.payload.power + self.ctrls.power

    def compute_subsys_mass(self):
        return self.comp.mass + self.comms.mass + self.flight_con.mass + self.ctrls.mass
