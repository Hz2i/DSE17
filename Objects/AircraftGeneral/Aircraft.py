import numpy as np
import ambiance as am


from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PowerSystem_sizing import power_storage, power_generation
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Constants import Constants


class Aircraft:
    def __init__(self, MTOW_guess=100.0, TAS=20.0, h=18000.0, gamma=0.0, lat=50.0, day_margin=0, DoD=0.8, wing=wing(), fus=fuselage(), emp=empennage(), comp=ComputerSystem(), comms=CommunicationSystem(), flight_con=FlightConditionsSystem(), payload=PayloadSystem(), ctrls=ControlSystem()):
        self.MTOW = MTOW_guess
        self.const = Constants()

        self.wing = wing
        self.fus = fus
        self.emp = emp
        self.comp = comp
        self.comms = comms
        self.flight_con = flight_con
        self.payload = payload
        self.ctrls = ctrls
        self.pow_store = None
        self.solar = None
        self.prop = None

        self.TAS = TAS
        self.h = h
        self.lat = lat
        self.day_margin = day_margin
        self.DoD = DoD
        self.gamma = gamma

        self.CL_CD = None
        self.CD0 = None
        self.CL_opt = None

        self.size_wing()

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
        self.wing.S = self.wing.S / 5.0

        # Ensure current Cl and CD values are in scope:
        CL_current = None
        CD_current = None

        while surface_check:
            self.wing.compute_oswald_eff()
            self.wing.zero_lift_drag(rho_cruise=am.Atmosphere(self.h).density[0], V_cruise=self.TAS, M=0.1)
            self.fus.zero_lift_drag(rho_cruise=am.Atmosphere(self.h).density[0], V_cruise=self.TAS)
            self.emp.zero_lift_drag(rho_cruise=am.Atmosphere(self.h).density[0], V_cruise=self.TAS, M=0.1)

            self.CD0 = (self.fus.CD0 + self.wing.CD0 + self.emp.CD0)/self.wing.S
            self.CL_opt = (3.0 * self.CD0 * np.pi * self.wing.AR * self.wing.e)**0.5

            TAS_opt = (self.MTOW*self.const.g / (0.5 * am.Atmosphere(self.h).density[0] * self.wing.S * self.CL_opt))**0.5

            if TAS_opt > self.TAS:
                CL_current = self.CL_opt
                CD_current = self.CD0 + CL_current**2/(np.pi*self.wing.AR*self.wing.e)
                self.CL_CD = CL_current/CD_current
            else:
                CL_current = self.MTOW*self.const.g / (0.5 * am.Atmosphere(self.h).density[0] * self.TAS**2 * self.wing.S)
                CD_current = self.CD0 + CL_current**2/(np.pi*self.wing.AR*self.wing.e)
                self.CL_CD = CL_current/CD_current

            self.T_req = self.MTOW*self.const.g/self.CL_CD + self.MTOW*self.const.g * np.sin(np.radians(self.gamma))
            self.prop = PropulsionSystem(T=self.T_req, velocity=self.TAS, alt=self.h, rpm=1000.0, torque=4.0, motor_temp=-40.0)

            self.Pow_motor = self.prop.power_required
            self.Pow_req = self.compute_subsys_pow() + self.Pow_motor

            self.pow_store = power_storage(self.Pow_req, latitude=self.lat, days_from_solstice=self.day_margin, DOD=self.DoD)
            self.pow_store.compute_weight_volume()
            self.solar = power_generation(self.Pow_req, latitude=self.lat, days_from_solstice=self.day_margin)
            self.solar.compute_weight_surface()

            if self.solar.area < self.wing.S:
                surface_check = False
            else:
                self.wing.S += np.maximum(0.1, damping * (self.solar.area - self.wing.S))
                iterations += 1
                # print("Inner iteration:", iterations)
                # print("Surface difference:", self.solar.area - self.wing.S)

        print("Optimal CL:", self.CL_opt)
        print("CL:", CL_current)
        print("CL/CD:", self.CL_CD)
        print("Oswald efficiency:", self.wing.e)
        print("Propulsive efficiency:", self.prop.overall_eff)


    def compute_subsys_pow(self):
        return self.comp.comp_electrical_power_required + self.comms.comms_electrical_power_required + self.flight_con.FCS_power_required + self.payload.power_required + self.ctrls.power_required


    def compute_total_mass(self):
        pass

    def compute_total_volume(self):
        pass
