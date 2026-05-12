import numpy as np
from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PowerSystem import power_storage, solar
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Constants import Constants


class Aircraft:
    def __init__(self, MTOW_guess=100.0, TAS=20.0, h=18000.0, gamma=0.0, wing=wing(), fus=fuselage(), emp=empennage(), comp=ComputerSystem(), comms=CommunicationSystem(), flight_con=FlightConditionsSystem(), payload=PayloadSystem(), ctrls=ControlSystem()):
        self.wing = wing
        self.fus = fus
        self.emp = emp
        self.comp = comp
        self.comms = comms
        self.flight_con = flight_con
        self.payload = payload
        self.ctrls = ctrls

        self.TAS = TAS
        self.h = h
        self.gamma = gamma
        self.T_req = MTOW_guess*Constants.g/self.wing.CL_CD + MTOW_guess*Constants.g * np.sin(np.radians(gamma))

        self.pow_store = None
        self.solar = None
        self.prop = PropulsionSystem(T=T_req, gamma=self.gamma, velocity=self.TAS, alt=self.h, rpm=1000.0, torque=4.0, motor_temp=-40.0)

        self.Pow_motor = self.prop.power_required
        self.Pow_req = self.compute_subsys_pow() + self.Pow_motor
        self.Pow_prop = MTOW_guess*Constants.g*TAS / self.wing.CL_CD

    def update_flight_conditions(self, TAS_new=20.0, h_new=18000.0, gamma_new):
        self.TAS = TAS
        self.h = h
        self.gamma = gamma
        self.T_req = MTOW_guess*Constants.g/self.wing.CL_CD + MTOW_guess*Constants.g * np.sin(np.radians(gamma))
        self.prop = PropulsionSystem(T=T_req, gamma=self.gamma, velocity=self.TAS, alt=self.h)


    def compute_subsys_pow(self):
        return self.comp.Pow + self.comms.Pow + self.flight_con.Pow + self.payload.Pow + self.ctrls.Pow


    def compute_total_mass(self):
        pass

    def compute_total_volume(self):
        pass
