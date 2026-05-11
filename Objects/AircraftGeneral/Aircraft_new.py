import numpy as np
from Objects.Characteristics.Airframe import wing, fuselage, empennage
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PowerSystem import power_storage, solar
from Objects.Characteristics.PropulsionSystem import prop_eff_height
from Objects.Constants import Constants


class Aircraft:
    def __init__(self, MTOW_guess=100.0, TAS=20.0, wing=wing(), fus=fuselage(), emp=empennage(), comp=ComputerSystem(), comms=CommunicationSystem(), flight_con=FlightConditionsSystem(), payload=PayloadSystem(), ctrls=ControlSystem()):
        self.wing = wing
        self.fus = fus
        self.emp = emp
        self.comp = comp
        self.comms = comms
        self.flight_con = flight_con
        self.payload = payload
        self.ctrls = ctrls

        self.TAS = TAS

        self.pow_store = None
        self.solar = None

        self.Pow_req = None
        self.Pow_motor = None
        self.Pow_prop = MTOW_guess*Constants.g*TAS / self.wing.CL_CD




    def compute_subsys_pow(self):
        return self.comp.Pow + self.comms.Pow + self.flight_con.Pow + self.payload.Pow + self.ctrls.Pow

    def compute_motor_pow(self, h=18000.0, T=273.15):
        self.Pow_motor = Pow_prop / prop_eff_height(h, T)

    def compute_total_pow(self):
        self.Pow_req = self.Pow_motor + self.compute_subsys_pow()
