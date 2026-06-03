import numpy as np
import aerosandbox as asb

from objects_detailed.Characteristics.Airframe import airframe, fuselage, nacelles
from objects_detailed.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from objects_detailed.Characteristics.PropulsionSystem import PropulsionSystem
from objects_detailed.Characteristics.ReferenceGeometries import *
from objects_detailed.Constants import Constants

# from objects_detailed.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.5
payload_apprx_frac = 0.2

MTOW_initial = 120.0
TAS_initial = 25.0
gamma = 0.0
h_cruise = 18500.0
lat = 30.0
day_margin = 0
use_batt = True
energy_delta = 0.0
DoD = 0.7
night_time = 0.0

S = 36.0
Sh_S = 0 # 0.15
Sv_S = 0 # 0.07


# Flying wing planform:
fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
nac_geo = nacelles(nr_of_engines=4)
structure = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0,fus=fus_geo, nac=nac_geo, display=True)
