import numpy as np

from Objects.Characteristics.Airframe import wing, fuselage, empennage, nacelles
from Objects.Characteristics.GeneralSubsystems import ComputerSystem, CommunicationSystem, FlightConditionsSystem, PayloadSystem, ControlSystem
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import *
from Objects.Constants import Constants
from Objects.Performance.ScissorPlot import ScissorPlot

from Objects.AircraftGeneral.Aircraft import Aircraft


powM_frac_target = 0.6  # Target power system mass fraction
payload_apprx_frac = 0.2  # Approximate payload mass fraction

MTOW_initial = 120.0  # Initial guess for MTOW
TAS_initial = 25.0  # True airspeed
gamma = 0.0  # Flight path angle
h_cruise = 18500.0  # Cruise altitude
lat = 30.0  # Latitude
day_margin = 0  # Daylight margin
use_batt = True  # Use battery
energy_delta = 0.0  # Energy delta
DoD = 0.7  # Depth of discharge
night_time = 0.0  # Night time duration

S = 36.0  # Wing surface area
Sh_S = 0.15  # Horizontal tail area ratio
Sv_S = 0.1  # Vertical tail area ratio


def compute(d_AR=0, d_sweep=0, d_taper=0, d_thickness=0, d_x_thickness=0):
    
    class airfoil_e387: #TRADITIONAL
        def __init__(self):                         # Study literature and add foil parameters, if applicable
            self.clmax = 1.3
            self.max_thickness = 0.091+d_thickness
            self.thickness_pos = 0.311+d_x_thickness
            self.max_camber = 0.032
            self.camber_pos = 0.448
            self.cm_0 = -0.09

    # Traditional wing planform
    wing_geo = wing(S=S, A=25.0 + d_AR, qc_sweep=(0.0 +d_sweep)* np.pi / 180, taper=(1.0+d_taper), dihedral=0.0 * np.pi / 180.0, airfoil=airfoil_e387())
    fus_geo = fuselage()
    emp_geo = empennage(S_h=S * Sh_S, S_v=S * Sv_S)
    nac_geo = nacelles(nr_of_engines=4)
    '''
    class airfoil_e334: #FLYING WING
        def __init__(self):                         # Study literature and add foil parameters, if applicable
            self.clmax = 1.4
            self.max_thickness = 0.119+d_thickness
            self.thickness_pos = 0.303+d_x_thickness
            self.max_camber = 0.04
            self.camber_pos = 0.254
            self.cm_0 = -0.05
    # Flying wing planform:
    wing_geo = wing(S=S,A=25.0+d_AR, qc_sweep=(15.0+d_sweep)*np.pi/180, taper=1.0+d_taper, dihedral=0.0*np.pi/180.0, airfoil=airfoil_e334())
    fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
    emp_geo = empennage(S_h = 0.0, S_v = 0.0)
    nac_geo = nacelles(nr_of_engines=4)
    '''

    MTOW = MTOW_initial

    # Compute initial error
    AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD,
                     Sh_S=Sh_S, Sv_S=Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt, energy_delta=energy_delta)

    powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass) / MTOW
    payload_frac = AHAPS.payload.mass_payload / MTOW
    error = abs(powM_frac - powM_frac_target) / powM_frac_target + abs(payload_frac - payload_apprx_frac) / payload_apprx_frac

    error_vec = np.ones(10) * error
    monitoring_var = np.linalg.norm(abs(error_vec - np.mean(error_vec)))

    iterations = 0
    while monitoring_var > 5e-3 or iterations < 10:
        AHAPS = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=lat, day_margin=day_margin, DoD=DoD,
                         Sh_S=Sh_S, Sv_S=Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt, energy_delta=energy_delta)

        powM_frac = (AHAPS.pow_store.mass + AHAPS.solar.mass) / MTOW
        payload_frac = AHAPS.payload.mass_payload / MTOW
        error = abs(powM_frac - powM_frac_target) / powM_frac_target + abs(payload_frac - payload_apprx_frac) / payload_apprx_frac

        MTOW = (AHAPS.pow_store.mass + AHAPS.solar.mass + AHAPS.payload.mass_payload) / (powM_frac_target + payload_apprx_frac)

        iterations += 1
        error_vec = np.roll(error_vec, 1)
        error_vec[0] = error
        monitoring_var = np.linalg.norm(abs(error_vec - np.mean(error_vec)))

    # Return the final values as a list
    return [
        AHAPS.MTOW,
        AHAPS.Pow_req,
        AHAPS.wing.S,
        AHAPS.solar.area,
        AHAPS.pow_store.mass,
        AHAPS.pow_store.volume,
        AHAPS.MTOW * (1 - powM_frac) - AHAPS.payload.mass_payload,
        AHAPS.prop.lambda_adv
    ]


labels = [
    "MTOW",
    "Power Consumption",
    "Wing Surface Area",
    "Solar Panel Area",
    "Energy Storage System Mass",
    "Energy Storage System Volume",
    "Remaining Mass (MTOW - Power System Mass - Payload Mass)",
    "Lambda Advance Ratio"
]
delta=0.01
# Compute results for two different aspect ratios
results_0 = compute()
results_1 = compute(d_AR=delta)

# Compute the difference between the two results and calculate derivatives
differences = [r1 - r0 for r1, r0 in zip(results_1, results_0)]
derivatives = [diff / delta for diff in differences]

# Print the derivatives with labels
print("\nDerivative with respect to aspect ratio:")
for label, derivative in zip(labels, derivatives):
    print(f"{label}: {derivative}")


results_1 = compute(d_sweep=delta)

# Compute the difference between the two results and calculate derivatives
differences = [r1 - r0 for r1, r0 in zip(results_1, results_0)]
derivatives = [diff / delta for diff in differences]

# Print the derivatives with labels
print("\nDerivative with respect to sweep:")
for label, derivative in zip(labels, derivatives):
    print(f"{label}: {derivative}")



results_1 = compute(d_taper=delta)

# Compute the difference between the two results and calculate derivatives
differences = [r1 - r0 for r1, r0 in zip(results_1, results_0)]
derivatives = [diff / delta for diff in differences]

# Print the derivatives with labels
print("\nDerivative with respect to taper:")
for label, derivative in zip(labels, derivatives):
    print(f"{label}: {derivative}")


results_1 = compute(d_thickness=delta)

# Compute the difference between the two results and calculate derivatives
differences = [r1 - r0 for r1, r0 in zip(results_1, results_0)]
derivatives = [diff / delta for diff in differences]

# Print the derivatives with labels
print("\nDerivative with respect to airfoil max thicknesss:")
for label, derivative in zip(labels, derivatives):
    print(f"{label}: {derivative}")



results_1 = compute(d_x_thickness=delta)

# Compute the difference between the two results and calculate derivatives
differences = [r1 - r0 for r1, r0 in zip(results_1, results_0)]
derivatives = [diff / delta for diff in differences]

# Print the derivatives with labels
print("\nDerivative with respect to airfoil max thickness position:")
for label, derivative in zip(labels, derivatives):
    print(f"{label}: {derivative}")