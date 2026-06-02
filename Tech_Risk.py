import numpy as np
import matplotlib.pyplot as plt
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


# Compute normalized derivatives
delta = 0.01  # Small change in the input parameter

# Compute results for the baseline (no change)
results_0 = compute()

# Compute normalized derivatives
def compute_normalized_derivative(results_0, results_1, delta, x0, labels):
    normalized_derivatives = []
    for r0, r1 in zip(results_0, results_1):
        if r0 != 0 and x0 != 0:  # Avoid division by zero
            normalized_derivative = ((r1 - r0) / r0) / (delta / x0)  # Normalized as (dy/y0) / (dx/x0)
        else:
            normalized_derivative = 0  # If the baseline value is zero, set derivative to 0
        normalized_derivatives.append(normalized_derivative)
    return normalized_derivatives

# Baseline input values
aspect_ratio_baseline = 25.0
sweep_baseline = 0.0
taper_baseline = 1.0
thickness_baseline = 0.091
thickness_pos_baseline = 0.311

# Derivative with respect to aspect ratio
results_1 = compute(d_AR=delta)
normalized_derivatives1 = compute_normalized_derivative(results_0, results_1, delta, aspect_ratio_baseline, labels)

# Derivative with respect to sweep
results_1 = compute(d_sweep=30)
normalized_derivatives2 = compute_normalized_derivative(results_0, results_1, 30, sweep_baseline + 25, labels)  # Add small value to avoid division by zero

# Derivative with respect to taper
results_1 = compute(d_taper=-0.9)
normalized_derivatives3 = compute_normalized_derivative(results_0, results_1, -0.9, taper_baseline, labels)

# Derivative with respect to airfoil max thickness
results_1 = compute(d_thickness=10*delta)
normalized_derivatives4 = compute_normalized_derivative(results_0, results_1, 10*delta, thickness_baseline, labels)

# Derivative with respect to airfoil max thickness position
results_1 = compute(d_x_thickness=10*delta)
normalized_derivatives5 = compute_normalized_derivative(results_0, results_1, 10*delta, thickness_pos_baseline, labels)

# Print normalized derivatives
print("\nNormalized Derivative with respect to aspect ratio:")
for label, derivative in zip(labels, normalized_derivatives1):
    print(f"{label}: {derivative}")

print("\nNormalized Derivative with respect to sweep:")
for label, derivative in zip(labels, normalized_derivatives2):
    print(f"{label}: {derivative}")

print("\nNormalized Derivative with respect to taper:")
for label, derivative in zip(labels, normalized_derivatives3):
    print(f"{label}: {derivative}")

print("\nNormalized Derivative with respect to airfoil max thickness:")
for label, derivative in zip(labels, normalized_derivatives4):
    print(f"{label}: {derivative}")

print("\nNormalized Derivative with respect to airfoil max thickness position:")
for label, derivative in zip(labels, normalized_derivatives5):
    print(f"{label}: {derivative}")


sensitivity_data = np.array([
    normalized_derivatives1[:3],  # Sensitivity of MTOW, Power Consumption, Wing Surface Area to Aspect Ratio
    normalized_derivatives2[:3],  # Sensitivity to Sweep
    normalized_derivatives3[:3],  # Sensitivity to Taper
    normalized_derivatives4[:3],  # Sensitivity to Max Thickness
    normalized_derivatives5[:3]   # Sensitivity to Thickness Position
])

# Take the absolute value of the sensitivities
sensitivity_data = np.abs(sensitivity_data)



# Define a function to create a grouped bar plot for sensitivities
def create_combined_sensitivity_plot(sensitivity_data, input_labels, output_labels):
    num_inputs = len(input_labels)
    num_outputs = len(output_labels)
    bar_width = 0.15  # Width of each bar
    x = np.arange(num_inputs)  # X positions for the groups

    # Create the plot
    plt.figure(figsize=(12, 6))
    for i in range(num_outputs):
        plt.bar(x + i * bar_width, sensitivity_data[:, i], bar_width, label=output_labels[i])

    # Set axis labels and title
    plt.xlabel("Concept Design Inputs", fontsize=12)
    plt.ylabel("Sensitivity (|Normalized Derivative|)", fontsize=12)

    # Set x-axis ticks and labels
    plt.xticks(x + bar_width * (num_outputs - 1) / 2, input_labels, rotation=45, ha="right")

    # Add legend
    plt.legend(title="Outputs", fontsize=10)

    # Add grid for better readability
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Show the plot
    plt.tight_layout()
    plt.show()

# Input labels and output labels
input_labels = ["Aspect Ratio", "Sweep", "Taper", "Max Thickness", "Thickness Position"]
output_labels = ["MTOW", "Power Consumption", "Wing Surface Area"]

# Create the combined sensitivity plot
create_combined_sensitivity_plot(sensitivity_data, input_labels, output_labels)
