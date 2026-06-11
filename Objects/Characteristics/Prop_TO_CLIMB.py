import aerosandbox as asb
import aerosandbox.numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from pathlib import Path
import sys
import matplotlib.pyplot as plt

try:
	from Objects.Characteristics.PropulsionSystem import PropulsionSystem
except ModuleNotFoundError:
	# Allow running this file directly via: python Objects/Characteristics/Prop_TO_CLIMB.py
	repo_root = Path(__file__).resolve().parents[2]
	if str(repo_root) not in sys.path:
		sys.path.insert(0, str(repo_root))
	from Objects.Characteristics.PropulsionSystem import PropulsionSystem


# DESIGN INPUTS
D = 1.7045684187645866  # m (already optimized for cruise)
v_inf_cruise = 32.45494362484863  # m/s
required_thrust_cruise = 65.54910335953765  # N
m_TO = 278.1416060277118 + 10.0  # kg 10 for landing gear!!
S = 55.79535859750816  # m^2
CL_max = 1.0319550892283087  # -

v_initial = 5.0  # m/s human push-off speed

# Take-off model assumptions
f_LW = 1.2
mu_runway = 0.08 # raymer
C_D_TO = 0.024
P_TO_ELEC_PER_MOTOR = 1162.0  # W electrical per motor
DT = 0.05
MAX_TAKEOFF_TIME = 60.0


def build_takeoff_airfoil_interpolants(propulsion, v_ref=v_initial):
    mach_to = v_ref / propulsion.atmo_to.speed_of_sound()
    alphas_to = np.linspace(-30.0, 85.0, 300)
    aero_to = asb.Airfoil(propulsion.airfoil_name).get_aero_from_neuralfoil(
        alpha=alphas_to,
        Re=1_000_000,
        mach=mach_to,
    )
    cl_interp_to = interp1d(alphas_to, aero_to["CL"], kind="linear", fill_value="extrapolate")
    cd_interp_to = interp1d(alphas_to, aero_to["CD"], kind="linear", fill_value="extrapolate")
    return cl_interp_to, cd_interp_to


def solve_power_limited_takeoff_rpm(propulsion, diameter, cl_interp_to, cd_interp_to, v_eval=v_initial):
    # P_TO_ELEC_PER_MOTOR is motor input power (after ESC losses)
    # Shaft power = motor_input * eta_motor
    p_motor_input = P_TO_ELEC_PER_MOTOR
    p_to_mech_motor = p_motor_input * propulsion.eta_motor

    def power_residual(rpm_guess):
        _, _, p_mech_guess, _ = propulsion._evaluate_bemt(
            v_eval,
            rpm_guess,
            diameter,
            propulsion.rho_to,
            cl_interp_to,
            cd_interp_to,
        )
        return p_mech_guess - p_to_mech_motor

    rpm_min = 100.0
    rpm_max = 3000.0
    f_min = power_residual(rpm_min)
    f_max = power_residual(rpm_max)

    while f_max < 0.0 and rpm_max < 30000.0:
        rpm_max *= 1.5
        f_max = power_residual(rpm_max)

    if f_min > 0.0:
        return rpm_min
    if f_max < 0.0:
        raise ValueError("Power-limited RPM root not bracketed for take-off analysis.")

    return brentq(power_residual, rpm_min, rpm_max)


def simulate_takeoff_roll(propulsion, diameter, takeoff_rpm, cl_interp_to, cd_interp_to):
    c_l_to = 0.8 * propulsion.CL_max
    v_to = np.sqrt(f_LW * 2.0 * 9.81 * propulsion.m_TO / (c_l_to * propulsion.rho_to * propulsion.S))

    # Constant-thrust assumption: freeze thrust at initial roll speed.
    t_single_const, _, _, _ = propulsion._evaluate_bemt(
        v_initial,
        takeoff_rpm,
        diameter,
        propulsion.rho_to,
        cl_interp_to,
        cd_interp_to,
    )
    total_thrust_const = t_single_const * propulsion.num_engines

    v_ground = v_initial
    x_distance = 0.0
    t_time = 0.0

    hist_t = []
    hist_v = []
    hist_x = []
    hist_thrust = []
    hist_drag = []
    hist_friction = []
    hist_lift = []
    hist_accel = []

    failed = False
    failure_speed = None

    while v_ground < v_to and t_time < MAX_TAKEOFF_TIME:
        lift = 0.5 * propulsion.rho_to * (v_ground**2) * propulsion.S * c_l_to
        drag = 0.5 * propulsion.rho_to * (v_ground**2) * propulsion.S * C_D_TO
        wow = max(0.0, propulsion.m_TO * 9.81 - lift)
        friction = mu_runway * wow

        net_force = total_thrust_const - drag - friction

        if net_force <= 0.0:
            failed = True
            failure_speed = v_ground
            break

        acceleration = net_force / propulsion.m_TO

        hist_t.append(t_time)
        hist_v.append(v_ground)
        hist_x.append(x_distance)
        hist_thrust.append(total_thrust_const)
        hist_drag.append(drag)
        hist_friction.append(friction)
        hist_lift.append(lift)
        hist_accel.append(acceleration)

        v_ground += acceleration * DT
        x_distance += v_ground * DT
        t_time += DT

    t_to, m_to, p_mech_to, eta_to = propulsion._evaluate_bemt(
        v_to,
        takeoff_rpm,
        diameter,
        propulsion.rho_to,
        cl_interp_to,
        cd_interp_to,
    )
    # Battery power = shaft_power / (eta_motor * eta_esc)
    # But more directly: if motor_input was 1162 W, battery = 1162 / eta_esc
    p_motor_input_to = p_mech_to / propulsion.eta_motor
    p_battery_to = p_motor_input_to / propulsion.eta_esc

    n_rps_to = takeoff_rpm / 60.0
    omega_to = 2.0 * np.pi * n_rps_to
    tip_tangential_to = omega_to * (diameter / 2.0)
    tip_speed_to = np.sqrt(tip_tangential_to**2 + v_to**2)
    a_to = propulsion.atmo_to.speed_of_sound()
    tip_mach_to = tip_speed_to / a_to if a_to > 0 else 0.0

    return {
        "v_to": float(v_to),
        "distance": float(x_distance),
        "time": float(t_time),
        "failed": failed,
        "failure_speed": None if failure_speed is None else float(failure_speed),
        "thrust_total_constant": float(total_thrust_const),
        "thrust_per_prop_constant": float(t_single_const),
        "thrust_liftoff_per_prop": float(t_to),
        "torque_liftoff_per_prop": float(m_to),
        "power_mech_liftoff_per_prop": float(p_mech_to),
        "power_motor_input_per_prop": float(p_motor_input_to),
        "power_battery_per_prop": float(p_battery_to),
        "power_battery_total": float(p_battery_to * propulsion.num_engines),
        "eta_liftoff": float(eta_to),
        "tip_mach_liftoff": float(tip_mach_to),
        "hist_t": hist_t,
        "hist_v": hist_v,
        "hist_x": hist_x,
        "hist_thrust": hist_thrust,
        "hist_drag": hist_drag,
        "hist_friction": hist_friction,
        "hist_lift": hist_lift,
        "hist_accel": hist_accel,
    }


def plot_takeoff_dashboard(result, mass_kg):
    fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    ax1 = axs[0]
    ax1.plot(result["hist_t"], result["hist_v"], "b-", label="Velocity (m/s)", linewidth=2.2)
    ax1.set_ylabel("Velocity (m/s)", color="b", fontsize=11)
    ax1.tick_params(axis="y", labelcolor="b")
    ax1.grid(True, alpha=0.3)
    ax1_twin = ax1.twinx()
    ax1_twin.plot(result["hist_t"], result["hist_x"], "g--", label="Distance (m)", linewidth=2.2)
    ax1_twin.set_ylabel("Distance (m)", color="g", fontsize=11)
    ax1_twin.tick_params(axis="y", labelcolor="g")
    ax1.set_title("Take-Off Kinematics: Velocity and Ground Roll", fontsize=13, fontweight="bold")
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_1t, labels_1t = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_1t, labels_1 + labels_1t, loc="upper left")

    ax2 = axs[1]
    ax2.plot(result["hist_t"], result["hist_thrust"], "r-", label="Total Propeller Thrust", linewidth=2.2)
    ax2.plot(result["hist_t"], result["hist_drag"], "k-", label="Aerodynamic Drag", linewidth=1.8)
    ax2.plot(result["hist_t"], result["hist_friction"], "y-", label="Rolling Friction", linewidth=1.8)
    ax2.set_ylabel("Force (N)", fontsize=11)
    ax2.set_title("Forces During Take-Off", fontsize=13, fontweight="bold")
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    ax3 = axs[2]
    ax3.plot(result["hist_t"], result["hist_accel"], "m-", label="Acceleration (m/s^2)", linewidth=2.2)
    ax3.set_ylabel("Acceleration (m/s^2)", color="m", fontsize=11)
    ax3.tick_params(axis="y", labelcolor="m")
    ax3.set_xlabel("Time (s)", fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3_twin = ax3.twinx()
    ax3_twin.plot(result["hist_t"], result["hist_lift"], "c-", label="Aerodynamic Lift", linewidth=2.2)
    ax3_twin.axhline(mass_kg * 9.81, color="k", linestyle=":", label="Aircraft Weight")
    ax3_twin.set_ylabel("Force (N)", color="c", fontsize=11)
    ax3_twin.tick_params(axis="y", labelcolor="c")
    ax3.set_title("Acceleration and Lift Buildup", fontsize=13, fontweight="bold")
    lines_3, labels_3 = ax3.get_legend_handles_labels()
    lines_3t, labels_3t = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines_3 + lines_3t, labels_3 + labels_3t, loc="center right")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    propulsion = PropulsionSystem(
        v_inf_cruise=v_inf_cruise,
        required_thrust_cruise=required_thrust_cruise,
        m_TO=m_TO,
        S=S,
        CL_max=CL_max,
    )

    cruise_power_total, _ = propulsion.run_full_analysis()
    propulsion.D = D

    cl_interp_to, cd_interp_to = build_takeoff_airfoil_interpolants(propulsion)
    takeoff_rpm = solve_power_limited_takeoff_rpm(propulsion, D, cl_interp_to, cd_interp_to)
    result = simulate_takeoff_roll(propulsion, D, takeoff_rpm, cl_interp_to, cd_interp_to)

#     print("================ Take-off Computation ================")
#     print(f"Fixed Propeller Diameter      : {D:.4f} m")
#     print(f"Power-Limited Take-off RPM    : {takeoff_rpm:.0f} RPM")
#     print(f"Lift-off Speed                : {result['v_to']:.2f} m/s")
#     print(f"Ground Roll Distance          : {result['distance']:.1f} m")
#     print(f"Take-off Time                 : {result['time']:.2f} s")
#     print("------------------------------------------------------")
#     print(f"thrust per propeller at liftoff : {result['thrust_liftoff_per_prop']:.2f} N")
#     print(f"Liftoff Shaft Power per Prop  : {result['power_mech_liftoff_per_prop']:.2f} W")
#     print(f"Liftoff Motor Input per Prop  : {result['power_motor_input_per_prop']:.2f} W")
#     print(f"Liftoff Battery Power per Prop: {result['power_battery_per_prop']:.2f} W")
#     print(f"Liftoff Total Battery Power   : {result['power_battery_total'] / 1000.0:.3f} kW")
#     print(f"Liftoff Tip Mach              : {result['tip_mach_liftoff']:.3f}")
#     if result["failed"]:
#         print(f"WARNING: Acceleration stalled at {result['failure_speed']:.2f} m/s")
#     print("======================================================")
#     print(f"Reference Cruise Total Elec Power: {cruise_power_total / 1000.0:.3f} kW")

#     plot_takeoff_dashboard(result, m_TO)

# print(f"\nBattery power at liftoff: {result['power_battery_total']/1000.0:.3f} kW")


# =============================================================
# CLimb Performance Analysis
# =============================================================

def evaluate_climb_state(propulsion, diameter, v_climb, altitude_m, p_battery_per_motor):
    """
    I hope this works
    """
    # atmospheric conditions at this altitude
    atmo = asb.Atmosphere(altitude=altitude_m)
    rho_climb = atmo.density()
    a_climb = atmo.speed_of_sound()
    
    current_reynolds = np.interp(altitude_m, [0.0, 18288.0], [1_000_000, 100_000])
    current_mach = v_climb / a_climb

    alphas = np.linspace(-30.0, 85.0, 300)
    aero = asb.Airfoil(propulsion.airfoil_name).get_aero_from_neuralfoil(
        alpha=alphas,
        Re=current_reynolds,
        mach=current_mach,
    )
    cl_interp_local = interp1d(alphas, aero["CL"], kind="linear", fill_value="extrapolate")
    cd_interp_local = interp1d(alphas, aero["CD"], kind="linear", fill_value="extrapolate")

    # eff
    p_motor_input = p_battery_per_motor * propulsion.eta_esc
    p_mech_target = p_motor_input * propulsion.eta_motor
    
    # rpm determination
    def power_residual(rpm_guess):
        _, _, p_mech_guess, _ = propulsion._evaluate_bemt(
            v_climb,
            rpm_guess,
            diameter,
            rho_climb,
            cl_interp_local,
            cd_interp_local
        )
        return p_mech_guess - p_mech_target

    # 4. Bracket and solve for the optimal RPM
    rpm_min = 100.0
    rpm_max = 5000.0
    
    optimal_rpm = brentq(power_residual, rpm_min, rpm_max)
        
    # 5. Evaluate BEMT one last time at the optimal RPM to get final thrust
    t_single, m_single, p_mech_final, eta_prop = propulsion._evaluate_bemt(
        v_climb,
        optimal_rpm,
        diameter,
        rho_climb,
        cl_interp_local,
        cd_interp_local
    )
    
    total_thrust = t_single * propulsion.num_engines
    
    power_available = total_thrust * v_climb 

    # Calculate Tip Mach Number
    n_rps = optimal_rpm / 60.0
    omega = 2.0 * np.pi * n_rps
    tip_tangential = omega * (diameter / 2.0)
    tip_speed = np.sqrt(tip_tangential**2 + v_climb**2)
    tip_mach = tip_speed / a_climb


    return {
        "altitude_m": altitude_m,
        "rpm": optimal_rpm,
        "thrust_total": total_thrust,
        "power_available": power_available,
        "propeller_efficiency": eta_prop,
        "reynolds": current_reynolds, 
        "tip_mach": tip_mach          
    }

# --- Setup Inputs ---
# 80 percent of take-off
TO_BATTERY_PER_MOTOR = result['power_battery_total'] / 4
CLIMB_BATTERY_PER_MOTOR = TO_BATTERY_PER_MOTOR * 0.80 

# IAN change the speed!!!!
v_climb_target = 30.0  # m/s

# Altitudes from 100 to 60,000 ft (converted to meters)
altitudes_ft = np.linspace(100, 60000, 20)
altitudes_m = altitudes_ft * 0.3048

# --- Compute Sweep ---
climb_results = []

print(f"--- Climbing at {v_climb_target} m/s with {CLIMB_BATTERY_PER_MOTOR:.2f} W battery power per motor ---")
print(f"{'Alt (ft)':>8} | {'Alt (m)':>7} | {'RPM':>5} | {'Thrust (N)':>10} | {'Power (W)':>10} | {'Prop Eta':>8} | {'Reynolds':>8} | {'Tip Mach':>8}")
print("-" * 84)

for alt_m, alt_ft in zip(altitudes_m, altitudes_ft):
    res = evaluate_climb_state(
        propulsion=propulsion,
        diameter=D,
        v_climb=v_climb_target,
        altitude_m=alt_m,
        p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR
    )
    
    # Updated 
    print(f"{alt_ft:8.0f} | {alt_m:7.0f} | {res['rpm']:5.0f} | {res['thrust_total']:10.2f} | {res['power_available']:10.2f} | {res['propeller_efficiency']:8.3f} | {res['reynolds']:8.0f} | {res['tip_mach']:8.3f}")
    
    climb_results.append(res)