from pathlib import Path

import numpy as np
import ambiance as am
import matplotlib.pyplot as plt
import pandas as pd
#Get Data Points First


class PropulsionSystem:
    def __init__(self, plotdata=False, T=0.0, velocity=0.0, alt=0.0, rpm=0.0, torque=0.0, motor_temp=0.0, propeller_diameter=2.5):                                     # Initialise with proper values
        self.mass = None  # kg, estimated mass of the propulsion system
        self.alt = alt  # Altitude in meters
        self.velocity = velocity  # m/s, cruise airspeed
        self.T = T
        
        #Motor characteristics
        self.motor_temp = motor_temp
        
        # self.q = 0  # Dynamic pressure, to be calculated based on altitude and velocity
        # Propeller characteristics for a HAPS-scale low-speed propeller
        self.propeller_diameter = propeller_diameter  # meters
        self.propeller_area = np.pi * (self.propeller_diameter / 2) ** 2
        self.velocity = velocity  # m/s, cruise airspeed
        self.lambda_adv = None  # Advance ratio, to be calculated based on velocity and propeller characteristics
        
        #Efficiencies
        self.rpm_out = rpm  # Initialize RPM output (assuming direct drive, no gearbox)
        self.motor_eff = self.calc_motor_eff(motor_temp, rpm, torque, plotdata)
        self.gearbox_eff = self.calc_gearbox_eff()
        self.propeller_eff = self.calc_propeller_eff()
        self.overall_eff = self.calc_overall_eff()

        #Power Required
        self.power_required = self.Calc_Power_Req()
    
    def calc_motor_eff(self, motor_temp, rpm, torque, plotdata):                          # Compute all relevant characteristics of the propulsion system
        def convex_hull(points):
            pts = np.unique(points, axis=0)
            if len(pts) <= 2:
                return pts

            pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

            def cross(o, a, b):
                return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

            lower = []
            for p in pts:
                while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                    lower.pop()
                lower.append(p)

            upper = []
            for p in reversed(pts):
                while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                    upper.pop()
                upper.append(p)

            return np.array(lower[:-1] + upper[:-1])

        # Build closed contour polygons from datapoints and classify the input point.
        data_dir = Path(__file__).resolve().parent / 'MotorData'
        data_94_minus40C = pd.read_csv(data_dir / 'minus40_94.csv', delimiter=';')
        data_93_plus20C = pd.read_csv(data_dir / 'plus20_93.csv', delimiter=';')
        x_minus40C = data_94_minus40C['RPM'].values
        y_minus40C = data_94_minus40C['TORQUE'].values
        x_plus20C = data_93_plus20C['RPM'].values
        y_plus20C = data_93_plus20C['TORQUE'].values
        from matplotlib.path import Path as MplPath

        points_minus40C = np.column_stack((x_minus40C, y_minus40C))
        points_plus20C = np.column_stack((x_plus20C, y_plus20C))

        hull_minus40C = convex_hull(points_minus40C)
        hull_plus20C = convex_hull(points_plus20C)

        contour_minus40C = MplPath(hull_minus40C, closed=True)
        contour_plus20C = MplPath(hull_plus20C, closed=True)

        if motor_temp == -40:
            selected_contour = contour_minus40C
            selected_eff = 0.94
            selected_hull = hull_minus40C
            selected_points = points_minus40C
            selected_label = '-40C contour (0.94)'
            selected_color = 'blue'
        elif motor_temp == 20:
            selected_contour = contour_plus20C
            selected_eff = 0.93
            selected_hull = hull_plus20C
            selected_points = points_plus20C
            selected_label = '+20C contour (0.93)'
            selected_color = 'orange'
        else:
            raise ValueError('motor_temp must be either -40 or 20.')

        input_point = np.array([[rpm, torque]])
        in_selected_contour = selected_contour.contains_points(input_point)[0]
        motor_eff = selected_eff if in_selected_contour else 0.0

        if plotdata:
            plt.figure()
            plt.fill(selected_hull[:, 0], selected_hull[:, 1], color=selected_color, alpha=0.2, label=selected_label)
            plt.plot(
                np.r_[selected_hull[:, 0], selected_hull[0, 0]],
                np.r_[selected_hull[:, 1], selected_hull[0, 1]],
                color=selected_color
            )
            plt.scatter(selected_points[:, 0], selected_points[:, 1], color=selected_color, s=18, alpha=0.7)
            plt.scatter([rpm], [torque], color='red', s=40, label=f'Input point -> eff={motor_eff:.2f}')
            plt.xlabel('RPM')
            plt.ylabel('Torque')
            plt.title(f'Motor Efficiency Contour at {motor_temp}C')
            plt.legend()
            plt.grid()
            plt.show()
        return motor_eff
    
    def calc_gearbox_eff(self):
        gearbox_eff = 0.95
        return gearbox_eff
    
    def calc_propeller_eff(self):
        """
        Calculate propeller efficiency using Truckenbrodt 1999 formula (eq. 3.5).
         
        η_prop ≈ 2·(1-λ²·ln(1 + 1/λ²)) / (1 + √(1 + T/(q·A) - 2·λ²·ln(1 + 1/λ²)))
        
        where:
            λ = advance ratio = V / (n·D)
            q = dynamic pressure = 0.5·ρ·V²
            T = thrust
            A = propeller disk area
        """
        # Calculate parameters
        n_hz = self.rpm_out / 60.0  # rotational frequency [Hz]
        self.lambda_adv = self.velocity / (n_hz * self.propeller_diameter * np.pi)  # advance ratio λ = V / (n·D)
        q = 0.5 * am.Atmosphere(self.alt).density[0] * self.velocity ** 2  # dynamic pressure [Pa]
        
        # Compute ln term (appears twice in formula)
        ln_term = np.log(1.0 + 1.0 / (self.lambda_adv ** 2))
        
        # Numerator: 2·(1 - λ²·ln(...))
        numerator = 2.0 * (1.0 - self.lambda_adv ** 2 * ln_term)
        
        # Denominator: 1 + √(1 + T/(q·A) - 2·λ²·ln(...))
        # Disc loading: T/(q·A) [Pa]
        disc_loading = self.T / (q * self.propeller_area) if q > 0 else 0.0
        inner = 1.0 + disc_loading - 2.0 * self.lambda_adv ** 2 * ln_term     
        denominator = 1.0 + np.sqrt(inner)
        
        eta_prop = numerator / denominator
        
        # Clamp to valid range [0, 1]
        return float(eta_prop)
    
    def calc_overall_eff(self):
        # Implement method to calculate overall efficiency based on motor, gearbox, and propeller efficiencies
        if self.motor_eff is not None and self.gearbox_eff is not None and self.propeller_eff is not None:
            return self.motor_eff * self.gearbox_eff * self.propeller_eff
        else:
            return None
        
    def Calc_Power_Req(self):
        # Implement method to calculate power required based on thrust and velocity
        return float(self.T * self.velocity / self.overall_eff)
    
if __name__ == "__main__":
    #Example Inputs
    velocity = 25.0  # m/s, cruise airspeed
    alt = 15000 # m, altitude
    propeller_diameter = 2.5  # m
    rpm = 1000
    torque = 4  # Nm, torque of the motor
    motor_temp = -40  # °C, motor temperature
    gamma = 0  # degrees, flight path angle
    W = 1000  # N, weight of the aircraft
    CD = 0.04 # Drag coefficient
    S = 36.0  # m², reference area for drag calculation

    # Example usage
    propulsion_system = PropulsionSystem(plotdata=True, gamma=gamma, W=W, velocity=velocity, alt=alt, rpm=rpm, torque=torque, motor_temp=motor_temp, CD=CD, S=S, propeller_diameter=propeller_diameter)
    
    # Print calculated efficiencies
    print("Propulsion System Characteristics:")
    print(f"Motor Efficiency: {propulsion_system.motor_eff:.2f}")
    print(f"Gearbox Efficiency: {propulsion_system.gearbox_eff:.2f}")
    print(f"Propeller Efficiency: {propulsion_system.propeller_eff:.2f}")
    print(f"Overall Efficiency: {propulsion_system.overall_eff:.2f}")
    print(f"Power Required: {propulsion_system.power_required:.2f} W")
    print(f"Advance Ratio (λ): {propulsion_system.lambda_adv:.2f}")
    
    #Plot Changes in Alt with same velocity
    altitudes = np.linspace(0, 20000, 5)  # Alt
    powers = []
    efficiencies = []
    for alt in altitudes:
        propulsion_system.alt = alt
        propulsion_system.T = propulsion_system.calc_Thrust()
        propulsion_system.propeller_eff = propulsion_system.calc_propeller_eff()
        propulsion_system.overall_eff = propulsion_system.calc_overall_eff()
        power_req = propulsion_system.Calc_Power_Req()
        powers.append(power_req)
        efficiencies.append(propulsion_system.propeller_eff)
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(altitudes, powers, marker='o')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power Required (W)')
    plt.title('Power Required vs Altitude at Constant Velocity')
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.plot(altitudes, efficiencies, marker='o')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Propeller Efficiency')
    plt.title('Propeller Efficiency vs Altitude at Constant Velocity')
    plt.grid()
    plt.show()

    #plot changes in velocity with same altitude
    velocities = np.linspace(10, 50, 5)  # m/s
    powers = []
    efficiencies = []
    for velocity in velocities:
        propulsion_system.velocity = velocity
        propulsion_system.T = propulsion_system.calc_Thrust()
        propulsion_system.propeller_eff = propulsion_system.calc_propeller_eff()
        propulsion_system.overall_eff = propulsion_system.calc_overall_eff()
        power_req = propulsion_system.Calc_Power_Req()
        powers.append(power_req)
        efficiencies.append(propulsion_system.propeller_eff)
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(velocities, powers, marker='o')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power Required (W)')
    plt.title('Power Required vs Velocity at Constant Altitude')
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.plot(velocities, efficiencies, marker='o')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Propeller Efficiency')
    plt.title('Propeller Efficiency vs Velocity at Constant Altitude')
    plt.grid()
    plt.show()

