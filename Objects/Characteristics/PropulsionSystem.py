import aerosandbox as asb
import aerosandbox.numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize_scalar
import warnings

warnings.filterwarnings("ignore")

class PropulsionSystem:
    def __init__(self, v_inf_cruise, required_thrust_cruise, m_TO, S, CL_max):
        # ---------------------------------------------------------
        # Input config
        # ---------------------------------------------------------
        self.v_inf_cr = v_inf_cruise
        self.thrust_req_cr = required_thrust_cruise
        self.m_TO = m_TO  # Take-Off Mass (kg)
        self.S = S        # Wing Reference Area (m^2)
        self.CL_max = CL_max  # Maximum Lift Coefficient
        # ---------------------------------------------------------
        # Trade-off parameters
        # ---------------------------------------------------------
        self.airfoil_name = "sd7037"
        self.Nb = 2
        self.num_engines = 4
        
        # Drivetrain Efficiencies
        self.eta_motor = 0.90
        self.eta_esc = 0.99
        self.eta_elec = self.eta_motor * self.eta_esc  # Combined Electrical Efficiency
        
        # Cruise (18,288 m) 60,000ft
        self.atmo_cr = asb.Atmosphere(altitude=18288)
        self.rho_cr = self.atmo_cr.density()
        self.mach_cr = self.v_inf_cr / self.atmo_cr.speed_of_sound()
        
        # Take-Off (Sea Level)
        self.atmo_to = asb.Atmosphere(altitude=0)
        self.rho_to = self.atmo_to.density()
        
        # ---------------------------------------------------------
        # Geometry
        # ---------------------------------------------------------
        self.r = np.linspace(0.1, 1.0, 100)
        self.beta_07 = 20.5
        # The chord polynomial (b) will be multiplied by the sized Diameter (D) later
        self.b_unscaled = (0.084241 - 0.85789*self.r + 4.7176*self.r**2 - 9.6225*self.r**3 + 8.5004*self.r**4 - 2.7959*self.r**5)
        self.beta_deg = self.beta_07 + (0.4387 + 0.3040*self.r - 3.9616*self.r**2 + 5.1180*self.r**3 - 1.6284*self.r**4 - 0.3244*self.r**5) * (180/np.pi)
        
        # State Variables initialized
        self.optimal_J = None
        self.D = None

    def _evaluate_bemt(self, v_inf, rpm, D, rho, cl_interp, cd_interp):
        """
        BEMT solver for cruise and take-off conditions. Returns Thrust, Torque, Mechanical Power, and Efficiency.
        """
        n_rps = rpm / 60.0
        omega = 2 * np.pi * n_rps
        J = v_inf / (n_rps * D)
        
        R_abs = D / 2.0
        r_abs = self.r * R_abs
        b_actual = self.b_unscaled * D
        
        phi_0_rad = np.arctan(v_inf / (omega * r_abs))
        dT_dr = np.zeros_like(r_abs)
        dM_dr = np.zeros_like(r_abs)
        
        for i in range(len(r_abs)):
            def equilibrium_equation(alpha_guess):
                cl_val = cl_interp(alpha_guess)
                phi_guess_rad = np.radians(self.beta_deg[i] - alpha_guess)
                lhs = cl_val * (self.Nb * b_actual[i]) / (2 * np.pi * r_abs[i])
                rhs = 4 * np.sin(phi_guess_rad) * np.tan(phi_guess_rad - phi_0_rad[i])
                return lhs - rhs
            
            try:
                # Primary attempt: perfect mathematical root crossing
                true_alpha_deg = brentq(equilibrium_equation, -35, 80)
            except ValueError:
                # Fallback: Deep stall minimizer (Crucial for Take-Off)
                res = minimize_scalar(lambda a: abs(equilibrium_equation(a)), bounds=(-35, 80), method='bounded')
                true_alpha_deg = res.x
                
            try:
                true_phi_rad = np.radians(self.beta_deg[i] - true_alpha_deg)
                cl_local = cl_interp(true_alpha_deg)
                cd_local = cd_interp(true_alpha_deg)
                
                # Induction Parameter (No F)
                K = (cl_local * self.Nb * b_actual[i] * np.cos(true_phi_rad)) / (8 * np.pi * r_abs[i] * (np.sin(true_phi_rad)**2))
                K_effective = np.clip(K, -0.5, 0.5)
                one_plus_a = 1.0 / (1.0 - K_effective)
                velocity_term = (one_plus_a**2) / (np.sin(true_phi_rad)**2)
                
                # Prandtl Tip Loss (Applied directly to physical forces)
                f_prandtl = (self.Nb / 2.0) * (R_abs - r_abs[i]) / (R_abs * np.sin(true_phi_rad))
                f_prandtl = np.clip(f_prandtl, 0.0, 100.0)
                F_local = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))
                
                dT_dr[i] = b_actual[i] * velocity_term * (cl_local * np.cos(true_phi_rad) - cd_local * np.sin(true_phi_rad)) * F_local
                dM_dr[i] = b_actual[i] * velocity_term * (cl_local * np.sin(true_phi_rad) + cd_local * np.cos(true_phi_rad)) * r_abs[i] * F_local
                
            except ZeroDivisionError:
                dT_dr[i] = 0.0
                dM_dr[i] = 0.0
                
        # Integration
        F0 = 0.5 * rho * (v_inf**2) * self.Nb
        Total_Thrust = F0 * trapezoid(dT_dr, r_abs)
        Total_Torque = F0 * trapezoid(dM_dr, r_abs)
        Mech_Power = Total_Torque * omega
        
        C_T = Total_Thrust / (rho * n_rps**2 * D**4)
        C_P = Mech_Power / (rho * n_rps**3 * D**5)
        eta = (J * C_T) / C_P if C_P > 0 else 0.0
        
        return Total_Thrust, Total_Torque, Mech_Power, eta

    def run_full_analysis(self):
        """
        Cruise Optimization -> Propeller Sizing -> Take-Off Analysis.
        """
        
        # =========================================================
        # Highest Efficiency Cruise Advance Ratio Search
        # =========================================================
        alphas_cr = np.linspace(-30, 85, 250)
        aero_cr = asb.Airfoil(self.airfoil_name).get_aero_from_neuralfoil(alpha=alphas_cr, Re=100_000, mach=self.mach_cr)
        cl_interp_cr = interp1d(alphas_cr, aero_cr['CL'], kind='cubic', fill_value="extrapolate")
        cd_interp_cr = interp1d(alphas_cr, aero_cr['CD'], kind='cubic', fill_value="extrapolate")
        
        J_sweep = np.linspace(0.4, 0.9, 25)
        efficiencies = []
        D_ref = 2.0
        
        for J in J_sweep:
            rpm_sweep = (self.v_inf_cr / (J * D_ref)) * 60.0
            _, _, _, eta = self._evaluate_bemt(self.v_inf_cr, rpm_sweep, D_ref, self.rho_cr, cl_interp_cr, cd_interp_cr)
            efficiencies.append(eta)
            
        self.optimal_J = J_sweep[np.argmax(efficiencies)]
        
        # =========================================================
        # Propeller Sizing for Required Cruise Thrust
        # =========================================================
        thrust_target_per_prop = self.thrust_req_cr / self.num_engines
        rpm_ref = (self.v_inf_cr / (self.optimal_J * D_ref)) * 60.0
        
        # Reference values for scaling
        T_ref, M_ref, P_mech_ref, _ = self._evaluate_bemt(self.v_inf_cr, rpm_ref, D_ref, self.rho_cr, cl_interp_cr, cd_interp_cr)
        
        # Scaling Diameter according to Thrust
        self.D = D_ref * np.sqrt(thrust_target_per_prop / T_ref)
        
        # Re-evaluate
        cruise_rpm = (self.v_inf_cr / (self.optimal_J * self.D)) * 60.0
        T_cr, M_cr, P_mech_cr, eta_cr = self._evaluate_bemt(self.v_inf_cr, cruise_rpm, self.D, self.rho_cr, cl_interp_cr, cd_interp_cr)
        
        # motor + ESC losses
        P_elec_cr = P_mech_cr / self.eta_elec
        P_available_cr = P_mech_cr * eta_cr

        # Tip Mach check for cruise
        n_rps_cr = cruise_rpm / 60.0
        omega_cr = 2 * np.pi * n_rps_cr
        R = self.D / 2.0
        tip_tangential_cr = omega_cr * R
        tip_speed_cr = np.sqrt(tip_tangential_cr**2 + self.v_inf_cr**2)
        a_cr = self.atmo_cr.speed_of_sound()
        tip_mach_cr = tip_speed_cr / a_cr if a_cr > 0 else 0.0

        # ========================================================
        # mass estimate
        # ========================================================
        m_esc = 0.523 # esc kg
        m_motor = 0.8 # motor kg
        m_add = 0.4 # cables, rod, insulation, etc. 200 gram nacelle 100 gram cable 100 gram insulation, etc
        m_rod = 1200 * (0.025/2)**2*np.pi*0.5  # density * volume of a 0.5m long, 25mm diameter lightweight carbon rod, 300 grams
        m_hub = 0.4 # 20 percent of kg from CAD  -   assumption 400 g now
        m_blades = 2.89*0.7 * (self.D)/(2.1357)# 70 percent of kg from CAD
        m_total_per_engine = m_esc + m_motor + m_add + m_rod +  m_hub + m_blades
        m_total_all_engines = m_total_per_engine * self.num_engines


        # =========================================================
        # REPORT
        # =========================================================
        # print("\n=======================================================")
        # print(f"                        CRUISE")
        # print("=======================================================")
        # print(f"Design Airspeed          : {self.v_inf_cr} m/s")
        # print(f"Optimal Advance Ratio (J): {self.optimal_J:.3f}")
        # print(f"Sized Propeller Diameter : {self.D:.4f} m")
        # print(f"Cruise Rotational Speed  : {cruise_rpm:.0f} RPM")
        # print("-------------------------------------------------------")
        # print(f"Thrust per Prop          : {T_cr:.2f} N")
        # print(f"Torque per Prop          : {M_cr:.2f} Nm")
        # print(f"Mechanical Shaft Power   : {P_mech_cr:.2f} W")
        # print(f"Electrical Power Draw    : {P_elec_cr:.2f} W  (Motor={self.eta_motor}, ESC={self.eta_esc})")
        # print(f"Propeller Aerodynamic Eff: {eta_cr * 100:.2f} %")
        # print(f"Total Aircraft Elec Pwr  : {(P_elec_cr * self.num_engines) / 1000.0:.3f} kW")
        # print(f"Total Power Available at Propeller: {(P_available_cr * self.num_engines) / 1000.0:.3f} kW")
        # print(f"Tip Mach (cruise)        : {tip_mach_cr:.3f} {'OK' if tip_mach_cr < 0.7 else 'WARNING: >0.7'}")

        return P_elec_cr * self.num_engines, m_total_all_engines
        
if __name__ == "__main__":
    ahaps = PropulsionSystem(
        v_inf_cruise=32.45494362484863, 
        required_thrust_cruise=65.54910335953765, 
        m_TO=278.1416060277118, 
        S=55.79535859750816,
        CL_max= 1.0319550892283087
    )

    ahaps.run_full_analysis()
