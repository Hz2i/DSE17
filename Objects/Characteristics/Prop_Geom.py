# import library
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import warnings

warnings.filterwarnings("ignore")

# ===========================================================================
# same initialization as Prop_Airfoil.py but wrapped in a function for reuse
# ===========================================================================
v_inf = 27.6      # Freestream velocity in m/s
altitude = 18288  # 60,000 ft operating altitude for HAPS
atmo = asb.Atmosphere(altitude=altitude)
rho = atmo.density() 
speed_of_sound = atmo.speed_of_sound()
mach_number = v_inf / speed_of_sound

# Propeller Architecture
Nb = 2            # Number of blades
beta_07 = 20.5    # Blade pitch setting at 70% radius (degrees)
r = np.linspace(0.10, 1.0, 100) # Active aerodynamic span grid (r/R)

def size_propeller_diameter(J_target, airfoil_name, total_thrust_required, num_engines):
    """
    Analytically sizes the propeller diameter to meet a target thrust requirement
    at a specific advance ratio J using fully coupled BEMT physics.
    """
    # Calculate required thrust per individual propeller
    thrust_target_per_prop = total_thrust_required / num_engines
    
    # AF polars
    alphas_sweep = np.linspace(-30, 30, 250)
    airfoil = asb.Airfoil(airfoil_name)
    aero_data = airfoil.get_aero_from_neuralfoil(
        alpha=alphas_sweep,
        Re=100_000,
        mach=mach_number
    )
    cl_interp = interp1d(alphas_sweep, aero_data['CL'], kind='cubic', bounds_error=False, fill_value="extrapolate")
    cd_interp = interp1d(alphas_sweep, aero_data['CD'], kind='cubic', bounds_error=False, fill_value="extrapolate")

    # --- Performance at ref 2.0 m diameter ---
    D_ref = 2.0
    R_ref = D_ref / 2.0
    r_abs_ref = r * R_ref
    
    # Geometry scaled to reference dimensions
    b_ref = (0.084241 - 0.85789*r + 4.7176*r**2 - 9.6225*r**3 + 8.5004*r**4 - 2.7959*r**5) * D_ref
    beta_deg = beta_07 + (0.4387 + 0.3040*r - 3.9616*r**2 + 5.1180*r**3 - 1.6284*r**4 - 0.3244*r**5) * (180/np.pi)
    
    n_ref = v_inf / (D_ref * J_target)
    omega_ref = 2 * np.pi * n_ref
    phi_0_rad = np.arctan(v_inf / (omega_ref * r_abs_ref))
    
    dT_dr_ref = np.zeros_like(r_abs_ref)
    dM_dr_ref = np.zeros_like(r_abs_ref)
    
    for i in range(len(r_abs_ref)):
        def equilibrium_equation(alpha_guess):
            cl_val = cl_interp(alpha_guess)
            phi_guess_rad = np.radians(beta_deg[i] - alpha_guess)
            
            # Prandtl Tip Loss Calculation
            f_prandtl = (Nb / 2.0) * (R_ref - r_abs_ref[i]) / (R_ref * np.sin(phi_guess_rad))
            f_prandtl = np.clip(f_prandtl, 0.0, 100.0)
            F = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))
            
            lhs = cl_val * (Nb * b_ref[i]) / (2 * np.pi * r_abs_ref[i])
            rhs = 4 * F * np.sin(phi_guess_rad) * np.tan(phi_guess_rad - phi_0_rad[i])
            return lhs - rhs
        
        try:
            true_alpha_deg = brentq(equilibrium_equation, -25, 25)
            true_phi_rad = np.radians(beta_deg[i] - true_alpha_deg)
            
            cl_local = cl_interp(true_alpha_deg)
            cd_local = cd_interp(true_alpha_deg)
            
            
            # Compute Momentum parameter
            K = (cl_local * Nb * b_ref[i] * np.cos(true_phi_rad)) / (8 * np.pi * r_abs_ref[i] * (np.sin(true_phi_rad)**2))
            K_effective = np.clip(K ,-0.5, 0.5)
            one_plus_a = 1.0 / (1.0 - K_effective)
            
            velocity_term = (one_plus_a**2) / (np.sin(true_phi_rad)**2)

            # local tip loss factor for torque and thrust
            f_prandtl = (Nb / 2.0) * (R_ref - r_abs_ref[i]) / (R_ref * np.sin(true_phi_rad))
            f_prandtl = np.clip(f_prandtl, 0.0, 100.0)
            F_local = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_prandtl), 0.0, 1.0))
            
            # Cartesian integration elements
            dT_dr_ref[i] = b_ref[i] * velocity_term * (cl_local * np.cos(true_phi_rad) - cd_local * np.sin(true_phi_rad)) * F_local
            dM_dr_ref[i] = b_ref[i] * velocity_term * (cl_local * np.sin(true_phi_rad) + cd_local * np.cos(true_phi_rad)) * r_abs_ref[i] * F_local
        except (ValueError, ZeroDivisionError):
            dT_dr_ref[i] = 0.0
            dM_dr_ref[i] = 0.0

    F0_ref = 0.5 * rho * (v_inf**2) * Nb
    T_ref = F0_ref * trapezoid(dT_dr_ref, r_abs_ref)
    M_ref = F0_ref * trapezoid(dM_dr_ref, r_abs_ref)
    
    # --- Scaling according to d equations ---
    D_actual = D_ref * np.sqrt(thrust_target_per_prop / T_ref)
    
    # --- Scaling according to n equations ----
    n_actual = v_inf / (D_actual * J_target)
    rpm_actual = n_actual * 60.0
    omega_actual = 2 * np.pi * n_actual

    # --- Scaling according to M equations ---
    M_actual = M_ref * (D_actual / D_ref)**3
    Power_actual = M_actual * omega_actual
    
    # Compute Dimensionless Parameters
    C_T = thrust_target_per_prop / (rho * (n_actual**2) * (D_actual**4))
    C_P = Power_actual / (rho * (n_actual**3) * (D_actual**5))
    efficiency = (J_target * C_T) / C_P
    
    return {
        'Diameter (m)': D_actual,
        'RPM': rpm_actual,
        'Thrust per Prop (N)': thrust_target_per_prop,
        'Torque per Prop (Nm)': M_actual,
        'Power per Prop (W)': Power_actual,
        'Total Aircraft Power (kW)': (Power_actual * num_engines) / 1000.0,
        'C_T': C_T,
        'C_P': C_P,
        'Efficiency (%)': efficiency * 100.0
    }

# ===========================================================================
# !!!!!!SPECS!!!!!!!
# ===========================================================================
sizing_results = size_propeller_diameter(
    J_target=0.76, 
    airfoil_name="SD7037", 
    total_thrust_required=52.0, 
    num_engines=4
)

print("\n=======================================================")
print(f"      PROPULSION SIZING RESULTS: SD7037 PROFILE")
print("=======================================================")
for key, val in sizing_results.items():
    if '%' in key:
        print(f"{key:<28} : {val:.2f}%")
    elif 'kW' in key:
        print(f"{key:<28} : {val:.3f} kW")
    else:
        print(f"{key:<28} : {val:.4f}")
print("=======================================================")