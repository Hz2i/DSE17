import numpy as np
import aerosandbox as asb
import matplotlib.pyplot as plt
import sympy as sp

import sys
import os
# Add the folder containing characteristics_airframe.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

# Now you can import the Airframe module
from Airframe import airframe
import Components_Materials

points_loads = 100

def airfoil_properties(airfoil, chord_length=1.2):

    # Get the upper and lower coordinates
    upper_coords = airfoil.upper_coordinates()  # Upper surface (x, y)
    lower_coords = airfoil.lower_coordinates()  # Lower surface (x, y)
    x=np.linspace(0,1,51)
    interp_higher= np.interp(x,upper_coords[:,0][::-1],upper_coords[:,1][::-1])
    interp_lower = np.interp(x,lower_coords[:,0],lower_coords[:,1])
    y=interp_higher-interp_lower

    # Calculate the cross-sectional area using the Shoelace formula
    area=np.trapezoid(y,x)
    # Calculate the perimeter
    dx = np.diff(x)

    dy_upper = np.diff(interp_higher)
    dy_lower = np.diff(interp_lower)

    upper_length = np.sum(np.sqrt(dx**2 + dy_upper**2))
    lower_length = np.sum(np.sqrt(dx**2 + dy_lower**2))

    perimeter = upper_length + lower_length

    # Scale area and perimeter by the chord length
    scaled_area = area * chord_length**2
    scaled_perimeter = perimeter * chord_length
    # Return the results as a dictionary
    return scaled_area,scaled_perimeter
    


def internal_loading():             # Implement method for computing internal loading (lift distribution potentially taken from VLM analysis)
    pass

def pos_first_connection(airframe):
    dz = airframe.b/(points_loads*2)

    wingbox_length = 0.4329 #come from container
    sections_length = 4
    connection_length = 0.09 # is what bas told me, but the effect is probably negligible anyways due to resolution
    x_max_connection = int(np.round((wingbox_length +sections_length/2-connection_length/2)/dz,0))
    return x_max_connection

def bending_stress(Bending_distribution, airframe):    # Compute bending stresses from bending distribution
    CFRP = Components_Materials.CFRP()
    GLARE = Components_Materials.GLARE()
    safety_factor = 5
    yield_stress = CFRP.sigma
    chord=airframe.c_r
    max_thickness = airframe.foil.max_thickness()
    min_I = safety_factor*Bending_distribution[0]*chord*max_thickness/2/(yield_stress) # calcaulte max normal stress by using airfoil thickness/2 as max possible y 
    # find position where connection first starts taking (all) load
    x_max_connection=pos_first_connection(airframe)
    yield_stress = GLARE.sigma
    min_connection = safety_factor*Bending_distribution[x_max_connection]*chord*max_thickness/2/(yield_stress) # find I for connection
    return min_I, min_connection

def bending_deflection(Bending_distribution, airframe): # find bending deflection in either direction
    spar_I, connection_I = bending_stress(Bending_distribution, airframe) # get I from bending loads
    CFRP = Components_Materials.CFRP()
    I = min(spar_I, connection_I) # conservative estimate, difference is small as connections length are minimal, with small change in I
    # Compute deflection using beam theory

    dz = airframe.b/(points_loads*2*np.cos(airframe.qc_sweep))
    print(airframe.qc_sweep)
    dv2dz2 = Bending_distribution / (CFRP.E * I)
    plt.plot()
    dvdz = np.cumsum(dv2dz2) * dz
    z = np.cumsum(dvdz) * dz

    #plt.plot(np.linspace(0,len(theta)*dz,len(theta), z)
    #plt.ylim((0,len(z)*dz))
    #plt.show()
    # plot if 1:1 axes to see realistic deflection
    return z[-1]


def torsional_stress(Torsion_distribution, airframe): # Compute torsional stresses from torsion distribution
    PETa = Components_Materials.PET()
    safety_factor = 5
    max_shear = PETa.shear
    A_skin,_ = airfoil_properties(airframe.foil, airframe.c_r)
    t_skin = safety_factor*Torsion_distribution[0]/(2*A_skin*(max_shear))
    # 5x safety factor torque, min skin thickness calculated if skin carries all torque
    return t_skin

def twist_deflection(Torsion_distribution, airframe, r_thickness=0.002, r_spar=0.04): # Compute twist deflection from torsion distribution
    # i guessed r_thickness and r_spar but you need to get those from bas
    t_skin = torsional_stress(Torsion_distribution, airframe) # get area from torsional loads
    PET = Components_Materials.PET()
    CFRP = Components_Materials.CFRP()
    # considering all torsion carried by skin
    airfoil_max_thickness = airframe.foil.max_thickness()
    chord = airframe.c_r   
    # Compute twist deflection using torsion theory
    # symbolic solve dtheta_dz for skin and spar 
    T_spar = sp.Symbol('T_spar')  # Unknown torque carried by the skin
    A_skin, p_skin = airfoil_properties(airframe.foil, chord) # calculate area and perimeter of airfoil for skin torsion calculation
    skin_int_t_ds = p_skin/t_skin # 
    A_spar = np.pi*r_spar**2 # conservative estimation
    int_t_ds = np.pi*2*r_spar/r_thickness # calculate torsional constant by conservative approximation of min thickness all around, with assumption of 50% area efficiency (i.e. same area as square, twice the perimter)
    dtheta_dz = sp.Symbol('dtheta_dz')  # Common twist rate (not solved for)
    eq1 = sp.Eq(dtheta_dz, (Torsion_distribution[0] - T_spar)*skin_int_t_ds / (4 * PET.G * A_skin**2))  # Skin equation
    eq2 = sp.Eq(dtheta_dz,  T_spar*int_t_ds/ (4 * CFRP.G*A_spar**2))  # Spar equation
    solution = sp.solve([eq1, eq2], (T_spar,dtheta_dz))
    if len(solution) == 0:
        raise("No solution found for the system of equations.")
    compatibility_dtdz = solution[dtheta_dz]

    #reused variable name but can change this
    dtheta_dz = Torsion_distribution / (4*CFRP.G * A_spar**2)*int_t_ds
    compatibility_factor = compatibility_dtdz/dtheta_dz[0] # calculate compatibility factor by comparing dtheta_dz from spar and skin at root
    dtheta_dz=dtheta_dz*compatibility_factor
    dz = airframe.b/(points_loads*2*np.cos(airframe.qc_sweep))
    theta = np.cumsum(dtheta_dz) * dz # convert to degrees
    #plt.plot(np.linspace(0,len(theta)*dz,len(theta)), theta)
    #plt.ylim((0,len(theta)*dz))
    #plt.show()
    return theta[-1]#in degrees


def shear_force(Lift_distribution, airframe, t_spar, t_sleeve): #returns true if it passes this test #input drag works too
    safety_factor=5
    #root test stress/unit span
    sigma_root = safety_factor*Lift_distribution[0]/(2*t_spar) #lift stress at thinnest part of spar, with 5x safety factor
    sigma_sleeve = safety_factor*Lift_distribution[pos_first_connection(airframe)]/(2*t_sleeve) #lift stress at thinnest part of sleeve, with 5x safety factor
    root = Components_Materials.CFRP().sigma > sigma_root
    sleeve = Components_Materials.GLARE().sigma > sigma_sleeve
    return bool(sleeve*root)

def stress_analysis():              # Implement method to compute maximum stresses
    pass

def buckling_analysis():            # Analyze buckling performance
    pass

def deflection_analysis():          # Compute all relevant deflections
    pass

def dynamic_analysis():             # Check dynamic loads
    pass
