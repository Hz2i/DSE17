import numpy as np
from Objects.Characteristics.ReferenceGeometries import *


class wing:
    def __init__(self, S=36.0, A=25.0, qc_sweep=0.0, taper=1.0, dihedral=0.0 , airfoil=airfoil_e387()): # airfoil = foil()
        self.foil = airfoil
        self.AR = A
        self.qc_sweep = qc_sweep
        self.taper = taper
        self.dihedral = dihedral
        self.S = S                
        self.CL_grad = 0.0          # Currently initialised with 0; Add method to compute
        self.CL_max = 0.0           # Currently initialised with 0; Add method to compute
        self.e = 0.0
        self.x_ac = 0.25

        self.m = 0.0                # Currently initialised with 0; Class 2 estimation methods required!
        self.v_int = 0.0            # Currently initialised with 0; Stress calculations and internal design required!
        self.v_tot = 0.0            # Currently initialised with 0

    def compute_required_coefficients(self):
        self.c_sweep = np.arctan( np.tan( self.qc_sweep - 1/self.AR * (1-self.taper)/(1+self.taper) ) )
        self.geo_chord = np.sqrt(self.S/self.AR)
        self.root_chord = self.geo_chord * 2 / (1+self.taper)
        self.MAC = self.root_chord * 2/3 * ((1+self.taper+self.taper**2)/(1+self.taper))
        self.b = self.geo_chord*self.AR

    def compute_CL_grad(self):
        beta = 1.0
        eta = 0.95
        self.CL_grad = 2*np.pi*self.AR/(2+np.sqrt(4+(self.AR*beta/eta)**2 * (1+np.tan(self.c_sweep)/beta**2)))

    def compute_CL_max(self):       ## NOT COMPLETE ## cl_max needs to be replaced with effective cl_max of airfoil; If sweep =/= 0, scale needs to be adjusted according to graphs in formula sheet aircraft
        cl_max = self.foil.clmax
        scale = 0.9
        self.CL_max = scale*cl_max

    def compute_Cm_ac(self):
        self.Cm_ac = self.foil.cm_0 * (self.AR * np.cos(self.qc_sweep)/(self.AR + 2*np.cos(self.qc_sweep)))

    def compute_oswald_eff(self): # Paper: Estimating the Oswald Factor from Basic Aircraft Geometrical Parameters
        taper_adjusted = self.taper - (-0.357 + 0.45*np.exp(0.0375*self.qc_sweep))
        f_lambda = 0.0524*taper_adjusted**4 - 0.15*taper_adjusted**3 + 0.1659*taper_adjusted**2 - 0.0706*taper_adjusted + 0.0119
        k_WL = 2.83
        self.k_e_dihedral = (1+1/k_WL * (1/np.cos(self.dihedral)-1))**2
        self.e =  1/(1+f_lambda*self.AR)

    def zero_lift_drag(self,rho_cruise,V_cruise, M): # ADSEE 2 lectures; Requires cruise rho, cruise V, cruise M; Assumes average chord length based on surface area and AR.
        S_wet_w = 1.07*2*self.S
        chord = self.MAC
        mu = 1.4216e-5 # at 60,000 ft (18,500 m)
        k = 0.152e-5 # surface roughness of polished sheet metal

        Reynolds = np.minimum(rho_cruise*V_cruise*chord/mu,38.21*(chord/k)**1.053)

        Cf_l = 1.328/np.sqrt(Reynolds)
        Cf_t = 0.455/(np.log10(Reynolds)**2.58)

        f_l = 0.35

        Cf = f_l*Cf_l + (1-f_l)*Cf_t

        airfoil_x_maxthickness = self.foil.thickness_pos
        airfoil_thickness = self.foil.max_thickness

        FF = (1+0.6/airfoil_x_maxthickness * airfoil_thickness + 100 * airfoil_thickness**4)*(1.34*M**0.18*(np.cos(self.qc_sweep))**0.28)

        self.CD0 = Cf * FF * S_wet_w ## 1/S needs to be applied externally


'''

wing_properties = wing(taper=1.0,dihedral=np.pi/30)
wing_properties.compute_oswald_eff()
wing_properties.S = 40
print(wing_properties.e)
wing_properties.zero_lift_drag(0.07,25,0.1)
wing_properties.LD_computation(0.7)

Computed L/D ratio: 
    def LD_computation(self,CL):
        self.CL_CD = CL/(self.CD0 + CL**2 * 1/(np.pi*self.AR*self.e))
'''

class fuselage:
    def __init__(self,D=0.3,L1=0.5,L2=3,L3=1.5):
        self.D = D
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3

        self.m = 0.0            # Currently initialised with 0; Class 2 estimation methods required!
        self.v_total = 0.0      # Currently initialised with 0
    
    def zero_lift_drag(self, rho_cruise, V_cruise):
        self.Sw = np.pi*self.D/4 * (1/(3*self.L1**2)*((4*self.L1**2+self.D**2/4)**1.5-self.D**3/8)-self.D+4*self.L2+2*np.sqrt(self.L3**2+self.D**2/4))

        length = self.L1+self.L2+self.L3
        mu = 1.4216e-5 # at 60,000 ft (18,500 m)
        k = 0.152e-5 # surface roughness of polished sheet metal

        Reynolds = np.minimum(rho_cruise*V_cruise*length/mu,38.21*(length/k)**1.053)

        Cf_l = 1.328/np.sqrt(Reynolds)
        Cf_t = 0.455/(np.log10(Reynolds)**2.58)

        f_l = 0.1

        Cf = f_l*Cf_l + (1-f_l)*Cf_t

        if self.D == 0:
            f = 1
        else:
            f = (length)/self.D

        self.CD0 = self.Sw * (Cf*(1+60/f**3+f/400)) ## 1/S needs to be applied externally



class empennage:
    def __init__(self,S_h = 6, S_v = 4, h_AR = 4.0, v_AR = 1.5, lh = 3, vh = 0.0, taper_h = 1.0, taper_v = 1.0, airfoilh=airfoil_NACA0012(),airfoilv=airfoil_NACA0012(), qcsweep_h = 0.0, qcsweep_v = 0.0 ): # average values taken from Roelof Vos' ADSEE book
        self.foil_h = airfoilh
        self.foil_v = airfoilv
        self.AR_h = h_AR         # Currently initialised with 0
        self.AR_v = v_AR         # Currently initialised with 0
        self.Sh = S_h           # Currently initialised with 0
        self.Sv = S_v           # Currently initialised with 0
        self.qc_sweep_h = qcsweep_h   # Currently initialised with 0
        self.qc_sweep_v = qcsweep_v   # Currently initialised with 0
        self.taper_h = taper_h
        self.taper_v = taper_v

        self.m = 0.0            # Currently initialised with 0
        self.lh = lh            # Longitudinal distance between x_ac_w and x_ac_h
        self.vh = vh            # Vertical distance between x_ac_w and x_ac_h

    def compute_required_coefficients(self):
        self.geo_chord_h = np.sqrt(self.Sh/self.AR_h)
        self.root_chord_h = self.geo_chord_h * 2 / (1+self.taper_h)
        self.MAC_h = self.root_chord_h * 2/3 * ((1+self.taper_h+self.taper_h**2)/(1+self.taper_h))

        self.geo_chord_v = np.sqrt(self.Sv/self.AR_v)
        self.root_chord_v = self.geo_chord_v * 2 / (1+self.taper_v)
        self.MAC_v = self.root_chord_v * 2/3 * ((1+self.taper_v+self.taper_v**2)/(1+self.taper_v))

        self.c_sweep_h = np.arctan( np.tan( self.qc_sweep_h - 1/self.AR_h * (1-self.taper_h)/(1+self.taper_h) ) )

    def compute_CL_grad(self):
        beta = 1.0
        eta = 0.95
        self.CL_grad_h = 2*np.pi*self.AR_h/(2+np.sqrt(4+(self.AR_h*beta/eta)**2 * (1+np.tan(self.c_sweep_h)/beta**2)))

    def zero_lift_drag(self,rho_cruise,V_cruise, M): # ADSEE 2 lectures; Requires cruise rho, cruise V, cruise M; Assumes average chord length based on surface area and AR.
        S_wet_h = 1.05*2*self.Sh
        S_wet_v = 1.05*2*self.Sv
        
        chord_h = self.MAC_h
        chord_v = self.MAC_v
        mu = 1.4216e-5 # at 60,000 ft (18,500 m)
        k = 0.152e-5 # surface roughness of polished sheet metal

        Reynolds_h = np.minimum(rho_cruise*V_cruise*chord_h/mu,38.21*(chord_h/k)**1.053)
        Reynolds_v = np.minimum(rho_cruise*V_cruise*chord_v/mu,38.21*(chord_v/k)**1.053)

        if Reynolds_h <= 0:
            Reynolds_h = 10
        if Reynolds_v <= 0:
            Reynolds_v = 10

        Cf_l_h = 1.328/np.sqrt(Reynolds_h)
        Cf_t_h = 0.455/(np.log10(Reynolds_h)**2.58)

        Cf_l_v = 1.328/np.sqrt(Reynolds_v)
        Cf_t_v = 0.455/(np.log10(Reynolds_v)**2.58)

        f_l = 0.35 # laminar fraction

        Cf_h = f_l*Cf_l_h + (1-f_l)*Cf_t_h # coefficient of friction
        Cf_v = f_l*Cf_l_v + (1-f_l)*Cf_t_v # coefficient of friction

        FF_h = (1+0.6/self.foil_h.thickness_pos * self.foil_h.max_thickness + 100 * self.foil_h.max_thickness**4)*(1.34*M**0.18*(np.cos(self.qc_sweep_h))**0.28)
        FF_v = (1+0.6/self.foil_v.thickness_pos * self.foil_v.max_thickness + 100 * self.foil_v.max_thickness**4)*(1.34*M**0.18*(np.cos(self.qc_sweep_v))**0.28)

        self.CD0 = (Cf_h * FF_h * S_wet_h + Cf_v * FF_v * S_wet_v) ## 1/S needs to be applied externally


class nacelles:
    def __init__(self,nr_of_engines = 4,diameter = 0.5, length = 1):
        self.nr_of_engines = nr_of_engines
        self.D = diameter
        self.L = length

        
    def zero_lift_drag(self, rho_cruise, V_cruise):
        self.Sw = np.pi*self.D/4 * (1/(3*self.L**2)*((4*self.L**2+self.D**2/4)**1.5))


        mu = 1.4216e-5 # at 60,000 ft (18,500 m)
        k = 0.152e-5 # surface roughness of polished sheet metal

        Reynolds = np.minimum(rho_cruise*V_cruise*self.L/mu,38.21*(self.L/k)**1.053)

        Cf_l = 1.328/np.sqrt(Reynolds)
        Cf_t = 0.455/(np.log10(Reynolds)**2.58)

        f_l = 0.1

        Cf = f_l*Cf_l + (1-f_l)*Cf_t

        if self.D == 0:
            f = 1
        else:
            f = (self.L)/self.D

        self.CD0 = self.nr_of_engines*self.Sw * (Cf*(1+0.35/f)) ## 1/S needs to be applied externally

