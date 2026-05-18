import numpy as np
from matplotlib import pyplot as plt
from Objects.Characteristics.Airframe import wing, empennage, fuselage
class ScissorPlot:
    def __init__(self,wing=wing(),empennage=empennage(),fuselage=fuselage()):
        self.wing = wing
        self.empennage = empennage
        self.fuselage = fuselage                                 # Initialise with required values (add inputs after self, as per necessity
        self.Vh2 = 0.85     # SEAD Lecture 7 slide 42

    def compute_required_coefs(self):
        self.wing.compute_required_coefficients()
        self.wing.compute_CL_grad()
        self.wing.compute_CL_max()
        self.wing.compute_Cm_ac()
        self.empennage.compute_required_coefficients()
        self.empennage.compute_CL_grad()
        self.cbar = np.sqrt(self.wing.S/self.wing.AR)
        self.CL_h = -0.35*(self.empennage.AR_h)**(1/3) # SEAD Lecture 8

        ## wing downwash gradient computation ##
        r = self.empennage.lh * 2 / self.wing.b
        m_tv = self.empennage.vh * 2 / self.wing.b
        K_eps_0 = 0.1124/r**2 + 0.1024/r + 2 # SEAD lecture 7
        K_eps = (0.1124 + 0.1265*self.wing.qc_sweep + 0.1766*self.wing.qc_sweep**2)/r**2 + 0.1024/r + 2
        self.depsilon_dalpha = K_eps/K_eps_0 * (r/(r**2 + m_tv**2) * 0.4876/np.sqrt(r**2+0.6319+m_tv**2) + (1+ (r**2/(r**2+0.7915 + 5.0734 * m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2)))) * self.wing.CL_grad / (np.pi*self.wing.AR) ## for testing, compare with result of 4/(AR+2)

        ## x_ac positioning computation ##
        x_ac_fus = -1.8/self.wing.CL_grad * self.fuselage.D**2 * (self.fuselage.L1 + 0.5) / (self.wing.S*self.wing.MAC) # destabilising contribution (Length of fuselage before wing start is given with a 0.5 m extra margin)
        x_ac_fus += 0.273/(1+self.wing.taper) * (self.fuselage.D * self.wing.geo_chord * (self.wing.b-self.fuselage.D))/(self.wing.MAC**2 * (self.wing.b + 2.15*self.fuselage.D))

        self.x_ac = self.wing.x_ac + x_ac_fus

    def stability_min_Sh_S(self,x_cg,margin):
        return (x_cg + margin - self.x_ac) / (self.empennage.CL_grad_h/self.wing.CL_grad * (1-self.depsilon_dalpha) * self.empennage.lh / self.wing.MAC * self.Vh2)

    def controllability_min_Sh_S(self,x_cg):
        return (x_cg - self.x_ac + self.wing.Cm_ac/self.wing.CL_max) / (self.CL_h/self.wing.CL_max * self.empennage.lh/self.wing.MAC * self.Vh2)
    
    def minimum_Sh_S(self,x_cg_min=0.2,x_cg_max=0.4): ## iterate for CG and Sh size
        Sh_S_stability = self.stability_min_Sh_S(x_cg_max,margin=0.05)
        Sh_S_controllability = self.controllability_min_Sh_S(x_cg_min)

        Sh_S = np.maximum(Sh_S_stability,Sh_S_controllability)
        return Sh_S

    def plot_scissor_plot(self,x_cg_min=0.2,x_cg_max=0.4):
        x_cg = np.arange(-0.1,0.8+0.1,0.1)

        Sh_S = self.empennage.Sh/self.wing.S
        
        y_control = self.controllability_min_Sh_S(x_cg)
        y_stab_margin = self.stability_min_Sh_S(x_cg,margin=0.05)
        y_stab = self.stability_min_Sh_S(x_cg,margin=0.0)

        plt.hlines(Sh_S,x_cg_min,x_cg_max,colors=['black'],label="CG range")

        plt.plot(x_cg,y_control,color='blue',label="Controllability")
        plt.fill_between(x_cg,y_control,-1,color='mistyrose',alpha=1)
        plt.plot(x_cg,y_stab,color='red',label="Stability")
        plt.fill_between(x_cg,y_stab,-1,color='mistyrose',alpha=1)
        plt.plot(x_cg,y_stab_margin,color='green',label="Stability with margin")
        plt.hlines(0,x_cg[0],x_cg[-1],colors=['black'],linestyles=["dashed"])
        plt.ylim(-0.05,0.5)
        plt.xlabel(r"$x_{cg}/MAC$")
        plt.ylabel(r'$S_h / S $')
        plt.legend()
        plt.show()

scissor_plot = ScissorPlot()
scissor_plot.compute_required_coefs()
scissor_plot.plot_scissor_plot()