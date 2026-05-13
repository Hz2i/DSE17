import numpy as np
from matplotlib import pyplot as plt
from Objects.Characteristics.Airframe import wing, empennage, fuselage
class ScissorPlot:
    def __init__(self,wing=wing(),empennage=empennage(),fuselage=fuselage()):                                 # Initialise with required values (add inputs after self, as per necessity
        wing.compute_CL_grad()
        wing.compute_CL_max()
        wing.compute_Cm_ac()
        empennage.compute_CL_grad()
        self.MAC = wing.MAC
        self.Cm_ac = wing.Cm_ac
        self.CL_Ah = wing.CL_max
        self.CL_h = -0.35*(empennage.AR_h)**(1/3) # SEAD Lecture 8
        self.Sh = empennage.Sh
        self.lh = empennage.lh             # meter
        self.S = wing.S
        self.cbar = np.sqrt(wing.S/wing.AR)
        self.Vh2 = 0.85     # SEAD Lecture 7 slide 42

        self.CL_alphah = empennage.CL_grad_h
        self.CL_alpha = wing.CL_grad

        ## wing downwash gradient computation ##
        r = empennage.lh * 2 / wing.b
        m_tv = empennage.vh * 2 / wing.b
        K_eps_0 = 0.1124/r**2 + 0.1024/r + 2 # SEAD lecture 7
        K_eps = (0.1124 + 0.1265*wing.qc_sweep + 0.1766*wing.qc_sweep**2)/r**2 + 0.1024/r + 2
        self.depsilon_dalpha = K_eps/K_eps_0 * (r/(r**2 + m_tv**2) * 0.4876/np.sqrt(r**2+0.6319+m_tv**2) + (1+ (r**2/(r**2+0.7915 + 5.0734 * m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2)))) * self.CL_alpha / (np.pi*wing.AR) ## for testing, compare with result of 4/(AR+2)

        ## x_ac positioning computation ##
        x_ac_fus = -1.8/self.CL_alpha * fuselage.D**2 * (fuselage.L1 + 0.5) / (wing.S*wing.MAC) # destabilising contribution (Length of fuselage before wing start is given with a 0.5 m extra margin)
        x_ac_fus += 0.273/(1+wing.taper) * (fuselage.D * wing.geo_chord * (wing.b-fuselage.D))/(wing.MAC**2 * (wing.b + 2.15*fuselage.D))

        self.x_ac = wing.x_ac + x_ac_fus

    def stability_min_Sh_S(self,x_cg,margin):
        return (x_cg + margin - self.x_ac) / (self.CL_alphah/self.CL_alpha * (1-self.depsilon_dalpha) * self.lh / self.MAC * self.Vh2)

    def controllability_min_Sh_S(self,x_cg):
        return (x_cg - self.x_ac + self.Cm_ac/self.CL_Ah) / (self.CL_h/self.CL_Ah * self.lh/self.MAC * self.Vh2)

    def plot_scissor_plot(self,x_cg_min=0.2,x_cg_max=0.4):
        x_cg = np.arange(-0.1,0.8+0.1,0.1)
        
        y_control = self.controllability_min_Sh_S(x_cg)
        y_stab_margin = self.stability_min_Sh_S(x_cg,margin=0.05)
        y_stab = self.stability_min_Sh_S(x_cg,margin=0.0)

        Sh_S = self.Sh/self.S

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
