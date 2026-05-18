import numpy as np
import ambiance as am
import matplotlib.pyplot as plt
import Objects.Characteristics.PowerSystem_sizing as ps
import Objects.Characteristics.PropulsionSystem as prop

class MissionProfile:
    def __init__(self, guesslist):  
        #Guesses
        self.m_solar_guess = guesslist[0] #kg
        self.m_battery_guess = guesslist[1] #kg
        self.m_prop_guess = guesslist[2] #kg
        self.gamma_guess = guesslist[3] #deg
        self.S_guess = guesslist[4] #m^2

        #Given
        self.mass_subsys = 46.8 #kg
        self.g = 9.81                                     
        self.alt = 60000*0.3048
        self.CD0 = 0.010
        self.AR = 24
        self.e = 0.9
        self.Pavg_climb_subsys = 300 #W
        self.Pavg_cruise_subsys = 425 #W
        self.eta_prop = 0.8
        self.LD = 40
        self.V_cruise = 25 #m/s
        self.t_night = 827*60 #s
        self.specific_power = 400 #Wh/kg Battery specific power
        self.specific_mass_solar = 4 #kg/m^2

        #Extra Parameters
        self.m_total_guess = self.m_solar_guess + self.m_battery_guess + self.m_prop_guess + self.mass_subsys
        self.density_climb = self.Calc_Density_Climb()
        self.Cl_opt_climb = self.Calc_Cl_opt_climb()
        self.CD_total_climb = self.Calc_CD_total(self.Cl_opt_climb)
        self.gamma_rad = np.radians(self.gamma_guess)

        
        #Climb
        self.V_climb = self.Calc_V_climb()
        self.t_climb = self.Calc_t_climb()
        self.D_climb = self.Calc_D_climb()
        self.E_climb = self.Calc_E_climb()
        self.ROC = self.Calc_ROC()

        #Cruise
        self.Pprop_cruise = self.Calc_Pprop_cruise()
        self.Pavg_cruise = self.Pprop_cruise + self.Pavg_cruise_subsys
        self.E_cruise = self.Calc_E_cruise()

        #Solar Size
        self.S_solar = Calc_S_solar()


    def Calc_Density_Climb(self):
        h_range = np.linspace(0, self.alt, 100000)
        density = am.Atmosphere(h_range).density
        #integrate density to get mass of air in the climb envelope
        dens_climb = np.trapz(density, h_range)/(self.alt)
        return dens_climb

    def Calc_Cl_opt_climb(self):
        return np.sqrt(self.CD0 * np.pi * self.AR * self.e)
    
    def Calc_CD_total(self, CL):
        k = 1.0/(np.pi * self.AR * self.e)
        return self.CD0 + k * CL**2
    
    def Calc_V_climb(self):
        V_climb = np.sqrt(2*self.m_total_guess*self.g*np.cos(self.gamma_rad)/(self.density_climb*self.S_guess*self.Cl_opt_climb))
        return V_climb
    
    def Calc_D_climb(self):
        D_climb = 0.5*self.CD_total_climb*self.S_guess*self.density_climb*self.V_climb**2
        return D_climb
    
    def Calc_t_climb(self):
        t_climb = self.alt / (self.V_climb * np.sin(self.gamma_rad))
        return t_climb

    def Calc_E_climb(self):
        Epot = self.m_total_guess * self.g * self.alt
        Edrag = self.D_climb * (self.alt / np.sin(self.gamma_rad))
        Esubsys_climb = self.Pavg_climb_subsys * self.t_climb
        Eclimb = Epot + Edrag + Esubsys_climb
        return Eclimb
    
    def Calc_ROC(self):
        ROC = self.V_climb * np.sin(self.gamma_rad)
        return ROC

    def Calc_Pprop_cruise(self):
        Pprop_cruise = ((self.m_total_guess*self.g)/(self.LD))*self.V_cruise/self.eta_prop
        return Pprop_cruise

    def Calc_E_cruise(self):
        Ecruise = (self.Pavg_cruise) * self.t_night
        return Ecruise
    
    def Calc_S_solar(self):
        S_solar =  
class MassUpdate:

    
if __name__ == "__main__":
    #Guesses
    m_solar_guess = 9.2 #kg
    m_battery_guess = 80 #kg
    m_prop_guess = 10 #kg
    gamma_guess = 5 #deg
    S_guess = 36 #m^2
    guesslist = [m_solar_guess, m_battery_guess, m_prop_guess, gamma_guess, S_guess]

    #Optimization loop
    target_ROC = float(input("Target ROC: "))
    roc_tolerance = 0.01
    guess_tolerance = 1e-3
    max_iterations = 500
    step = 0.1

    initial_guess = DeploybilityWorstCase(guesslist=guesslist)
    previous_guesslist = guesslist.copy()
    iteration = 0

    while iteration < max_iterations:
        roc_error = target_ROC - initial_guess.ROC
        guess_delta = np.max(np.abs(np.array(guesslist) - np.array(previous_guesslist)))

        if iteration > 0 and abs(roc_error) <= roc_tolerance and guess_delta <= guess_tolerance:
            break

        previous_guesslist = guesslist.copy()

        direction = np.sign(roc_error)
        adjustment = max(step, 0.05 * abs(roc_error))

        # Reduce total mass and wing loading when ROC is too low, and relax them when ROC is too high.
        guesslist[0] = max(0.1, guesslist[0] - direction * adjustment)
        guesslist[1] = max(0.1, guesslist[1] - direction * adjustment)
        guesslist[2] = max(0.1, guesslist[2] - direction * adjustment)
        guesslist[3] = np.clip(guesslist[3] + direction * adjustment, 1.0, 20.0)
        guesslist[4] = max(1.0, guesslist[4] - direction * adjustment)

        initial_guess = DeploybilityWorstCase(guesslist=guesslist)
        iteration += 1
        print(
            f"Iteration {iteration}: ROC = {initial_guess.ROC:.3f} m/s, "
            f"Masses = {guesslist[0]:.2f}, {guesslist[1]:.2f}, {guesslist[2]:.2f} kg, "
            f"gamma = {guesslist[3]:.2f} deg, S = {guesslist[4]:.2f} m^2"
        )

    if iteration >= max_iterations:
        print("Warning: optimizer reached max iterations before converging.")

    print(f"Final ROC: {initial_guess.ROC:.3f} m/s")
    print(f"Final guesses: {guesslist}")