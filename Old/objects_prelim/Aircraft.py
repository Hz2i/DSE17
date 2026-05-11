import numpy as np

class aircraft:
    def __init__(self):                             # Initialise with required values (add inputs after self, as per necessity)
        self.volume = self.compute_volume()
        self.weight = self.compute_weight()
        self.cg = self.compute_cg() #From Nose
        self.Pr_tot, self.Pr_gen, self.Pr_stor = self.power_budget()
        self.Linkbudget = self.link_budget()
        self.velo = self.velocity()
        self.pos = self.position()
        self.aero = self.aero_properties()

    def compute_weight(self):                      # Compute weight of the aircraft
        pass

    def compute_volume(self):                      # Compute volume of the aircraft
        pass

    def compute_cg(self):                          # Compute CG of the aircraft
        pass
    
    def power_budget(self):                        # Compute power budget of the aircraft
        pass

    def link_budget(self):                         # Compute link budget of the aircraft
        pass

    def velocity(self):                            # Compute velocity of the aircraft
        pass

    def position(self):                            # Compute position of the aircraft
        pass

    def aero_properties(self):                       # Compute aerodynamic properties of the aircraft
        pass