import sys
import os
import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

import Airframe

#Find n sections and slanted span start of first section and section length

class Sections:
    def __init__(self, airframe=Airframe.airframe(qc_sweep=np.radians(15)), Plot=False):
        #Initial Guess
        self.main_section_inside = 0.5
        self.airframe = airframe
        self.max_length = 5.6#From Container #5.6
        self.slanted_span_b_2 = (self.airframe.b/2) / np.cos(self.airframe.qc_sweep)
        self.n_sections = 4
        self.optimized_main_section_inside,_ = self.optimize_sections()
        self.l_main_section, self.l_main_section_spar, self.l_wingtip_section, self.l_total, self.l_middle_wing_section = self.find_sections(self.optimized_main_section_inside)
        if Plot:
            self.plot_sections(self.optimized_main_section_inside)

        

    def find_sections(self, main_section_inside):
        n_sections = 4
        main_section_length = (self.max_length/2) / np.cos(self.airframe.qc_sweep)
        main_section_spar = main_section_length - main_section_inside
        wingtip_section_length = self.slanted_span_b_2 - main_section_inside - (n_sections-1)*main_section_spar*2
        total_length = main_section_inside + wingtip_section_length + (n_sections-1)*main_section_spar*2
        middle_wing_section = main_section_spar*2
        return main_section_length, main_section_spar, wingtip_section_length, total_length, middle_wing_section
    
    def optimize_sections(self):
        opti = asb.Opti()
        
        # Define optimization variables
        main_section_inside = opti.variable(init_guess=self.main_section_inside, lower_bound=0.5, upper_bound=(self.max_length/2-1))
        #Optimize Function
        function = self.find_sections(main_section_inside)
        #Feasibility
        opti.subject_to(function[1] <= (self.max_length/2))
        # opti.subject_to((function[1]) >= (self.max_length/2))
        
        
        opti.maximize(function[1])
        # Define the objective function (e.g., minimize the difference between main section length and wingtip section length)
        sol = opti.solve(verbose=False)
        return sol.value(main_section_inside), self.find_sections(sol.value(main_section_inside))
    
    def plot_sections(self, main_section_inside):
        main_section_length, main_section_spar, wingtip_section_length, total_length, middle_wing_section = self.find_sections(main_section_inside)

        section_names = [
            "Main Inside",
            "Spar 1",
            "Spar 2",
            "Spar 3",
            "Wingtip"
        ]

        section_lengths = [
            main_section_inside,
            middle_wing_section,
            middle_wing_section,
            middle_wing_section,
            wingtip_section_length
        ]

        print(f"Main Section Length: {main_section_length:.3f} m")
        print(f"Main Section Spar: {main_section_spar:.3f} m")
        print(f"Wingtip Section Length: {wingtip_section_length:.3f} m")
        print(f"Total Length: {total_length:.3f} m")
        print(f"Middle Wing Section Length: {middle_wing_section:.3f} m")

        plt.figure(figsize=(8, 4))
        bars = plt.bar(section_names, section_lengths)

        for bar, length in zip(bars, section_lengths):
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{length:.2f}",
                ha="center",
                va="bottom"
            )

        plt.ylabel("Length [m]")
        plt.xlabel("Wing Section")
        plt.title("Wing Section Lengths Along Slanted Span")
        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    sections = Sections(Plot=True, airframe=Airframe.airframe(qc_sweep=np.radians(15), S=69.2, A=20, init_polar=False))