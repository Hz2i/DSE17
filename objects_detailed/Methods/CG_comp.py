import sys
import os
import aerosandbox.numpy as np
import plotly.graph_objects as go
import aerosandbox as asb
import StructuralAnalysis as sa

import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

import Airframe
import Components_Materials
import GeneralSubsystems
import PowerSystem_sizing
import SparGeometryParam 

class Available_BatteryCrossSection:
    def __init__(self, AirGEO = SparGeometryParam.AirfoilGeometry(), width_clamp=0.125, plot = False):
        self.airfoil_geometry = AirGEO
        self.width_clamp = width_clamp
        self.battery_cross_section_area, self.x_centroid_batt, self.z_centroid_batt = self.Calculate_Battery_CrossSection_Area()
        if plot:
            self.plot_airfoil_cross_section()

    def Calculate_Battery_CrossSection_Area(self):
        self.x_wall_battery = self.airfoil_geometry.x_max_thickness - self.width_clamp/2
        #Find Area Before X_wall
        x_u = self.airfoil_geometry.xu_off
        y_u = self.airfoil_geometry.yu_off
        x_l = self.airfoil_geometry.xl_off
        y_l = self.airfoil_geometry.yl_off
        
        x_wall = self.x_wall_battery
        #Reverse order of x_u and y_u to ensure correct ordering for interpolation and area calculation
        x_u = x_u[::-1]
        y_u = y_u[::-1]

        # Find insertion index
        idx_u = np.searchsorted(x_u, x_wall)
        idx_l = np.searchsorted(x_l, x_wall)

        # Interpolate y-values at wall
        y_wall_u = np.interp(x_wall, x_u, y_u)
        y_wall_l = np.interp(x_wall, x_l, y_l)

        

        # Build truncated arrays including wall point
        self.x_u_cut = np.concatenate([x_u[:idx_u], [x_wall]])
        self.y_u_cut = np.concatenate([y_u[:idx_u], [y_wall_u]])

        self.x_l_cut = np.concatenate([x_l[:idx_l], [x_wall]])
        self.y_l_cut = np.concatenate([y_l[:idx_l], [y_wall_l]])

        #interpolate
        x = np.linspace(0, self.x_wall_battery, 200)

        interp_upper = np.interp(x, x_u, y_u)
        interp_lower = np.interp(x, x_l, y_l)

        height = interp_upper - interp_lower

        # Area
        area = np.trapezoid(height, x)
        x_centroid_batt = np.trapezoid(x * height, x) / area
        z_centroid_batt = np.trapezoid(0.5*(interp_upper + interp_lower) * height, x) / area
        return area, x_centroid_batt, z_centroid_batt         

    def plot_airfoil_cross_section(self):
        plt.figure(figsize=(10, 5))
        plt.scatter(self.x_u_cut, self.y_u_cut, color='blue', label="Upper Surface (Truncated)", s=10)
        plt.scatter(self.x_l_cut, self.y_l_cut, color='red', label="Lower Surface (Truncated)", s=10)
        plt.scatter(self.x_centroid_batt, self.z_centroid_batt, color='black', label="Battery Centroid", s=50)
        plt.plot(self.airfoil_geometry.xu_off, self.airfoil_geometry.yu_off, color='blue', label="Upper Surface Offset", linestyle='--')
        plt.plot(self.airfoil_geometry.xl_off, self.airfoil_geometry.yl_off, color='red', label="Lower Surface Offset", linestyle='--')
        plt.vlines(self.airfoil_geometry.x_max_thickness, ymin=-0.5, ymax=0.5, color='orange', label="Max Thickness", linestyle='-.')
        plt.vlines(self.x_wall_battery, ymin=-0.5, ymax=0.5, color='green', label="Battery Cross-Section Limit", linestyle='-.')
        plt.title("Airfoil Cross-Section with Available Region")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.legend()
        plt.axis('equal')
        plt.grid()
        plt.show()

class CG_comp:
    def __init__(self, x_cg_goal=2.4, airframe=Airframe.airframe(),
                #  Sol=PowerSystem_sizing.power_generation(),
                #  Batt=PowerSystem_sizing.power_storage(),
                 Fcs = GeneralSubsystems.FlightConditionsSystem(),
                 Computer = GeneralSubsystems.ComputerSystem(),
                 Comms = GeneralSubsystems.CommunicationSystem(), 
                 Payload = GeneralSubsystems.PayloadSystem(),
                 Batt = Available_BatteryCrossSection(AirGEO = SparGeometryParam.AirfoilGeometry(), width_clamp=0.125),
                 t_skin_airfoil=0.0002):
        
        #section
        self.l_main_section = 0.5
        self.l_main_section_spar = 2.28 - self.l_main_section
        #Guesses
        self.l_batt_dist_guess = self.l_main_section_spar*2
        print(f"Initial Battery Distance Guess: {self.l_batt_dist_guess:.3f} m")
        self.Batt = Batt
        self.airframe = airframe
        self.airfoil = airframe.foil
        self.x_cg_goal = x_cg_goal
        self.panel_size = 0.01
        self.t_skin_airfoil = t_skin_airfoil
        # self.plot_mesh_plotly()
        #Centroid of Airfoil
        self.x_centroid, self.z_centroid = self.Calculate_Airfoil_Centroid()
        '''Constraint Geometry'''
        self.x_cg_skin()
        self.x_cg_wingtips()
        self.x_cg_solar_panel()
        # self.x_cg_spar()
        '''Changable Positions'''
        #Masses All below 1kg is omitted
        self.mass_payload = Payload.mass_payload + Payload.mass_mounting
        self.mass_fcs_pitot = Fcs.pitot_mass * Fcs.redundancy
        self.mass_computer = 0 #No components above 1kg, so not included
        self.mass_comms_ssr = Comms.ssr_mass*Comms.redundancy
        self.mass_comms_elt = Comms.elt_mass*Comms.redundancy
        self.mass_battery_1 = 0.104 #[kg]
        self.volume_battery_tot = 0.0368 #[m^3]
        self.volume_battery_1 = 0.140*0.054*0.0065
        self.prop_eng = 10 #[kg] Combining 2!!
        self.prop_sticklength = 0.5 #[m]
        #Positions
        self.x_payload = self.x_centroid
        self.x_fcs_pitot = self.x_centroid
        self.x_computer = self.x_centroid
        self.x_comms_ssr = self.x_centroid
        self.x_comms_elt = self.x_centroid
        self.x_prop_2_3, self.x_prop_1_4 = self.calc_x_prop_positions()

        self.y_payload = 0
        self.y_fcs_pitot = 0
        self.y_computer = 0
        self.y_comms_ssr = 0
        self.y_comms_elt = 0
        self.y_prop2_3 = self.airframe.b/6
        self.y_prop1_4 = self.airframe.b/6*2
        #Make CGlists
        self.mass_list = self.gen_mass_list()
        self.x_list = self.gen_x_list()
        self.y_list = self.gen_y_list()
        self.name_list = self.gen_name_list()

        #Plot CG
        print(f"Battery x: {self.Batt_Distribution(self.l_batt_dist_guess)[0]}")
        print(f"Battery y: {self.Batt_Distribution(self.l_batt_dist_guess)[1]}")
        print(f"Battery areafilled: {self.Batt_Distribution(self.l_batt_dist_guess)[2]}")
        print(f"Battery total mass: {self.Batt_Distribution(self.l_batt_dist_guess)[3]}")
        self.plot_CG()


    def calc_x_prop_positions(self):
        #Equally spaced along span.
        d_eng = self.airframe.b/6
        x_eng_2_3 = d_eng*np.tan(self.airframe.le_sweep) - self.prop_sticklength
        x_eng_1_4 = x_eng_2_3 + d_eng*np.tan(self.airframe.le_sweep)
        return x_eng_2_3, x_eng_1_4

    def Batt_Distribution(self, l_batt_dist):
        n_batteries = 10*int(np.ceil((self.volume_battery_tot / self.volume_battery_1)/10))
        batt_total_mass = n_batteries * self.mass_battery_1
        batt_total_vol_half= n_batteries * self.volume_battery_1/2
        #Find Datum Point
        x_datum_bat = (self.l_main_section + 3.56*2)*np.sin(self.airframe.le_sweep)
        y_datum_bat = (self.l_main_section + 3.56*2)*np.cos(self.airframe.le_sweep) 
        print(f"Battery Datum Point: ({x_datum_bat:.3f}, {y_datum_bat:.3f}) m")
        half_batt_l = l_batt_dist/2
        print(f"Battery Distribution: {n_batteries} batteries, total mass {batt_total_mass:.2f} kg, total half volume {batt_total_vol_half:.6f} m^3, distributed over {l_batt_dist:.3f} m")
        batt_centroid_x = x_datum_bat + self.Batt.x_centroid_batt + half_batt_l*np.sin(self.airframe.le_sweep)
        batt_centroid_y = y_datum_bat + half_batt_l*np.cos(self.airframe.le_sweep)
        batt_cross_section_filled = batt_total_vol_half/l_batt_dist
        returnlist = np.asarray([batt_centroid_x, batt_centroid_y, batt_cross_section_filled, batt_total_mass]).astype(float).tolist()
        return returnlist

    def Calculate_Airfoil_Centroid(self):
        # Datum: leading edge, x along chord, z upward

        upper_coords = self.airfoil.upper_coordinates()  # (x, z)
        lower_coords = self.airfoil.lower_coordinates()  # (x, z)

        # nondimensional x grid
        x = np.linspace(0, 1, 200)

        # interpolate surfaces (make sure ordering is correct)
        z_upper = np.interp(x, upper_coords[:, 0][::-1], upper_coords[:, 1][::-1])
        z_lower = np.interp(x, lower_coords[:, 0], lower_coords[:, 1])

        thickness = z_upper - z_lower

        # differential area
        dA = thickness

        # total area (nondimensional)
        A = np.trapezoid(dA, x)

        # centroid x (nondimensional)
        x_centroid = np.trapezoid(x * dA, x) / A

        # centroid z (nondimensional)
        z_mid = 0.5 * (z_upper + z_lower)
        z_centroid = np.trapezoid(z_mid * dA, x) / A

        # scale to physical coordinates
        x_centroid *= self.airframe.c_r
        z_centroid *= self.airframe.c_r   # assumes geometric similarity scaling
        return x_centroid, z_centroid
    
    '''Constraint Geometry'''
    def x_cg_skin(self):
        area_out, perimeter_out = sa.airfoil_properties(self.airfoil, chord_length=self.airframe.c_r)
        area_skin = perimeter_out * self.t_skin_airfoil  # Approximate skin area as perimeter times thickness, which is reasonable for thin skins around complex shapes like airfoils
        #Account for sweep and dihedral in z-centroid by projecting it along the local normal direction of the airfoil surface, which is tilted due to sweep and dihedral
        self.x_skin = ((self.airframe.b/2)*np.tan(self.airframe.le_sweep))/2 + self.x_centroid  # Shift x-centroid by half the span times tan(sweep) to account for sweep
        self.y_skin = ((self.airframe.b/2)/2)
        self.mass_skin = area_skin * Components_Materials.Mylar().rho*self.airframe.b/np.cos(self.airframe.le_sweep)  # Mass = area * density, scaled by half the span to account for distribution along the wing
    
    def x_cg_wingtips(self):
        #Assume Point For on Centroid of Edge Airfoil
        self.x_wingtips = ((self.airframe.b/2)*np.tan(self.airframe.le_sweep))+self.x_centroid
        self.y_wingtips = self.airframe.b/2  
        self.mass_wingtips = 2*15.0 #[kg] ()

    def x_cg_solar_panel(self):
        self.solar_panel_area = 36
        self.mass_solar_panel = self.solar_panel_area * Components_Materials.solar_panel().surfRho
        self.x_solar_panel = ((self.airframe.b/2)*np.tan(self.airframe.le_sweep))/2 + self.x_centroid
        self.y_solar_panel = ((self.airframe.b/2)/2)
    # def x_cg_spar(self):
    #     #Clamps + Rubber
    #     #Sleeve
    #     #Spar
    #     return x_cg_spar, mass_spar
    
    def gen_mass_list(self):
        mass_list = [
            self.mass_payload,
            self.mass_fcs_pitot,
            self.mass_computer,
            self.mass_comms_ssr,
            self.mass_comms_elt,
            self.mass_solar_panel,
            self.prop_eng,
            self.prop_eng,
            self.mass_wingtips,
            self.mass_skin
            #self.mass_spar
        ]
        return mass_list
    
    def gen_x_list(self):
        x_list = [
            self.x_payload,
            self.x_fcs_pitot,
            self.x_computer,
            self.x_comms_ssr,
            self.x_comms_elt,
            self.x_solar_panel,
            self.x_prop_2_3, #Assuming propellers are distributed around their CGs
            self.x_prop_1_4,
            self.x_wingtips, #Assuming wingtips are distributed around their CG
            self.x_skin #Assuming skin mass is distributed around its centroid
            #self.x_cg_spar
        ]
        return x_list
    
    def gen_y_list(self):
        y_list = [
            self.y_payload,
            self.y_fcs_pitot,
            self.y_computer,
            self.y_comms_ssr,
            self.y_comms_elt,
            self.y_solar_panel,
            self.y_prop2_3,
            self.y_prop1_4,
            self.y_wingtips,
            self.y_skin
            #self.y_cg_spar
        ]
        return y_list
    
    def gen_name_list(self):
        name_list = [
            "Payload",
            "FCS Pitot",
            "Computer",
            "Comms SSR",
            "Comms ELT",
            "Solar Panel",
            "Propeller 2/3",
            "Propeller 1/4",
            "Wingtips",
            "Skin",
            # "Spar"
        ]
        return name_list

    '''Calculate CG'''
    # def calculate_cg(self):
    #     mass_list = self.mass_list
    #     x_list = self.x_list
    #     y_list = self.y_list
    #     x_cg = sum(m * x for m, x in zip(mass_list, x_list)) / sum(mass_list)
    #     y_cg = sum(m * y for m, y in zip(mass_list, y_list)) / sum(mass_list)
    #     return x_cg, y_cg

    def calculate_cg_batt(self, l_batt_dist):

        batt_x, batt_y, batt_area, batt_mass = self.Batt_Distribution(l_batt_dist)

        mass_list = np.append(np.asarray(self.gen_mass_list(), dtype=float), batt_mass)
        x_list = np.append(np.asarray(self.gen_x_list(), dtype=float), batt_x)
        y_list = np.append(np.asarray(self.gen_y_list(), dtype=float), batt_y)

        name_list = self.gen_name_list() + ["Battery"]

        # CG calculation (fully vectorized)
        x_cg = np.sum(mass_list * x_list) / np.sum(mass_list)
        y_cg = np.sum(mass_list * y_list) / np.sum(mass_list)
        print(name_list)

        return mass_list, x_list, y_list, name_list, float(x_cg), float(batt_area), float(l_batt_dist) # Return CG and battery parameters for optimization feedback
    

    '''Optimizer'''

    # def Optimizer_CG():
    #     opti = asb.Opti()
       
    #     # Define optimization variables
    #     l_batt_dist = opti.variable(init_guess=(self.main_section_spar*2), lower_bound=0.5, upper_bound=(self.max_length/2-1))
    #     #Optimize Function
    #     function = self.find_sections(main_section_inside)
    #     #Feasibility
    #     opti.subject_to(function[1] <= (self.max_length/2))
    #     # opti.subject_to((function[1]) >= (self.max_length/2))
        
        
    #     opti.maximize(function[1])
    #     # Define the objective function (e.g., minimize the difference between main section length and wingtip section length)
    #     sol = opti.solve()
    #     return sol.value(main_section_inside), self.find_sections(sol.value(main_section_inside))
    def plot_wing_contour(self):
        cr = self.airframe.c_r
        b = self.airframe.b
        sweep = self.airframe.le_sweep  # radians

        y = np.array([0, b/2])

        x_le = np.array([
            0,
            (b/2) * np.tan(sweep)
        ])

        x_te = x_le + cr

        x_wing = np.array([
            x_le[0], x_le[1],
            x_te[1], x_te[0],
            x_le[0]
        ])

        y_wing = np.array([
            y[0], y[1],
            y[1], y[0],
            y[0]
        ])

        return go.Scatter(
            x=x_wing,
            y=y_wing,
            mode="lines",
            name="Wing Contour",
            line=dict(color="black", width=2)
        )

    def plot_CG(self):
        import plotly.io as pio
        pio.renderers.default = "browser"
        fig = go.Figure()
        fig.update_layout(
            xaxis=dict(range=[-1, 6]),
            yaxis=dict(range=[-1, 20])
        )
        mass_list, x_list, y_list, name_list, x_cg, batt_area, l_batt_dist = self.calculate_cg_batt(self.l_batt_dist_guess)
        # CG points
        for i in range(len(mass_list)):
            fig.add_trace(
                go.Scatter(
                    x=[x_list[i]],
                    y=[y_list[i]],
                    mode="markers",
                    name=f"{name_list[i]} (m={mass_list[i]:.2f})",
                    marker=dict(
                        size=mass_list[i]*5,
                        opacity=0.7
                    )
                )
            )

        # Target CG
        fig.add_trace(
            go.Scatter(
                x=[self.x_cg_goal],
                y=[0],
                mode="markers",
                name="Target CG",
                marker=dict(symbol="x", size=12, color="black")
            )
        )
        # Computed CG
        fig.add_trace(
            go.Scatter(
                x=[x_cg],
                y=[0],
                mode="markers",
                name="Computed CG",
                marker=dict(symbol="x", size=24, color="red")
            )
        )

        # Wing contour
        fig.add_trace(self.plot_wing_contour())

        fig.update_layout(
            title="Center of Gravity Distribution",
            xaxis_title="X Coordinate",
            yaxis_title="Y Coordinate",
            showlegend=True
        )

        fig.update_layout(
            xaxis=dict(range=[-1, 6]),
            yaxis=dict(range=[-1, 20])
        )
        fig.write_html("cg_plot.html", auto_open=True)
        # fig.show()

class Fitting:
    def __init__(self, x_data, y_data):
        self.x_data = x_data
        self.y_data = y_data

    def linear_fit(self):
        coeffs = np.polyfit(self.x_data, self.y_data, 1)
        return coeffs
    
    '''Check Fitting'''
    def AC_layout(self):
        b = self.airframe.b
        c_r = self.airframe.c_r
        c_t = self.airframe.c_t
        sweep = self.airframe.le_sweep

        edge_margin = 0.005
        span = b / 2

        # -------------------------
        # Panel grid (centers)
        # -------------------------
        y_edges = np.arange(-span + edge_margin,
                             span - edge_margin + self.panel_size,
                             self.panel_size)

        x_edges = np.arange(edge_margin,
                            c_r + edge_margin,
                            self.panel_size)

        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        x_centers = (x_edges[:-1] + x_edges[1:]) / 2

        Xc, Yc = np.meshgrid(x_centers, y_centers)
        Xc = Xc + np.abs(Yc) * np.tan(sweep)  # Shift x-coordinates based on sweep

        # -------------------------
        # Wing geometry (FIXED SWEEP)
        # -------------------------
        chord = c_r - (np.abs(Yc) / span) * (c_r - c_t)
        chord = np.maximum(chord, 1e-6)

        x_le = np.abs(Yc) * np.tan(sweep)  

        # -------------------------
        # Airfoil thickness
        # -------------------------
        x_over_c = (Xc - x_le) / chord
        x_over_c = np.clip(x_over_c, 0.0, 1.0)

        t_over_c = self.airfoil.local_thickness(x_over_c)

        Z = t_over_c * chord  # meters

        print("Xc shape: ", Xc.shape)
        print("Yc shape: ", Yc.shape)
        print("Z shape: ", Z.shape)
        return Xc, Yc, Z

    def plot_mesh_plotly(self):

        X, Y, Z = self.AC_layout()

        b = self.airframe.b
        c_r = self.airframe.c_r
        c_t = self.airframe.c_t
        sweep = self.airframe.le_sweep

        span = b / 2

        # -------------------------
        # Wing outline (2D reference)
        # -------------------------
        y = np.linspace(-b/2, b/2, 300)

        chord = c_r - (np.abs(y) / span) * (c_r - c_t)
        x_le = np.abs(y) * np.tan(sweep)   # FIXED
        x_te = x_le + chord
        print(X)
        # -------------------------
        # Create figure
        # -------------------------
        fig = go.Figure()

        # -------------------------
        # 1. Wing surface
        # -------------------------
        fig.add_trace(go.Surface(
            x=X,
            y=Y,
            z=Z,
            colorscale='Viridis',
            showscale=True,
            name="Wing surface",
            opacity=1.0
        ))

        # -------------------------
        # 2. Leading edge
        # -------------------------
        fig.add_trace(go.Scatter3d(
            x=x_le,
            y=y,
            z=np.zeros_like(y),
            mode='lines',
            line=dict(color='black', width=5),
            name="Leading edge"
        ))

        # -------------------------
        # 3. Trailing edge
        # -------------------------
        fig.add_trace(go.Scatter3d(
            x=x_te,
            y=y,
            z=np.zeros_like(y),
            mode='lines',
            line=dict(color='black', width=5),
            name="Trailing edge"
        ))

        # -------------------------
        # 4. Chordwise guide lines (IMPORTANT FOR SWEEP VISIBILITY)
        # -------------------------
        for yi in np.linspace(-b/2, b/2, 10):

            chord_i = c_r - (np.abs(yi) / span) * (c_r - c_t)
            x_le_i = np.abs(yi) * np.tan(sweep)

            fig.add_trace(go.Scatter3d(
                x=[x_le_i, x_le_i + chord_i],
                y=[yi, yi],
                z=[0, 0],
                mode='lines',
                line=dict(color='gray', width=2),
                showlegend=False
            ))

        # -------------------------
        # Layout
        # -------------------------
        fig.update_layout(
            title=f"HAPS Wing Panel Mesh ({self.panel_size} m × {self.panel_size} m) with Sweep + Thickness",
            scene=dict(
                xaxis_title="x (m)",
                yaxis_title="y (m)",
                zaxis_title="thickness (m)",
                aspectmode="data"
            ),
            margin=dict(l=0, r=0, t=40, b=0)
        )

        fig.show()

if __name__ == "__main__":
    airframe = Airframe.airframe(qc_sweep=np.radians(15), init_polar=False, display=False)
    airfoil_geometry = SparGeometryParam.AirfoilGeometry(Airframe = airframe)
    battery_cross_section = Available_BatteryCrossSection(AirGEO=airfoil_geometry, width_clamp=0.05)
    print(f"Battery Cross-Section Area: {battery_cross_section.battery_cross_section_area:.6f} m^2")
    print(f"Battery Centroid (x): {battery_cross_section.x_centroid_batt:.6f} m")
    print(f"Battery Centroid (z): {battery_cross_section.z_centroid_batt:.6f} m")
    cg_calculator = CG_comp(x_cg_goal=2.4, airframe=airframe)
