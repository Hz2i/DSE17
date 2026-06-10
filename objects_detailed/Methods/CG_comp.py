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
import Sections


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
    def __init__(self, x_cg_goal=2.4, batt_section=1, airframe=Airframe.airframe(),
                #  Sol=PowerSystem_sizing.power_generation(),
                #  Batt=PowerSystem_sizing.power_storage(),
                 Fcs = GeneralSubsystems.FlightConditionsSystem(),
                 Computer = GeneralSubsystems.ComputerSystem(),
                 Comms = GeneralSubsystems.CommunicationSystem(), 
                 Payload = GeneralSubsystems.PayloadSystem(),
                 Batt = Available_BatteryCrossSection(AirGEO = SparGeometryParam.AirfoilGeometry(), width_clamp=0.125),
                 AirGEO = SparGeometryParam.AirfoilGeometry(),
                 Wing_sections = Sections.Sections(),
                 battery_char = Components_Materials.battery(),
                 t_skin_airfoil=0.0002):
        
        #Selection
        self.selected_section = batt_section # [0,1,2,3,4,5,6]
        #section
        self.l_main_section = Wing_sections.l_main_section
        self.l_main_section_spar = Wing_sections.l_main_section_spar
        self.wing_section_length = Wing_sections.l_middle_wing_section
        self.l_main_section_inside = self.l_main_section-self.l_main_section_spar
        self.wing_tip_section_length = Wing_sections.l_wingtip_section
        self.AirGEO = AirGEO
        self.Batt = Batt
        self.airframe = airframe
        self.battery_char = battery_char
        self.airfoil = airframe.foil
        self.x_cg_goal = x_cg_goal
        self.t_skin_airfoil = t_skin_airfoil
        #Define Centers of Section
        self.x_edge_first_wing_section = self.l_main_section_inside*np.sin(self.airframe.le_sweep)
        self.y_edge_first_wing_section = self.l_main_section_inside*np.cos(self.airframe.le_sweep)
        self.x_center_1 = (self.l_main_section_inside + self.wing_section_length*(0.5))*np.sin(self.airframe.le_sweep)
        self.y_center_1 = (self.l_main_section_inside + self.wing_section_length*(0.5))*np.cos(self.airframe.le_sweep)
        self.x_edge_second_wing_section = (self.l_main_section_inside + self.wing_section_length*(1))*np.sin(self.airframe.le_sweep)
        self.y_edge_second_wing_section = (self.l_main_section_inside + self.wing_section_length*(1))*np.cos(self.airframe.le_sweep)
        self.x_center_2 = (self.l_main_section_inside + self.wing_section_length*(1.5))*np.sin(self.airframe.le_sweep)
        self.y_center_2 = (self.l_main_section_inside + self.wing_section_length*(1.5))*np.cos(self.airframe.le_sweep)
        self.x_edge_third_wing_section = (self.l_main_section_inside + self.wing_section_length*(2))*np.sin(self.airframe.le_sweep)
        self.y_edge_third_wing_section = (self.l_main_section_inside + self.wing_section_length*(2))*np.cos(self.airframe.le_sweep)
        self.x_center_3 = (self.l_main_section_inside + self.wing_section_length*(2.5))*np.sin(self.airframe.le_sweep)
        self.y_center_3 = (self.l_main_section_inside + self.wing_section_length*(2.5))*np.cos(self.airframe.le_sweep)
        self.x_edge_wing_tip_section = (self.l_main_section_inside + self.wing_section_length*(3))*np.sin(self.airframe.le_sweep)
        self.y_edge_wing_tip_section = (self.l_main_section_inside + self.wing_section_length*(3))*np.cos(self.airframe.le_sweep)
        self.x_center_4 = (self.l_main_section_inside + self.wing_section_length*(3)+self.wing_tip_section_length*(0.5))*np.sin(self.airframe.le_sweep)
        self.y_center_4 = (self.l_main_section_inside + self.wing_section_length*(3)+self.wing_tip_section_length*(0.5))*np.cos(self.airframe.le_sweep)
        #Spar Structure 
        '''PUT IN BY HAND'''
        #Lengths
        self.l_clamp = 0.09
        self.l_sleeve = 1
        self.x_cg_spar_structure_wrt_LE = self.AirGEO.x_max_thickness
        #CrossSectionWeight
        self.cs_w_clamp = 2.78 #Hollow Clamp Assumption
        self.cs_w_rubber = 0.068
        self.cs_w_sleeve = 1.623
        self.cs_w_spar = 0.18
        # self.plot_mesh_plotly()
        #Centroid of Airfoil
        self.x_centroid, self.z_centroid = self.Calculate_Airfoil_Centroid()
        '''Constraint Geometry'''
        self.x_cg_skin()
        self.x_cg_wingtips()
        self.x_cg_solar_panel()
        self.x_clamp, self.y_clamp = self.x_cg_spar_structure()
        '''Changable Positions'''
        #Masses All below 1kg is omitted
        self.mass_payload = Payload.mass_payload + Payload.mass_mounting
        self.mass_fcs_pitot = Fcs.pitot_mass * Fcs.redundancy
        self.mass_computer = 0 #No components above 1kg, so not included
        self.mass_comms_ssr = Comms.ssr_mass*Comms.redundancy
        self.mass_comms_elt = Comms.elt_mass*Comms.redundancy
        '''PUT IN BY HAND'''
        self.mass_battery_total = 129.679 #[kg]
        self.volume_battery_tot = self.mass_battery_total / self.battery_char.massRho
        '''PUT IN BY HAND'''
        self.prop_eng = 16.129/2 #[kg] Combining 2!!
        self.prop_sticklength = 0.5 #[m]
        '''PUT IN BY HAND'''
        self.total_skid_mass = 9.15 #[kg]
        self.N_skids = 8
        self.mass_skid = self.total_skid_mass / self.N_skids
        self.skid_chord_pos = 0.4*self.airframe.c_r
        '''PUT IN BY HAND'''
        self.mass_actuator = 1.95 #[kg]
        self.actuator_chord_pos = 0.7*self.airframe.c_r
        self.outer_elevon_offset = 0.9*self.airframe.b/2
        self.outer_elevon_length = 0.185*self.airframe.b/2
        self.inner_elevon_length = 0.051*self.airframe.b/2
        self.elevon_chord = 0.3*self.airframe.c_r

        #Positions
        

        self.y_payload = 0
        self.y_fcs_pitot = 0
        self.y_computer = 0
        self.y_comms_ssr = 0
        self.y_comms_elt = 0
        self.y_prop2_3 = self.airframe.b/6
        self.y_prop1_4 = self.airframe.b/6*2
        self.y_skids_1 = self.y_edge_first_wing_section + 0.5*np.cos(self.airframe.le_sweep)
        self.y_skids_2 = self.airframe.b/6
        self.y_skids_3 = self.airframe.b/6*2
        self.y_skids_4 = self.y_wingtips - 0.5*np.cos(self.airframe.le_sweep)#Offset of 0.5m
        self.y_outer_elevons = self.outer_elevon_offset
        self.y_inner_elevons = self.y_outer_elevons - self.outer_elevon_length

        self.x_payload = self.x_centroid
        self.x_fcs_pitot = self.x_centroid
        self.x_computer = self.x_centroid
        self.x_comms_ssr = self.x_centroid
        self.x_comms_elt = self.x_centroid
        self.x_prop_2_3, self.x_prop_1_4 = self.calc_x_prop_positions()
        self.x_skids_1 = self.x_edge_first_wing_section + 0.5*np.sin(self.airframe.le_sweep)+self.skid_chord_pos
        self.x_skids_2 = self.x_prop_2_3 + self.prop_sticklength + self.skid_chord_pos
        self.x_skids_3 = self.x_prop_1_4 + self.prop_sticklength + self.skid_chord_pos
        self.x_skids_4 = self.x_wingtips - 0.5*np.sin(self.airframe.le_sweep)#Offset of 0.5m
        self.x_outer_elevons = (self.outer_elevon_offset)*np.tan(self.airframe.le_sweep) + self.actuator_chord_pos
        self.x_inner_elevons = (self.y_inner_elevons)*np.tan(self.airframe.le_sweep) + self.actuator_chord_pos
        #Make CGlists
        self.mass_list = self.gen_mass_list()
        self.x_list = self.gen_x_list()
        self.y_list = self.gen_y_list()
        self.name_list = self.gen_name_list()
        #Guesses
        self.l_batt_dist_guess = self.wing_section_length
        
        #Skip y_cg as assumption is symmetric
        self.x_cg_no_batt, self.y_cg_no_batt = self.calculate_cg()
        self.distance_to_goal = abs(self.x_cg_no_batt - self.x_cg_goal)
        print(f'Goal CG: {self.x_cg_goal:.3f} m')
        print(f"Initial CG without battery: {self.x_cg_no_batt:.3f} m, Distance to Goal: {self.distance_to_goal:.3f} m")
        #Print all starting points of the batteries
        print(f'Point 0: {self.x_edge_first_wing_section}')
        print(f'Point 1/2: {self.x_edge_second_wing_section}')
        print(f'Point 3/4: {self.x_edge_third_wing_section}')
        print(f'Point 5/6: {self.x_edge_wing_tip_section}')
        
        #Optimize Battery Distribution
        optimized_batt_dist, resulting_cg, batt_area = self.Optimizer_CG()
        #Plot CG
        self.plot_CG(optimized_batt_dist)
        print(f"Optimized Battery Distribution Length: {optimized_batt_dist:.3f} m, Resulting CG: {resulting_cg:.3f} m, Battery Area Used: {batt_area:.6f} m^2")

    def calc_x_prop_positions(self):
        #Equally spaced along span.
        d_eng = self.airframe.b/6
        x_eng_2_3 = d_eng*np.tan(self.airframe.le_sweep) - self.prop_sticklength
        x_eng_1_4 = x_eng_2_3 + d_eng*np.tan(self.airframe.le_sweep)
        return x_eng_2_3, x_eng_1_4

    def Batt_Distribution(self, l_batt_dist):
        batt_total_vol_half= self.volume_battery_tot/2
        batt_total_mass = self.mass_battery_total
        #Start Either at Beginning of Main Section or at End of Main Section, depending on selected section
        if self.selected_section % 2 == 0:
            #Find Datum Point
            x_datum_bat = (self.l_main_section_inside + self.wing_section_length*self.selected_section/2)*np.sin(self.airframe.le_sweep)
            y_datum_bat = (self.l_main_section_inside + self.wing_section_length*self.selected_section/2)*np.cos(self.airframe.le_sweep)
            half_batt_l = l_batt_dist/2
            batt_centroid_x = x_datum_bat + self.Batt.x_centroid_batt + half_batt_l*np.sin(self.airframe.le_sweep)
            batt_centroid_y = y_datum_bat + half_batt_l*np.cos(self.airframe.le_sweep)
            batt_cross_section_filled = batt_total_vol_half/l_batt_dist

        if self.selected_section % 2 == 1:
        #Find Datum Point
            x_datum_bat = (self.l_main_section_inside + self.wing_section_length*(self.selected_section+1)/2)*np.sin(self.airframe.le_sweep)
            y_datum_bat = (self.l_main_section_inside + self.wing_section_length*(self.selected_section+1)/2)*np.cos(self.airframe.le_sweep)
            half_batt_l = l_batt_dist/2
            batt_centroid_x = x_datum_bat + self.Batt.x_centroid_batt - half_batt_l*np.sin(self.airframe.le_sweep)
            batt_centroid_y = y_datum_bat - half_batt_l*np.cos(self.airframe.le_sweep)
            batt_cross_section_filled = batt_total_vol_half/l_batt_dist
        return np.array([
                            batt_centroid_x,
                            batt_centroid_y,
                            batt_cross_section_filled,
                            batt_total_mass,
                        ], dtype=object)

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
        self.half_span_wing_tips = 1.9
        self.sweep_wing_tips = self.airframe.le_sweep  #Assuming more sweep for wingtips
        self.x_wingtips = ((self.airframe.b/2)*np.tan(self.airframe.le_sweep))+self.x_centroid + 0.5*self.half_span_wing_tips*np.tan(self.sweep_wing_tips)  # Shift x-centroid by half the span times tan(sweep) to account for sweep, plus additional shift for wingtip geometry
        self.y_wingtips = self.airframe.b/2  
        self.mass_wingtips = 3.7 #[kg] (Assumptions)

    def x_cg_solar_panel(self):
        self.solar_panel_area = 54.2
        self.mass_solar_panel = self.solar_panel_area * Components_Materials.solar_panel().surfRho
        self.x_solar_panel = ((self.airframe.b/2)*np.tan(self.airframe.le_sweep))/2 + self.x_centroid
        self.y_solar_panel = ((self.airframe.b/2)/2)

    def x_cg_spar_structure(self):
        #Clamps + Rubber
        n_clamps = 3
        x_datum_clamps_and_rubber = np.array([self.x_center_1, self.x_center_2, self.x_center_3, self.x_center_4]) + self.x_cg_spar_structure_wrt_LE
        y_datum_clamps_and_rubber = np.array([self.y_center_1, self.y_center_2, self.y_center_3, self.y_center_4])
        mass_clamps_and_rubber = np.ones(len(x_datum_clamps_and_rubber)) * (n_clamps * self.l_clamp * self.cs_w_clamp+ n_clamps * self.l_clamp * self.cs_w_rubber)
        #Sleeve
        n_sleeves = 1
        x_datum_sleeve = np.array([self.x_center_1, self.x_center_2, self.x_center_3, self.x_center_4])+ self.x_cg_spar_structure_wrt_LE
        y_datum_sleeve = np.array([self.y_center_1, self.y_center_2, self.y_center_3, self.y_center_4])
        mass_sleeve = np.ones(len(x_datum_sleeve))*n_sleeves * self.l_sleeve * self.cs_w_sleeve
        #Spar, Along total length
        x_datum_spar = np.array([self.airframe.b / 2 * np.tan(self.airframe.le_sweep) / 2 + self.x_cg_spar_structure_wrt_LE])
        y_datum_spar = np.array([self.airframe.b / 2 / 2])
        mass_spar = np.array([(self.airframe.b/2)/np.cos(self.airframe.le_sweep) * self.cs_w_spar])
        #Combine CG of Spar Structure
        x_spar_list = np.concatenate([x_datum_clamps_and_rubber, x_datum_sleeve, x_datum_spar])
        y_spar_list = np.concatenate([y_datum_clamps_and_rubber, y_datum_sleeve, y_datum_spar])
        mass_spar_structure_list = np.concatenate([mass_clamps_and_rubber, mass_sleeve, mass_spar])
        #Final X,Y, CG for Spar Structure
        self.x_spar = np.sum(x_spar_list * mass_spar_structure_list) / np.sum(mass_spar_structure_list)
        self.y_spar = np.sum(y_spar_list * mass_spar_structure_list) / np.sum(mass_spar_structure_list)
        self.mass_spar = np.sum(mass_spar_structure_list)
        return x_datum_clamps_and_rubber, y_datum_clamps_and_rubber
    
    def gen_mass_list(self):
        mass_list = [
            self.mass_payload,
            self.mass_fcs_pitot,
            self.mass_computer,
            self.mass_comms_ssr,
            self.mass_comms_elt,
            self.mass_spar,
            self.mass_solar_panel,
            self.prop_eng,
            self.prop_eng,
            self.mass_skid*2,
            self.mass_skid*2,
            self.mass_skid*2,
            self.mass_skid*2,
            self.mass_actuator*2, #Outer
            self.mass_actuator*2, #Inner
            self.mass_wingtips,
            self.mass_skin,
            0, #mass_total_battery
        ]
        return mass_list
    
    def gen_x_list(self):
        x_list = [
            self.x_payload,
            self.x_fcs_pitot,
            self.x_computer,
            self.x_comms_ssr,
            self.x_comms_elt,
            self.x_spar,
            self.x_solar_panel,
            self.x_prop_2_3, #Assuming propellers are distributed around their CGs
            self.x_prop_1_4,
            self.x_skids_1, #Assuming skids are distributed around their CGs
            self.x_skids_2,
            self.x_skids_3,
            self.x_skids_4,
            self.x_outer_elevons, #Assuming elevons are distributed around their CGs
            self.x_inner_elevons,
            self.x_wingtips, #Assuming wingtips are distributed around their CG
            self.x_skin, #Assuming skin mass is distributed around its centroid
            0, #x_total_battery
        ]
        return x_list
    
    def gen_y_list(self):
        y_list = [
            self.y_payload,
            self.y_fcs_pitot,
            self.y_computer,
            self.y_comms_ssr,
            self.y_comms_elt,
            self.y_spar,
            self.y_solar_panel,
            self.y_prop2_3,
            self.y_prop1_4,
            self.y_skids_1,
            self.y_skids_2,
            self.y_skids_3,
            self.y_skids_4,
            self.y_outer_elevons,
            self.y_inner_elevons,
            self.y_wingtips,
            self.y_skin,
            0, #y_total_battery
        ]
        return y_list
    
    def gen_name_list(self):
        name_list = [
            "Payload",
            "FCS Pitot",
            "Computer",
            "Comms SSR",
            "Comms ELT",
            "Spar",
            "Solar Panel",
            "Propeller 2/3",
            "Propeller 1/4",
            "Skids 1",
            "Skids 2",
            "Skids 3",
            "Skids 4",
            "Outer Elevon Actuators",
            "Inner Elevon Actuators",
            "Wingtips",
            "Skin",
            "Battery"
        ]
        return name_list

    '''Calculate CG'''
    def calculate_cg(self):
        mass_list = np.array(self.mass_list)
        x_list = np.array(self.x_list)
        y_list = np.array(self.y_list)
        x_cg = np.sum(mass_list * x_list) / np.sum(mass_list)
        y_cg = np.sum(mass_list * y_list) / np.sum(mass_list)
        return x_cg, y_cg

    def calculate_cg_batt(self,l_batt_dist):
        batt_list = self.Batt_Distribution(l_batt_dist)
        batt_area = batt_list[2]
        
        # Create lists and append battery values at the end
        mass_list = list(self.mass_list)
        x_list = list(self.x_list)
        y_list = list(self.y_list)
        
        # Replace the last element (battery)
        mass_list[-1] = batt_list[3]
        x_list[-1] = batt_list[0]
        y_list[-1] = batt_list[1]
        
        # Convert to numpy arrays
        mass_list = np.array(mass_list)
        x_list = np.array(x_list)
        y_list = np.array(y_list)
        
        name_list = self.gen_name_list()
        
        # CG calculation (fully vectorized)
        x_cg = np.sum(mass_list * x_list) / np.sum(mass_list)
        y_cg = np.sum(mass_list * y_list) / np.sum(mass_list)
        print(name_list)

        return mass_list, x_list, y_list, name_list, x_cg, batt_area, l_batt_dist
    
    def optimize_function(self, l_batt_dist): 
        _, _, _, _, x_cg_with_batt, batt_area, _ = self.calculate_cg_batt(l_batt_dist)
        return l_batt_dist, x_cg_with_batt, batt_area


    '''Optimizer'''

    def Optimizer_CG(self):
        opti = asb.Opti()
       
        # Define optimization variables
        l_batt_dist = opti.variable(init_guess=self.l_batt_dist_guess, lower_bound=0.5, upper_bound=self.wing_section_length)  # Battery distribution length along the wing
        #Optimize Function
        function = self.optimize_function(l_batt_dist)
        #Feasibility
        opti.subject_to(function[1] == self.x_cg_goal)  # Enforce CG to be exactly at the target (can be relaxed to a range if needed)
        # opti.subject_to(function[1] >= self.x_cg_goal - 0.001)  # Allow for a small tolerance around the target CG
        # opti.subject_to(function[1] <= self.x_cg_goal + 0.001)
        opti.subject_to(function[2] <= 0.9 * self.Batt.battery_cross_section_area)  #Add Safety Factor to Allow for insulation 

        opti.minimize(function[2])#Allow for most insulation
        # Define the objective function (e.g., minimize the difference between main section length and wingtip section length)
        sol = opti.solve(verbose=False)
        return sol.value(l_batt_dist), sol.value(function[1]), sol.value(function[2]) # Return optimized section, battery distribution length, resulting CG, and battery area for feedback
    
    def plot_wing_contour(self, left_half=False):
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
        if left_half:
            y_wing = -y_wing

        return go.Scatter(
            x=x_wing,
            y=y_wing,
            mode="lines",
            name="Wing Contour",
            line=dict(color="black", width=2)
        )
    
    def plot_outer_elevon_contour(self, left_half=False):
        sweep = self.airframe.le_sweep

        # Outer elevon span stations
        y_outer_end = self.outer_elevon_offset
        y_outer_start = y_outer_end - self.outer_elevon_length

        # Leading edge of elevon
        x_le = np.array([
            y_outer_start * np.tan(sweep) + self.actuator_chord_pos,
            y_outer_end   * np.tan(sweep) + self.actuator_chord_pos
        ])

        # Trailing edge
        x_te = x_le + self.elevon_chord

        # Polygon
        x_elevon = np.array([
            x_le[0],
            x_le[1],
            x_te[1],
            x_te[0],
            x_le[0]
        ])

        y_elevon = np.array([
            y_outer_start,
            y_outer_end,
            y_outer_end,
            y_outer_start,
            y_outer_start
        ])

        if left_half:
            y_elevon = -y_elevon

        return go.Scatter(
            x=x_elevon,
            y=y_elevon,
            mode="lines",
            name="Outer Elevon",
            line=dict(color="red", width=2)
        )


    def plot_inner_elevon_contour(self, left_half=False):
        sweep = self.airframe.le_sweep

        

        # Outer elevon
        y_outer_end = self.outer_elevon_offset
        y_outer_start = y_outer_end - self.outer_elevon_length

        # Inner elevon starts where outer elevon ends
        y_inner_end = y_outer_start
        y_inner_start = y_inner_end - self.inner_elevon_length

        # Leading edge
        x_le = np.array([
            y_inner_start * np.tan(sweep) + self.actuator_chord_pos,
            y_inner_end   * np.tan(sweep) + self.actuator_chord_pos
        ])

        # Trailing edge
        x_te = x_le + self.elevon_chord

        # Polygon
        x_elevon = np.array([
            x_le[0],
            x_le[1],
            x_te[1],
            x_te[0],
            x_le[0]
        ])

        y_elevon = np.array([
            y_inner_start,
            y_inner_end,
            y_inner_end,
            y_inner_start,
            y_inner_start
        ])

        if left_half:
            y_elevon = -y_elevon

        return go.Scatter(
            x=x_elevon,
            y=y_elevon,
            mode="lines",
            name="Inner Elevon",
            line=dict(color="blue", width=2)
        )

    def group_by_position(self, x_list, y_list, name_list, mass_list, ndigits=6):
        from collections import defaultdict
        groups = defaultdict(list)

        # force clean numpy-safe conversion
        x_list = np.asarray(x_list).flatten()
        y_list = np.asarray(y_list).flatten()
        name_list = list(name_list)
        mass_list = np.asarray(mass_list).flatten()

        for x, y, name, m in zip(x_list, y_list, name_list, mass_list):
            key = (round(float(x), ndigits), round(float(y), ndigits))
            groups[key].append((name, float(m)))

        return groups
    
    def BendingRelief(self, mass_list=None, y_list=None):
        #Mass * Y * g= Bending Moment
        bending_moment_list = (mass_list) * y_list * 9.81 * 0.5
        #Calculate Bending Relief as the sum of bending moments divided by the total mass
        bending_relief = sum(bending_moment_list)
        return bending_relief
    
    def plot_CG(self, l_batt_dist=None):
        import plotly.io as pio
        pio.renderers.default = "browser"
        fig = go.Figure()
        fig.update_layout(
            xaxis=dict(range=[-1, 6]),
            yaxis=dict(range=[-1, 20])
        )
        x_batt = self.Batt_Distribution(l_batt_dist)[0]
        y_batt = self.Batt_Distribution(l_batt_dist)[1]
        mass_list, x_list, y_list, name_list, x_cg, batt_area, l_batt_dist = self.calculate_cg_batt(l_batt_dist)
        print(f'mass list: {mass_list}')
        bending_relief = self.BendingRelief(mass_list, y_list)
        print(f"Bending Relief: {bending_relief:.2f} N*m")
        # Update Name List so same x       
        # CG points
        for i in range(len(mass_list)):
            #alternate text position 
            # positions = [
            #     "top left",
            #     "top right",
            #     "bottom left",
            #     "bottom right",
            #     "middle right"
            # ]
            # text_position = positions[i % len(positions)]
            fig.add_trace(
                go.Scatter(
                    x=[x_list[i], x_list[i]],
                    y=[y_list[i], -y_list[i]],
                    mode="markers+text",
                    # text = f"{name_list[i]}",
                    # textposition=text_position,
                    name=f"{name_list[i]} - (total mass={mass_list[i]:.2f}), positions(x,y)=({x_list[i]:.2f}, {y_list[i]:.2f}) + ({x_list[i]:.2f}, {-y_list[i]:.2f})",
                    marker=dict(
                        size=mass_list[i]*5/2,
                        opacity=0.7
                    )
                )
            )
           
        #Annotations
        groups = self.group_by_position(x_list, y_list, name_list, mass_list)
        print(f"Groups for annotation: {groups}")
        for i, ((x, y), items) in enumerate(groups.items()):

            # merged label
            label = "<br>".join([
                f"{name} ({m/2:.2f} kg)" for name, m in items
            ])

            scale_arrow = ((sum(m for _, m in items))**(0.5)) * 65

            # flip direction every other element
            direction = 1 if i % 2 == 0 else -1
            #Adjust Arrows of Actuators
            if "Outer Elevon Actuators" in label or "Inner Elevon Actuators" in label:
                print(f"scale_arrow before: {scale_arrow}")
                print('TEST')
                print(f'direction: {direction}')
                fig.add_annotation(
                    x=x,
                    y=y,
                    text=label,
                    showarrow=True,
                    arrowhead=1,
                    font=dict(size=14, family="Arial Black"),
                    ax=direction * 300,
                    ay=0
                )

                # mirrored annotation only if not on centerline
                if y != 0:
                    fig.add_annotation(
                        x=x,
                        y=-y,
                        text=label,
                        showarrow=True,
                        arrowhead=1,
                        font=dict(size=14, family="Arial Black"),
                        ax=direction * 300,
                        ay=0
                    )
            else:
                fig.add_annotation(
                    x=x,
                    y=y,
                    text=label,
                    showarrow=True,
                    arrowhead=1,
                    font=dict(size=14, family="Arial Black"),
                    ax=direction * scale_arrow,
                    ay=0
                )

                # mirrored annotation only if not on centerline
                if y != 0:
                    fig.add_annotation(
                        x=x,
                        y=-y,
                        text=label,
                        showarrow=True,
                        arrowhead=1,
                        font=dict(size=14, family="Arial Black"),
                        ax=direction * scale_arrow,
                        ay=0
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
        # CG NO BATTERY
        x_cg_no_batt = self.calculate_cg()[0]
        fig.add_trace(
            go.Scatter(
                x=[x_cg_no_batt],
                y=[0],
                mode="markers",
                name="CG without Battery",
                marker=dict(symbol="x", size=12, color="blue"),
                opacity=0.5
            )
        )
        # Computed CG
        fig.add_trace(
            go.Scatter(
                x=[x_cg],
                y=[0],
                mode="markers",
                name="Computed CG",
                marker=dict(symbol="x", size=24, color="red"),
                opacity=0.5
            )
        )
        #Battery Distribution Annotation
        fig.add_annotation(
            x=x_batt,
            y=y_batt,
            text=f"Middle Point of Battery Distribution<br>Length: {l_batt_dist:.2f} m<br>Area Used: {batt_area:.4f} m²",
            showarrow=True,
            arrowhead=1,
            font=dict(size=15, family="Arial Black"),
            ax=0,
            ay=-500
        )
        #Battery Distribution Annotation
        fig.add_annotation(
            x=x_batt,
            y=-y_batt,
            text=f"Middle Point of Battery Distribution<br>Length: {l_batt_dist:.2f} m<br>Area Used: {batt_area:.4f} m²",
            showarrow=True,
            arrowhead=1,
            ax=0,
            ay=500
        )
        sweep = self.airframe.le_sweep

        n_points = 100

        # spanwise distribution (CENTERED at y_batt)
        y_batt_dist = np.linspace(-l_batt_dist/2, l_batt_dist/2, n_points) + y_batt

        # LE sweep projects x from y
        x_batt_dist =y_batt_dist * np.tan(sweep) + self.Batt.x_centroid_batt
        fig.add_scatter(
            x=x_batt_dist,
            y=y_batt_dist,
            mode="lines",
            name="Battery Distribution Right Wing",
            line=dict(color="red", width=2)
        )
        fig.add_scatter(
            x=x_batt_dist,
            y=-y_batt_dist,
            mode="lines",
            name="Battery Distribution Left Wing",
            line=dict(color="red", width=2)
        )
        # Wing contour
        fig.add_trace(self.plot_wing_contour())
        fig.add_trace(self.plot_wing_contour(left_half=True))
        # Elevon Contours
        fig.add_trace(self.plot_outer_elevon_contour())
        fig.add_trace(self.plot_outer_elevon_contour(left_half=True))
        fig.add_trace(self.plot_inner_elevon_contour())
        fig.add_trace(self.plot_inner_elevon_contour(left_half=True))
        #
        fig.update_layout(
            title=f"Center of Gravity Distribution - Bending Relief: {bending_relief:.2f} N*m - Available Battery Area: {self.Batt.battery_cross_section_area:.4f} m²",
            xaxis_title="X Coordinate",
            yaxis_title="Y Coordinate",
            showlegend=True
        )

        #Wing Section Lines
        section_positions = [
        self.l_main_section_inside,
        self.l_main_section_inside + self.wing_section_length,
        self.l_main_section_inside + self.wing_section_length * 2,
        self.l_main_section_inside + self.wing_section_length * 3
        ]
        section_names = [
            "Middle Wing Section 1",
            "Middle Wing Section 2",
            "Middle Wing Section 3",
            "Wingtip Section"
        ]
        y_sec_list = []
        for sec in section_positions:
            #Generate Horizontal Lines for Sections
            #Find y coordinate at section
            y_span = sec * np.cos(self.airframe.le_sweep)
            y_sec_list.append(y_span)
            name_sec = section_names[section_positions.index(sec)]
            fig.add_hline(y=y_span, line=dict(color="gray", width=1, dash="dash"), annotation_text=f"{name_sec}", annotation_position="top left")
            fig.add_hline(y=-y_span, line=dict(color="gray", width=1, dash="dash"), annotation_text=f"{name_sec}", annotation_position="bottom left")
        
        #Plot Clamp Positions
        self.x_clamp, self.y_clamp = self.x_cg_spar_structure()
        for i in range(len(self.x_clamp)):
            fig.add_trace(
                go.Scatter(
                    x=[self.x_clamp[i]],
                    y=[self.y_clamp[i]],
                    mode="markers",
                    name=f"Clamp {i+1}",
                    marker=dict(symbol="diamond", size=10, color="purple")
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=[self.x_clamp[i]],
                    y=[-self.y_clamp[i]],
                    mode="markers",
                    name=f"Clamp {i+1} (Mirrored)",
                    marker=dict(symbol="diamond", size=10, color="purple")
                )
            )
        #Layout
        fig.update_layout(
            xaxis=dict(range=[-1, 6]),
            yaxis=dict(range=[-1.5, 20])
        )
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
)
        fig.write_html("cg_plot.html", auto_open=True)
        # fig.show()
    
    


if __name__ == "__main__":
    airframe = Airframe.airframe(S=55.8, A=20, qc_sweep=15.0/180*np.pi, taper=1.0, dihedral=0.0 , airfoil=asb.Airfoil("mh91"), display=False, init_polar=False)
    airfoil_geometry = SparGeometryParam.AirfoilGeometry(Airframe = airframe)
    '''INPUT SELF'''
    battery_cross_section = Available_BatteryCrossSection(AirGEO=airfoil_geometry, width_clamp=0.0569)
    cg_calculator = CG_comp(x_cg_goal=2.275, batt_section=2, airframe=airframe, Batt=battery_cross_section, Wing_sections=Sections.Sections(airframe=airframe, Plot=False), t_skin_airfoil=0.0002)
