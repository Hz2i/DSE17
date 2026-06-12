import sys
import os
import aerosandbox.numpy as np
import plotly.graph_objects as go
import aerosandbox as asb
import StructuralAnalysis as sa

import matplotlib.pyplot as plt
import matplotlib.patches as pth
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

import Airframe
import Components_Materials
import GeneralSubsystems
import PowerSystem_sizing
import SparGeometryParam
import Sections
import CG_comp

def Plot_Cross_Section(Clamp=False, Control_Surface=False, Skid=False, Battery=False, cut='A'):
    airframe = Airframe.airframe(S=62.94, A=20, qc_sweep=15.0/180*np.pi, taper=1.0, dihedral=0.0 , airfoil=asb.Airfoil("mh91"), display=False, init_polar=False)
    airfoil_geometry = SparGeometryParam.AirfoilGeometry(Airframe = airframe, t_skin_airfoil=0.00015, plot=False)
    spar_optimizer = SparGeometryParam.SparGeometryOptimization(
            I_xx_spar_req=1e-7,  # Placeholder values for required inertia, these should be calculated based on load cases
            I_yy_spar_req=5e-8,
            I_xx_sleeve_req=1e-7,
            I_yy_sleeve_req=5e-8,
            n_sections=8,
            min_eccentricity_factor=1.5,
            airframe=airframe,
            airfoil_geometry=airfoil_geometry,
            Plot=False,
            Optimize=False
        )
    
    r_top = 0.010466944193631492
    t_spar = 0.002499999975015052
    t_sleeve = 0.005728992976302554
    Clamp_width = 0.0768877032096635
    eccentricity_factor = 5
    batt_section = CG_comp.Available_BatteryCrossSection(AirGEO=airfoil_geometry, width_clamp=Clamp_width, plot=False)

    geometry = spar_optimizer.calc_geometry_H_clamp(t_spar, t_sleeve, Clamp_width, eccentricity_factor)
    xcentroid = airfoil_geometry.x_centroid
    ycentroid = airfoil_geometry.y_centroid
    print(f"x,y centroid: {xcentroid:.4f}, {ycentroid:.4f}")
    plt.figure(figsize=(10, 5))
    #PLOTTING
    plt.plot(airfoil_geometry.xu, airfoil_geometry.yu, color='black', label="Upper Surface")
    plt.plot(airfoil_geometry.xl, airfoil_geometry.yl, color='black', label="Lower Surface")
    plt.plot(airfoil_geometry.xu_off, airfoil_geometry.yu_off, color='black', label="Upper Surface Offset", linestyle='--')
    plt.plot(airfoil_geometry.xl_off, airfoil_geometry.yl_off, color='black', label="Lower Surface Offset", linestyle='--')
    # Plot Centroid
    plt.scatter(xcentroid, ycentroid, color='black', s=10, label="Spar Centroid Location, Conservative", zorder=5)

    # Available Region
    y_upper = airfoil_geometry.y_upper
    y_lower = airfoil_geometry.y_lower
    x_min = airfoil_geometry.x_min
    x_max = airfoil_geometry.x_max
    print(f"x_min, x_max: {x_min:.4f}, {x_max:.4f}")
    print(f"y_lower, y_upper: {y_lower:.4f}, {y_upper:.4f}")
    # plt.fill_betweenx([y_lower, y_upper], x_min, x_max, color='gray', alpha=0.3, label='Available Region', zorder=0)
    if Clamp:
        # Clamp (centered at centroid)
        clamp = pth.Rectangle((xcentroid - geometry["Clamp_width"]/2, ycentroid - geometry["Clamp_height"]), 
                            geometry["Clamp_width"], geometry["Clamp_height"]*2, 
                            edgecolor='blue', facecolor='blue', alpha=0.1, label='Clamp')
        plt.gca().add_patch(clamp)
        cutout_half_circ_top = pth.Circle((xcentroid, ycentroid + geometry["Clamp_height"]), r_top, color='white', zorder=3)
        cutout_half_circ_bottom = pth.Circle((xcentroid, ycentroid - geometry["Clamp_height"]), r_top, color='white', zorder=3)
        plt.gca().add_patch(cutout_half_circ_top)
        plt.gca().add_patch(cutout_half_circ_bottom)

        # Sleeve Outer (centered at centroid)
        sleeve_out = pth.Ellipse((xcentroid, ycentroid), 
                            geometry["b_sleeve_out"]*2, geometry["a_sleeve_out"]*2, 
                            edgecolor='orange', facecolor='orange', alpha=0.5, label='Sleeve Outer')
        plt.gca().add_patch(sleeve_out)

        # Sleeve Inner (centered at centroid)
        sleeve_in = pth.Ellipse((xcentroid, ycentroid), 
                            geometry["b_sleeve_in"]*2, geometry["a_sleeve_in"]*2, 
                            edgecolor='orange', facecolor='white', linestyle=':', alpha=0.5, label='Sleeve Inner')
        plt.gca().add_patch(sleeve_in)

        # Rubber (centered at centroid)
        rubber_out = pth.Ellipse((xcentroid, ycentroid), 
                            geometry["b_rubber_out"]*2, geometry["a_rubber_out"]*2, 
                            edgecolor='purple', facecolor='purple', alpha=0.5, label='Rubber Outer Layer')
        rubber_in = pth.Ellipse((xcentroid, ycentroid), 
                            geometry["b_rubber_in"]*2, geometry["a_rubber_in"]*2, 
                            edgecolor='purple', facecolor='white', linestyle=':', alpha=0.5, label='Rubber Inner Layer')
        plt.gca().add_patch(rubber_out)
        plt.gca().add_patch(rubber_in)

    # Spar Outer (centered at centroid)
    spar_out = pth.Ellipse((xcentroid, ycentroid), 
                    geometry["b_spar_out"]*2, geometry["a_spar_out"]*2, 
                    edgecolor='green', facecolor='green', alpha=0.5, label='Spar Outer')
    plt.gca().add_patch(spar_out)

    # Spar Inner (centered at centroid)
    spar_in = pth.Ellipse((xcentroid, ycentroid), 
                        geometry["b_spar_in"]*2, geometry["a_spar_in"]*2, 
                        edgecolor='green', facecolor='white', linestyle=':', alpha=1, label='Spar Inner')
    plt.gca().add_patch(spar_in)













    # plt.axvline(x=airfoil_geometry.x_max_thickness, color='green', linestyle='--', label='Max Thickness Location')
    # plt.axvline(x=airfoil_geometry.x_min, color='green', linestyle=':', label='Available Width Boundaries')
    # plt.axvline(x=airfoil_geometry.x_max, color='green', linestyle=':', label='Available Width Boundaries')
    # plt.axhline(y=airfoil_geometry.y_upper, color='orange', linestyle='--', label='Straight Line Region')
    # plt.axhline(y=airfoil_geometry.y_lower, color='yellow', linestyle='--', label='Straight Line Region')
    # plt.scatter(airfoil_geometry.x_centroid, airfoil_geometry.y_centroid, color='black', s=10, label="Centroid Location, Conservative")
    # plt.fill_betweenx([airfoil_geometry.y_lower, airfoil_geometry.y_upper], airfoil_geometry.x_min, airfoil_geometry.x_max, color='gray', alpha=0.3, label='Available Region')



    # plt.scatter(batt_section.x_u_cut, batt_section.y_u_cut, color='black', label="Upper Surface (Truncated)", s=5, alpha=0.5)
    # plt.scatter(batt_section.x_l_cut, batt_section.y_l_cut, color='black', label="Lower Surface (Truncated)", s=5, alpha=0.5)
    x_left_wall_clamp = batt_section.x_wall_battery
    x_right_wall_clamp = airfoil_geometry.x_max_thickness + (airfoil_geometry.x_max_thickness - batt_section.x_wall_battery)
    x_u_wall = airfoil_geometry.xu[::-1]
    y_u_wall = airfoil_geometry.yu[::-1]
    x_l_wall = airfoil_geometry.xl
    y_l_wall = airfoil_geometry.yl
    y_left_wall_clamp_upper = np.interp(x_left_wall_clamp, x_u_wall, y_u_wall)
    y_left_wall_clamp_lower = np.interp(x_left_wall_clamp, x_l_wall, y_l_wall)
    y_right_wall_clamp_upper = np.interp(x_right_wall_clamp, x_u_wall, y_u_wall)
    y_right_wall_clamp_lower = np.interp(x_right_wall_clamp, x_l_wall, y_l_wall)
    plt.plot([x_left_wall_clamp, x_left_wall_clamp], [y_left_wall_clamp_lower, y_left_wall_clamp_upper], color='black', label='Clamp Wall Location')
    plt.plot([x_right_wall_clamp, x_right_wall_clamp], [y_right_wall_clamp_lower, y_right_wall_clamp_upper], color='black', label='Clamp Wall Location')
    if Battery:
        plt.scatter(batt_section.x_centroid_batt, batt_section.z_centroid_batt, color='black', label="Battery Centroid", s=50)

        #Update Upper and lower surface for filling
        x_batt_section = np.linspace(0, max(batt_section.x_u_cut), 100)
        y_upper_fill = np.interp(x_batt_section, batt_section.x_u_cut, batt_section.y_u_cut)
        y_lower_fill = np.interp(x_batt_section, batt_section.x_l_cut, batt_section.y_l_cut)
        plt.fill_between(x_batt_section, y_lower_fill, y_upper_fill, hatch='///', edgecolor='red', facecolor='none', alpha=0.5, label='Available Battery Region')
        #Rectangle for battery area
        area_batt_used = 0.0093
        l_batt_used = np.sqrt(area_batt_used)
        battery_rect = pth.Rectangle((batt_section.x_centroid_batt - l_batt_used/2, batt_section.z_centroid_batt - l_batt_used/2), 
                                    l_batt_used, l_batt_used, 
                                    edgecolor='black', facecolor='gray', linestyle=':', alpha=1, label='Battery Area Used')
        
        #Insulation for battery area both material and air, at least 1.5 cm insulation and air in between START WITH AIR LAYER THEN INSULATION LAYER
    
        n_layers = 5

        air_total = 0.015
        ins_total = 0.015

        air_step = air_total / n_layers
        ins_step = ins_total / n_layers

        # Total thickness of all layers
        total_offset = air_total + ins_total

        for i in range(n_layers):

            # Outer insulation layer
            outer_offset = total_offset - i*(air_step + ins_step)

            insulation_rect = pth.Rectangle(
                (
                    batt_section.x_centroid_batt - l_batt_used/2 - outer_offset,
                    batt_section.z_centroid_batt - l_batt_used/2 - outer_offset,
                ),
                l_batt_used + 2*outer_offset,
                l_batt_used + 2*outer_offset,
                facecolor='lightgray',
                edgecolor='black',
                hatch='xx',
                label='Insulation' if i == 0 else None,
            )
            plt.gca().add_patch(insulation_rect)

            # Next air layer inside it
            air_offset = outer_offset - ins_step

            air_rect = pth.Rectangle(
                (
                    batt_section.x_centroid_batt - l_batt_used/2 - air_offset,
                    batt_section.z_centroid_batt - l_batt_used/2 - air_offset,
                ),
                l_batt_used + 2*air_offset,
                l_batt_used + 2*air_offset,
                facecolor='white',
                edgecolor='blue',
                hatch='\\\\',
                label='Air' if i == 0 else None,
            )
            plt.gca().add_patch(air_rect)
        plt.gca().add_patch(battery_rect)
        # plt.axvline(airfoil_geometry.x_max_thickness, color='orange', label="Max Thickness", linestyle='-.')


    if Skid:
        #Skid
        x_skid_center = airframe.c_r * 0.4
        x_skid_1 = x_skid_center - airframe.c_r * 0.1
        x_skid_2 = x_skid_center + airframe.c_r * 0.1


        #Skid Strut
        l_strut = 0.3
        d_strut = 0.02

        left_strut = pth.Rectangle((x_skid_1, 0), d_strut, -l_strut, edgecolor='brown', facecolor='brown', alpha=0.5, label='Skid Strut')
        plt.gca().add_patch(left_strut)
        right_strut = pth.Rectangle((x_skid_2, 0), -d_strut, -l_strut, edgecolor='brown', facecolor='brown', alpha=0.5, label='Skid Strut')
        plt.gca().add_patch(right_strut)
        #Ski_surface
        l_ski = 0.9
        t_ski = 0.004

        # Find insertion index
        x_left = x_skid_center - l_ski/2
        x_right = x_skid_center + l_ski/2
        idx_left = np.searchsorted(airfoil_geometry.xl, x_left)
        idx_right = np.searchsorted(airfoil_geometry.xl, x_right)

        # Interpolate y-values at wall
        y_wall_left = np.interp(x_left, airfoil_geometry.xl, airfoil_geometry.yl)
        y_wall_right = np.interp(x_right, airfoil_geometry.xl, airfoil_geometry.yl)

        x_l_cut = airfoil_geometry.xl[idx_left:idx_right]
        y_l_cut = airfoil_geometry.yl[idx_left:idx_right]

        x_skid_top = x_l_cut
        y_skid_top = y_l_cut - l_strut + 0.1
        y_skid_bottom = y_l_cut - l_strut - t_ski + 0.1
        plt.fill_between(x_skid_top, y_skid_bottom, y_skid_top, color='brown', alpha=0.5, label='Skid Surface')


    x_control_surface_start = airframe.c_r * 0.7
    x_u_wall_control = airfoil_geometry.xu[::-1]
    y_u_wall_control = airfoil_geometry.yu[::-1]
    x_l_wall_control = airfoil_geometry.xl
    y_l_wall_control = airfoil_geometry.yl

    y_control_surface_upper = np.interp(x_control_surface_start, x_u_wall_control, y_u_wall_control)
    y_control_surface_lower = np.interp(x_control_surface_start, x_l_wall_control, y_l_wall_control)
    plt.plot([x_control_surface_start, x_control_surface_start], [y_control_surface_lower, y_control_surface_upper], color='black', label='Control Surface Start Location')

    if Control_Surface:
        #Control Surface
        
        x_control_rod = x_control_surface_start + airframe.c_r * 0.05
        x_u = airfoil_geometry.xu[::-1]
        y_u = airfoil_geometry.yu[::-1]
        #Interpolate y-values at control surface rod location
        idx_control_rod_l = np.searchsorted(airfoil_geometry.xl, x_control_rod)
        idx_control_rod_u = np.searchsorted(x_u, x_control_rod)

        y_control_rod_lower = airfoil_geometry.yl[idx_control_rod_l]
        y_control_rod_upper = y_u[idx_control_rod_u]

        y_control_rod = (y_control_rod_lower+y_control_rod_upper)/2
        radius_control_surface = abs(y_control_rod - y_control_rod_upper)

        circle_control_surface = pth.Circle((x_control_rod, y_control_rod), radius_control_surface, edgecolor='cyan', facecolor='cyan', alpha=0.5, label='Control Surface Area')
        plt.gca().add_patch(circle_control_surface)
        plt.scatter(x_control_rod, y_control_rod, color='black', label='Control Rod Location', s=50)

        #Contour left side of the circle
        theta = np.linspace(-np.pi/2, np.pi/2, 100)
        x_control_surface_contour = x_control_rod - radius_control_surface * np.cos(theta)
        y_control_surface_contour = y_control_rod - radius_control_surface * np.sin(theta)
        plt.plot(x_control_surface_contour, y_control_surface_contour, color='red', label='Control Surface Contour')
        #Update Upper and lower surface for filling
        x_control_surface_array = np.linspace(x_control_rod, max(airfoil_geometry.xu), 100)
        y_upper_fill = np.interp(x_control_surface_array, x_u, y_u)
        y_lower_fill = np.interp(x_control_surface_array, airfoil_geometry.xl, airfoil_geometry.yl)
        plt.fill_between(x_control_surface_array, y_lower_fill, y_upper_fill, hatch='///', edgecolor='cyan', facecolor='none', alpha=0.5, label='Available Battery Region')
        plt.plot(x_control_surface_array, y_upper_fill, color='red', label='Control Surface Contour')
        plt.plot(x_control_surface_array, y_lower_fill, color='red', label='Control Surface Contour')

    #Central Connection
    beam = pth.Rectangle((x_right_wall_clamp, 0),  (x_control_surface_start - x_right_wall_clamp),d_strut, edgecolor='brown', facecolor='brown', alpha=0.5, label='Skid Strut')
    plt.gca().add_patch(beam)
    plt.title(f"Airfoil Cross-Section at cut {cut}")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    # plt.legend()
    plt.axis('equal')
    plt.grid()
    plt.show()

if __name__ == "__main__":
    Plot_Cross_Section(Clamp=True, Control_Surface=True, Skid=True, Battery=True, cut = 'A')