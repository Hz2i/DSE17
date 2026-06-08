import numpy as np
import aerosandbox as asb
import matplotlib.patches as pth
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

# import sys
# import os
# # Add the folder containing characteristics_airframe.py to the Python path
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

# # Now you can import the Airframe module
# airfoil = asb.Airfoil("e344")
# from StructuralAnalysis import airfoil_properties

class AirfoilGeometry:
    def __init__(self, airfoil_geometry, plot=False):
        self.Airfoil = airfoil_geometry.get("airfoil")
        self.chord_length = airfoil_geometry.get("chord_length")
        self.Available_width = airfoil_geometry.get("Available width")
        self.t_skin_airfoil = airfoil_geometry.get("t_skin_airfoil")
        self.Available_Height_Spar = airfoil_geometry.get("Available_Height_Spar")
        self.xu, self.yu, self.xl, self.yl, self.xu_off, self.yu_off, self.xl_off, self.yl_off, self.x_max_thickness, self.x_min, self.x_max, self.y_upper, self.y_lower, self.h_region, self.x_centroid, self.y_centroid = self.find_available_region()

        if plot:
            self.airfoil_cross_section_plot = self.plot_airfoil_cross_section()
    def find_available_region(self):
        # Implement method to find available region of the airfoil
        upper_surface = self.chord_length * self.Airfoil.upper_coordinates()
        lower_surface = self.chord_length * self.Airfoil.lower_coordinates()
        #Get curves
        def parametric_spline(surface_xy, n=300):
            pts = np.asarray(surface_xy, float)

            # remove consecutive duplicate points
            d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
            keep = np.hstack([[True], d > 1e-12])
            pts = pts[keep]

            # cumulative arc length parameter s
            ds = np.linalg.norm(np.diff(pts, axis=0), axis=1)
            s = np.hstack([[0.0], np.cumsum(ds)])

            # if there are repeated points, s can have duplicates; make it strictly increasing
            s_unique, idx = np.unique(s, return_index=True)
            pts = pts[idx]
            s = s_unique

            sx = interpolate.CubicSpline(s, pts[:, 0], bc_type="natural")
            sy = interpolate.CubicSpline(s, pts[:, 1], bc_type="natural")

            s_new = np.linspace(s.min(), s.max(), n)
            x_new = sx(s_new)
            y_new = sy(s_new)
            return x_new, y_new
        
        xu, yu = parametric_spline(upper_surface, n=40000)
        xl, yl = parametric_spline(lower_surface, n=40000)

        #Offset splines inwards by t_skin_airfoil to account for skin thickness
        def offset_spline(x, y, t):
            dx = np.gradient(x)
            dy = np.gradient(y)
            length = np.sqrt(dx**2 + dy**2)
            nx = -dy / length
            ny = dx / length
            x_offset = x + t * nx
            y_offset = y + t * ny
            return x_offset, y_offset
        
        xu_off, yu_off = offset_spline(xu, yu, self.t_skin_airfoil)
        xl_off, yl_off = offset_spline(xl, yl, self.t_skin_airfoil)

        #Update offsetted splines so it stops at intersection point of upper and lower surface

        def trim_remove_points_right_of_intersection_TE_to_LE(xu, yu, xl, yl, n_grid=4000):
            """
            Inputs are TE -> LE (x typically decreasing).
            Finds x_stop where upper and lower intersect (yu(x)=yl(x)) approximately,
            then removes all points with x > x_stop (i.e. trims away TE-side / right side).
            Returns trimmed arrays (still TE -> LE) and x_stop.
            """
            xu = np.asarray(xu, float); yu = np.asarray(yu, float)
            xl = np.asarray(xl, float); yl = np.asarray(yl, float)

            # For interpolation we need increasing x
            xu_i, yu_i = (xu[::-1], yu[::-1]) if xu[0] > xu[-1] else (xu, yu)
            xl_i, yl_i = (xl[::-1], yl[::-1]) if xl[0] > xl[-1] else (xl, yl)

            # Shared x-grid
            x_min = max(xu_i.min(), xl_i.min())
            x_max = min(xu_i.max(), xl_i.max())
            if x_max <= x_min:
                raise ValueError("Upper and lower x ranges do not overlap.")

            xg = np.linspace(x_min, x_max, n_grid)
            yu_g = np.interp(xg, xu_i, yu_i)
            yl_g = np.interp(xg, xl_i, yl_i)
            diff = yu_g - yl_g

            # Intersection ~ where diff is closest to 0
            i = int(np.argmin(np.abs(diff)))
            x_stop = float(xg[i])

            # Remove all points right of x_stop => remove higher x => keep x <= x_stop
            mask_u = xu <= x_stop
            mask_l = xl <= x_stop

            # Safety: ensure not empty
            if np.count_nonzero(mask_u) < 2 or np.count_nonzero(mask_l) < 2:
                # If this happens, x_stop is too close to LE/TE or data doesn't actually intersect.
                # Return untrimmed but still give x_stop.
                return xu, yu, xl, yl, x_stop

            return xu[mask_u], yu[mask_u], xl[mask_l], yl[mask_l], x_stop

        xu_off, yu_off, xl_off, yl_off, x_stop = trim_remove_points_right_of_intersection_TE_to_LE(xu_off, yu_off, xl_off, yl_off)

        #Find Coordinates Based on Available Width
        x_max_thickness = self.Airfoil.local_thickness(np.linspace(0, 1, 1001)).argmax()/1000*self.chord_length
        print(f"Max thickness location: {x_max_thickness:.4f} m from LE")
        print(f"Thickness at that location: {self.Airfoil.max_thickness()/self.chord_length:.4f} fraction of chord")
        #Define interval around max thickness location
        interval = self.Available_width * self.chord_length / 2
        x_min = x_max_thickness - interval
        x_max = x_max_thickness + interval

        # Make x increasing for interpolation
        xu_i, yu_i = (xu_off[::-1], yu_off[::-1]) if xu_off[0] > xu_off[-1] else (xu_off, yu_off)
        xl_i, yl_i = (xl_off[::-1], yl_off[::-1]) if xl_off[0] > xl_off[-1] else (xl_off, yl_off)

        y_min_upper = np.interp(x_min, xu_i, yu_i)
        y_max_upper = np.interp(x_max, xu_i, yu_i)
        y_min_lower = np.interp(x_min, xl_i, yl_i)
        y_max_lower = np.interp(x_max, xl_i, yl_i)

        #Use Value closest to y = 0 boundary
        y_upper = y_min_upper if abs(y_min_upper) < abs(y_max_upper) else y_max_upper
        y_lower = y_min_lower if abs(y_min_lower) < abs(y_max_lower) else y_max_lower

        h_region = np.abs(y_upper - y_lower)
        #Find Centroid of Available Region
        x_centroid = (x_min + x_max) / 2
        y_centroid = (y_upper + y_lower) / 2
        return xu, yu, xl, yl, xu_off, yu_off, xl_off, yl_off,x_max_thickness, x_min, x_max, y_upper, y_lower, h_region, x_centroid, y_centroid
    
    def plot_airfoil_cross_section(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.xu, self.yu, color='blue', label="Upper Surface")
        plt.plot(self.xl, self.yl, color='red', label="Lower Surface")
        plt.plot(self.xu_off, self.yu_off, color='blue', label="Upper Surface Offset", linestyle='--')
        plt.plot(self.xl_off, self.yl_off, color='red', label="Lower Surface Offset", linestyle='--')
        plt.axvline(x=self.x_max_thickness, color='green', linestyle='--', label='Max Thickness Location')
        plt.axvline(x=self.x_min, color='green', linestyle=':', label='Available Width Boundaries')
        plt.axvline(x=self.x_max, color='green', linestyle=':', label='Available Width Boundaries')
        plt.axhline(y=self.y_upper, color='orange', linestyle='--', label='Straight Line Region')
        plt.axhline(y=self.y_lower, color='yellow', linestyle='--', label='Straight Line Region')
        plt.scatter(self.x_centroid, self.y_centroid, color='black', s=100, label="Centroid Location, Conservative")
        plt.fill_betweenx([self.y_lower, self.y_upper], self.x_min, self.x_max, color='gray', alpha=0.3, label='Available Region')
        plt.title("Airfoil Cross-Section with Available Region")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.legend()
        plt.axis('equal')
        plt.grid()
        plt.show()

class SparGeometryOptimization:
    def __init__(self, min_eccentricity_factor, required_geometry, wing_properties, airfoil_geometry, optimize_variables, determined_geometry, material_properties):
        self.I_xx_spar_req = required_geometry.get("I_xx_spar")
        self.I_yy_spar_req = required_geometry.get("I_yy_spar")
        self.I_xx_sleeve_req = required_geometry.get("I_xx_sleeve")
        self.I_yy_sleeve_req = required_geometry.get("I_yy_sleeve")

        #Min Eccentricity Factor for Ribs
        self.min_eccentricity_factor = min_eccentricity_factor        

        #Initial Guess Geometry
        self.t_spar = optimize_variables.get("t_spar")
        self.t_sleeve = optimize_variables.get("t_sleeve")
        self.t_connection = optimize_variables.get("clamp_width")
        self.eccentricity_factor = optimize_variables.get("eccentricity_factor")

        #Determined Geometry
        self.M6_Bolt = determined_geometry.get("M6_Bolt")
        self.t_bolts = self.M6_Bolt * determined_geometry.get("Bolt_safety_factor")  # Bolt diameter times safety factor for structural integrity
        self.l_clamp = self.t_bolts * determined_geometry.get( "Clamp_width_safety_factor")  # Clamp width must be large enough to accommodate bolts with safety factor
        self.t_rubber = determined_geometry.get("t_rubber")
        self.clamp_minimum_thickness = self.t_bolts * determined_geometry.get("Clamp_thickness_fraction_of_t_bolt")  # Minimum clamp thickness to ensure structural integrity and accommodate bolts
        self.l_sleeve = determined_geometry.get("l_sleeve")

        #Constraint Geometry
        self.airfoil_geometry = airfoil_geometry
        self.t_skin = self.airfoil_geometry.t_skin_airfoil
        self.Available_width = self.airfoil_geometry.Available_width * self.airfoil_geometry.chord_length
        self.Available_height = self.airfoil_geometry.h_region

        #Wing Properties
        self.span = wing_properties.get("span")
        self.sweep = wing_properties.get("sweep")
        self.slanted_span = self.span / np.cos(np.radians(self.sweep))

        #Min Limits of Thicknesses based on manufacturability and material limits
        self.CFRP_LIMIT_t = material_properties.get("CFRP_LIMIT_t")
        self.GLARE_LIMIT_t = material_properties.get("GLARE_LIMIT_t")
        self.GLARE_rho = material_properties.get("density_GLARE")
        self.CFRP_rho = material_properties.get("density_carbon_fiber_composite")
        self.Alu_rho = material_properties.get("density_aluminum")
        self.Rubber_rho = material_properties.get("density_rubber")
        self.PET_rho = material_properties.get("density_PET")

    def calc_I_square(self, b, h):
        return (b * h**3) / 12
    
    def calc_I_circular(self, d):
        return (np.pi * d**4) / 64
    
    def calc_I_xx_ellipse(self, a, b):
        return (np.pi * b * a**3) / 4
    
    def calc_I_yy_ellipse(self, a, b):
        return (np.pi * a * b**3) / 4
    
    def calc_parallel_axis_I(self, I_centroid, A, pos):
        return I_centroid + A * pos**2
    
    def calc_A_square(self, b, h):
        return b * h
    
    def calc_A_circular(self, d):
        return (np.pi * d**2) / 4
    
    def calc_A_ellipse(self, a, b):
        return np.pi * a * b
    
    def calc_geometry_H_clamp(self, t_spar, t_sleeve, Clamp_width, eccentricity_factor):
        #Dims
        h_clamp = self.Available_height/2
        r_top = (Clamp_width-2*self.t_bolts)/2

        #Semi Major
        a_sleeve_out = h_clamp - r_top - self.clamp_minimum_thickness  #Outer dimension of sleeve must fit within clamp height minus bolt radius and clearance
        a_sleeve_in = a_sleeve_out - t_sleeve
        a_rubber_out = a_sleeve_in
        a_rubber_in = a_rubber_out - self.t_rubber
        a_spar_out = a_rubber_in
        a_spar_in = a_spar_out - t_spar
        #Semi Minor
        b_sleeve_out = a_sleeve_out / eccentricity_factor
        b_sleeve_in = a_sleeve_in / eccentricity_factor
        b_rubber_out = a_rubber_out / eccentricity_factor
        b_rubber_in = a_rubber_in / eccentricity_factor
        b_spar_out = a_spar_out / eccentricity_factor
        b_spar_in = a_spar_in / eccentricity_factor
        
        # Areas
        A_sleeve = self.calc_A_ellipse(a_sleeve_out, b_sleeve_out) - self.calc_A_ellipse(a_sleeve_in, b_sleeve_in)
        A_rubber = self.calc_A_ellipse(a_rubber_out, b_rubber_out) - self.calc_A_ellipse(a_rubber_in, b_rubber_in)
        A_spar = self.calc_A_ellipse(a_spar_out, b_spar_out) - self.calc_A_ellipse(a_spar_in, b_spar_in)
        #There are 2 clamps, so multiply by 2
        A_clamp_1 = self.calc_A_square(Clamp_width, h_clamp) - self.calc_A_circular(r_top*2) - self.calc_A_ellipse(a_sleeve_out, b_sleeve_out)/2
        A_clamp_2x = A_clamp_1 * 2

        #Inertia
        I_xx_sleeve = self.calc_I_xx_ellipse(a_sleeve_out, b_sleeve_out) - self.calc_I_xx_ellipse(a_sleeve_in, b_sleeve_in)
        I_yy_sleeve = self.calc_I_yy_ellipse(b_sleeve_out, a_sleeve_out) - self.calc_I_yy_ellipse(b_sleeve_in, a_sleeve_in)
        I_xx_spar = self.calc_I_xx_ellipse(a_spar_out, b_spar_out) - self.calc_I_xx_ellipse(a_spar_in, b_spar_in)
        I_yy_spar = self.calc_I_yy_ellipse(b_spar_out, a_spar_out) - self.calc_I_yy_ellipse(b_spar_in, a_spar_in)

        #Cross Sectional Weight
        weight_clamp = A_clamp_2x * self.Alu_rho
        weight_sleeve = A_sleeve * self.GLARE_rho
        weight_spar = A_spar * self.CFRP_rho
        weight_rubber = A_rubber * self.Rubber_rho
        total_sl_and_sp_weight = weight_sleeve + weight_spar + weight_rubber
        total_weight = weight_clamp + weight_sleeve + weight_spar + weight_rubber

        return {
            #Echoed Inputs
            "t_spar": t_spar,
            "t_sleeve": t_sleeve,
            "Clamp_width": Clamp_width,
            "eccentricity_factor": eccentricity_factor,
            # Dims
            "r_top": r_top,
            "Clamp_height": h_clamp,
            "a_sleeve_out": a_sleeve_out,
            "a_sleeve_in": a_sleeve_in,
            "a_rubber_out": a_rubber_out,
            "a_rubber_in": a_rubber_in,
            "a_spar_out": a_spar_out,
            "a_spar_in": a_spar_in,
            "b_sleeve_out": b_sleeve_out,
            "b_sleeve_in": b_sleeve_in,
            "b_rubber_out": b_rubber_out,
            "b_rubber_in": b_rubber_in,
            "b_spar_out": b_spar_out,
            "b_spar_in": b_spar_in,
            # Areas
            "A_sleeve": A_sleeve,
            "A_rubber": A_rubber,
            "A_spar": A_spar,
            "A_clamp_1": A_clamp_1,
            "A_clamp_2x": A_clamp_2x,
            # Inertia
            "I_xx_sleeve": I_xx_sleeve,
            "I_yy_sleeve": I_yy_sleeve,
            "I_xx_spar": I_xx_spar,
            "I_yy_spar": I_yy_spar,
            # Weight
            "weight_clamp": weight_clamp,
            "weight_sleeve": weight_sleeve,
            "weight_spar": weight_spar,
            "weight_rubber": weight_rubber,
            "total_sl_and_sp_weight": total_sl_and_sp_weight,
            "total_weight": total_weight,
        }

    def optimize_geometry_H_clamp_with_asb(self):
        import aerosandbox as asb
        import aerosandbox.numpy as np
        
        opti = asb.Opti()

        t_spar = opti.variable(init_guess=float(self.t_spar), lower_bound=self.CFRP_LIMIT_t, upper_bound=0.02)
        t_sleeve = opti.variable(init_guess=float(self.t_sleeve), lower_bound=self.GLARE_LIMIT_t, upper_bound=0.02)
        clamp_width = opti.variable(init_guess=float(self.t_connection), lower_bound=0.02, upper_bound=min(self.Available_width, 0.3))
        eccentricity_factor = opti.variable(init_guess=float(self.eccentricity_factor), lower_bound=self.min_eccentricity_factor, upper_bound=10.0)

        g = self.calc_geometry_H_clamp(t_spar, t_sleeve, clamp_width, eccentricity_factor)

        # ===== FEASIBILITY CONSTRAINTS =====
        # Basic bounds
        opti.subject_to(g["Clamp_width"] <= self.Available_width)
        opti.subject_to(g["Clamp_height"]*2 <= self.Available_height)
        opti.subject_to(g["t_spar"] >= self.CFRP_LIMIT_t)
        opti.subject_to(g["t_sleeve"] >= self.GLARE_LIMIT_t)
        
        # Ensure clamp can physically fit the sleeve (loosen this constraint)
        opti.subject_to(g["Clamp_width"] >= (2*self.t_bolts + 2*g["b_sleeve_out"]))  # Clamp width must accommodate the sleeve and bolts

        #Required inertia constraints
        opti.subject_to(g["I_xx_spar"] >= float(self.I_xx_spar_req))
        opti.subject_to(g["I_yy_spar"] >= float(self.I_yy_spar_req))
        opti.subject_to(g["I_xx_sleeve"] >= float(self.I_xx_sleeve_req))
        opti.subject_to(g["I_yy_sleeve"] >= float(self.I_yy_sleeve_req))

        # ===== OBJECTIVE =====
        opti.minimize(g["total_weight"])

        # ===== SOLVE =====
        sol = opti.solve()  # Add verbose to see solver output
        
        result = {}
        for k, v in g.items():
            try:
                result[k] = float(sol(v))
            except (RuntimeError, TypeError):
                result[k] = float(v)
        
        return result
    
    def Plot_optimized_geometry(self, r_top, t_spar, t_sleeve, Clamp_width, eccentricity_factor):
        geometry = self.calc_geometry_H_clamp(t_spar, t_sleeve, Clamp_width, eccentricity_factor)
        xcentroid = self.airfoil_geometry.x_centroid
        ycentroid = self.airfoil_geometry.y_centroid
        plt.figure(figsize=(10, 5))
        
        # Plot Centroid
        plt.scatter(xcentroid, ycentroid, color='black', s=10, label="Centroid Location, Conservative", zorder=5)
        
        # Available Region
        y_upper = self.airfoil_geometry.y_upper
        y_lower = self.airfoil_geometry.y_lower
        x_min = self.airfoil_geometry.x_min
        x_max = self.airfoil_geometry.x_max
        plt.fill_betweenx([y_lower, y_upper], x_min, x_max, color='gray', alpha=0.3, label='Available Region', zorder=0)
        
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
        
        plt.title(f"Optimized H-Clamp Geometry with Eccentricity: {eccentricity_factor:.2f}")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.legend()
        plt.axis('equal')
        plt.grid()

    def calc_Mass_structure_span(self, n_sections, optimized_geometry):
        total_length_inc_clamp = n_sections*(self.l_clamp*3)
        total_length_sleeved = self.l_sleeve*n_sections - total_length_inc_clamp
        total_span_exc_sleeve_clamp = self.slanted_span - total_length_inc_clamp - total_length_sleeved
        mass_per_clamped_cross_section = optimized_geometry["total_weight"]
        mass_per_unclamped_sleeved_cross_section = optimized_geometry["total_sl_and_sp_weight"]  # Only include sleeve weight for unclamped sections, as clamp occupies the space of the spar in those sections
        mass_per_unclamped_unsleeved_cross_section = optimized_geometry["weight_spar"]  # Exclude clamp weight for unclamped sections
        total_mass_spar = (total_length_inc_clamp * mass_per_clamped_cross_section +
                          total_length_sleeved * mass_per_unclamped_sleeved_cross_section +
                          total_span_exc_sleeve_clamp * mass_per_unclamped_unsleeved_cross_section)
        return total_mass_spar, total_length_inc_clamp, total_span_exc_sleeve_clamp
    
    def closed_contour_area(self,xu, yu, xl, yl):
        xu = np.asarray(xu)
        yu = np.asarray(yu)
        xl = np.asarray(xl)
        yl = np.asarray(yl)

        # Determine orientation of each surface
        upper_le_to_te = xu[-1] > xu[0]
        lower_le_to_te = xl[-1] > xl[0]

        # We want:
        # upper: LE -> TE
        # lower: TE -> LE

        if not upper_le_to_te:
            xu = xu[::-1]
            yu = yu[::-1]

        if lower_le_to_te:
            xl = xl[::-1]
            yl = yl[::-1]

        # Remove duplicate TE point if present
        if np.isclose(xu[-1], xl[0]) and np.isclose(yu[-1], yl[0]):
            xl = xl[1:]
            yl = yl[1:]

        # Remove duplicate LE point if present
        if np.isclose(xu[0], xl[-1]) and np.isclose(yu[0], yl[-1]):
            xl = xl[:-1]
            yl = yl[:-1]

        # Closed contour
        x = np.concatenate([xu, xl])
        y = np.concatenate([yu, yl])
        
        # Shoelace formula
        area = 0.5 * np.abs(
            np.dot(x, np.roll(y, -1))
            - np.dot(y, np.roll(x, -1))
        )

        return area, x, y
    def calculate_airfoil_skin_weight(self, optimized_geometry, total_length_inc_clamp, total_span_exc_clamp):
        xu, yu = self.airfoil_geometry.xu, self.airfoil_geometry.yu
        xl, yl = self.airfoil_geometry.xl, self.airfoil_geometry.yl

        area_airfoil_total,x_total, y_total = self.closed_contour_area(xu, yu, xl, yl)

        xu_off, yu_off = self.airfoil_geometry.xu_off, self.airfoil_geometry.yu_off
        xl_off, yl_off = self.airfoil_geometry.xl_off, self.airfoil_geometry.yl_off

        area_airfoil_offset, x_offset, y_offset = self.closed_contour_area(
            xu_off, yu_off,
            xl_off, yl_off
        )

        area_airfoil_skin = area_airfoil_total - area_airfoil_offset
        
        # plt.figure(figsize=(10, 5))
        # plt.plot(x_total, y_total, color='blue', label="Outer Airfoil Contour")
        # plt.plot(x_offset, y_offset, color='red', label="Inner Airfoil Contour (offset by skin thickness)")
        # plt.fill_between(x_total, y_total, color='blue', alpha=0.3, label=f'Outer Airfoil Area: {area_airfoil_total:.4f} m^2')
        # plt.fill_between(x_offset, y_offset, color='red', alpha=0.3, label=f'Inner Airfoil Area: {area_airfoil_offset:.4f} m^2')
        # plt.title(f"Airfoil Contours and Skin Area, total skin area = {area_airfoil_skin:.4f} m^2")
        # plt.xlabel("x (m)")
        # plt.ylabel("y (m)")
        # plt.legend()
        # plt.axis('equal')
        # plt.grid()
        # plt.show()

        clamp_occupied_volume = (optimized_geometry["Clamp_width"] - optimized_geometry["r_top"]*2) * self.t_skin * total_length_inc_clamp  # Approximate volume occupied by clamp within the skin area
        volume_airfoil_skin = area_airfoil_skin * self.slanted_span - clamp_occupied_volume
        weight_airfoil_skin = volume_airfoil_skin * self.PET_rho

        return weight_airfoil_skin

def ShowEcc_vs_Mass(required_geometry, wing_properties, airfoil_geometry, optimize_variables, determined_geometry, material_properties):
    eccentricity_factors = np.linspace(1.1, 6.0, 100)
    Total_weights = []
    Skin_weights = []
    Spar_Structure_weights = []
    for ef in eccentricity_factors:
        geometry = SparGeometryOptimization(ef, required_geometry, wing_properties, airfoil_geometry, optimize_variables, determined_geometry, material_properties)
        optimized_geometry = geometry.optimize_geometry_H_clamp_with_asb()
        weight_spar_structure, l_inc_clamp, l_exc_clamp = geometry.calc_Mass_structure_span(n_sections=4, optimized_geometry=optimized_geometry)
        weight_airfoil_skin = geometry.calculate_airfoil_skin_weight(optimized_geometry, total_length_inc_clamp=l_inc_clamp, total_span_exc_clamp=l_exc_clamp)
        Total_weights.append(weight_spar_structure + weight_airfoil_skin)
        Spar_Structure_weights.append(weight_spar_structure)
        Skin_weights.append(weight_airfoil_skin)
        geometry.Plot_optimized_geometry(optimized_geometry["r_top"], optimized_geometry["t_spar"], optimized_geometry["t_sleeve"], optimized_geometry["Clamp_width"], optimized_geometry["eccentricity_factor"])
    # plt.show()
    plt.figure(figsize=(10, 5))

    # Three subplots: Total Weight, Spar Structure Weight, Skin Weight

    plt.subplot(1, 3, 1)
    plt.plot(eccentricity_factors, Total_weights, marker='o', label="Total Weight")
    plt.title("Total Mass vs Eccentricity Factor")
    plt.xlabel("Eccentricity Factor (a/b)")
    plt.ylabel("Total Mass of Clamp + Sleeve + Spar (kg/m)")
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(eccentricity_factors, Spar_Structure_weights,
            marker='o', color='orange', label="Spar Structure Weight")
    plt.title("Spar Structure Mass vs Eccentricity Factor")
    plt.xlabel("Eccentricity Factor (a/b)")
    plt.ylabel("Mass of Spar + Sleeve (kg/m)")
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(eccentricity_factors, Skin_weights,
            marker='o', color='green', label="Skin Weight")
    plt.title("Airfoil Skin Mass vs Eccentricity Factor")
    plt.xlabel("Eccentricity Factor (a/b)")
    plt.ylabel("Mass of Airfoil Skin (kg/m)")
    plt.grid()
    plt.legend()

    plt.show()
        



if __name__ == "__main__":
    minimum_eccentricity_factor = 1.1  # Minimum eccentricity factor to ensure the sleeve is not too close to circular, which would reduce the difference between I_xx and I_yy and limit design flexibility
    required_geometry = {
        "I_xx_spar": 1e-7, #m^4
        "I_yy_spar": 0.5e-7, #m^4
        "I_xx_sleeve": 1e-7, #m^4
        "I_yy_sleeve": 0.5e-7, #m^4
    }
    airfoil_geometry = {
        "airfoil": asb.Airfoil("e344"),  # Use NACA 2412 as an example
        "chord_length": 1.2,  #m
        "Available width": 0.2,  #Fraction of chord length
        "t_skin_airfoil": 0.0002,  #m
    }
    optimize_variables = {
        "t_spar": 0.005,  #m
        "t_sleeve": 0.005,  #m
        "clamp_width": 0.1,  #m
        "eccentricity_factor": 2, #dimensionless, ratio of semi major to semi minor axis of ellipses
    }

    determined_geometry = {
        "M6_Bolt": 0.006, #m
        "Bolt_safety_factor": 3.0, #dimensionless
        "Clamp_width_safety_factor": 5.0, #dimensionless
        "Clamp_thickness_fraction_of_t_bolt": 0.33, #dimensionless
        "t_rubber": 0.001, #m
        "l_sleeve": 1.0, #Sleeve length as a factor of clamp width to ensure proper load transfer and structural integrity
    }
    material_properties = {
        "density_titanium": 4500.0, #kg/m^3
        "density_aluminum": 2700.0, #kg/m^3
        "density_carbon_fiber_composite": 1800.0, #kg/m^3 (carbon fiber composite)
        "density_GLARE": 2500.0, #kg/m^3 (for the rubber layer in the clamp)
        "density_rubber": 1300.0, #kg/m^3 (for the rubber layer in the clamp)
        "density_PET": 1200.0, #kg/m^3 (for comparison)
        "CFRP_LIMIT_t": 0.001, #m, minimum manufacturable thickness of CFRP
        "GLARE_LIMIT_t": 0.001, #m, minimum manufacturable thickness of GLARE
    }
    wing_properties = {
        "span": 30.1, #m
        "sweep": 15.0, #degrees
    }
    constraint_geometry = AirfoilGeometry(airfoil_geometry, plot=True)
    #Assume Centroid in middle of main spar
    geometry_optimizer = SparGeometryOptimization(minimum_eccentricity_factor, required_geometry, wing_properties, constraint_geometry, optimize_variables, determined_geometry, material_properties)
    optimized_geometry = geometry_optimizer.optimize_geometry_H_clamp_with_asb()
    print("Optimized Geometry:", optimized_geometry)
    geometry_optimizer.Plot_optimized_geometry(
        optimized_geometry["r_top"],
        optimized_geometry["t_spar"],
        optimized_geometry["t_sleeve"],
        optimized_geometry["Clamp_width"],
        optimized_geometry["eccentricity_factor"]
    )
    plt.show()
    weight_spar_structure, l_inc_clamp, l_exc_clamp = geometry_optimizer.calc_Mass_structure_span(n_sections=4, optimized_geometry=optimized_geometry)
    weight_airfoil_skin = geometry_optimizer.calculate_airfoil_skin_weight(optimized_geometry, total_length_inc_clamp=l_inc_clamp, total_span_exc_clamp=l_exc_clamp)
    print("Estimated Mass of Structural Components Along Span:", weight_spar_structure, "kg")
    print("Estimated Mass of Airfoil Skin:", weight_airfoil_skin, "kg")
    print(f"Total Estimated Mass of Spar Structure + Airfoil Skin: {weight_spar_structure + weight_airfoil_skin:.2f} kg")
    # ShowEcc_vs_Mass(required_geometry, wing_properties, constraint_geometry, optimize_variables, determined_geometry, material_properties)