import numpy as np
import aerosandbox as asb
import matplotlib.patches as pth
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

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

class GeometryOptimization:
    def __init__(self, required_geometry, airfoil_geometry, optimize_variables_square, optimize_variables_circular, material_properties):
        self.required_thickness = required_geometry.get("Minimim Thickness of Spar and Connector")
        self.I_xx_spar_req = required_geometry.get("I_xx_spar")
        self.I_yy_spar_req = required_geometry.get("I_yy_spar")
        self.I_xx_connection_req = required_geometry.get("I_xx_connection")
        self.I_yy_connection_req = required_geometry.get("I_yy_connection")
        self.I_xx_clamp_req = required_geometry.get("I_xx_clamp")
        self.I_yy_clamp_req = required_geometry.get("I_yy_clamp")

        #Initial Guess Geometry Square
        self.t_spar = optimize_variables_square.get("t_spar")
        self.t_connection = optimize_variables_square.get("t_connection")
        self.Available_Height_Spar = optimize_variables_square.get("Available_Height_Spar")

        #Initial Guess Geometry Circular
        self.t_spar = optimize_variables_circular.get("t_spar")
        self.t_bolts = optimize_variables_circular.get("t_bolts")
        self.t_clamp = optimize_variables_circular.get("t_clamp")

        #Constraint Geometry
        self.airfoil_geometry = airfoil_geometry
        self.width = self.airfoil_geometry.Available_width * self.airfoil_geometry.chord_length

        #Min Limits of Thicknesses based on manufacturability and material limits
        self.CFRP_LIMIT_t = material_properties.get("CFRP_LIMIT_t")
        self.Bolt_LIMIT_t = material_properties.get("Bolt_LIMIT_t")
        self.Clamp_LIMIT_t = material_properties.get("Clamp_LIMIT_t")

        '''Square'''
        self.h_connection = self.airfoil_geometry.h_region * ((1 - self.Available_Height_Spar)/2)
        self.h_spar = self.airfoil_geometry.h_region * self.Available_Height_Spar  

        '''Circular'''
        self.t_rubber = optimize_variables_circular.get("t_rubber") 
    
    def calc_I_square(self, b, h):
        return (b * h**3) / 12
    
    def calc_I_circular(self, d):
        return (np.pi * d**4) / 64
    
    def calc_J_circular(self, d):
        return (np.pi * d**4) / 32
    
    def calc_parallel_axis_I(self, I_centroid, A, pos):
        return I_centroid + A * pos**2
    
    def calc_A_square(self, b, h):
        return b * h
    
    def calc_A_circular(self, d):
        return (np.pi * d**2) / 4
    
    '''Square'''
    def calc_geometry_square_outputs(self, t_spar, t_connection, f_spar):
        """
        Pure calculator: works with either floats OR asb.Opti symbolic variables.
        Returns all useful outputs in a dict: dims, areas, inertias, totals.
        """
        b = self.width
        h_region = self.airfoil_geometry.h_region

        h_spar_out = h_region * f_spar
        h_conn_out = h_region * (1 - f_spar) / 2

        # Outer dims
        b_spar_out = b
        b_conn_out = b

        # Inner dims
        b_spar_in = b_spar_out - 2 * t_spar
        h_spar_in = h_spar_out - 2 * t_spar

        b_conn_in = b_conn_out - 2 * t_connection
        h_conn_in = h_conn_out - 2 * t_connection

        # Areas
        A_spar = self.calc_A_square(b_spar_out, h_spar_out) - self.calc_A_square(b_spar_in, h_spar_in)
        A_conn = self.calc_A_square(b_conn_out, h_conn_out) - self.calc_A_square(b_conn_in, h_conn_in)

        # Spar inertias (tube = outer - inner)
        I_xx_spar = self.calc_I_square(b_spar_out, h_spar_out) - self.calc_I_square(b_spar_in, h_spar_in)
        I_yy_spar = self.calc_I_square(h_spar_out, b_spar_out) - self.calc_I_square(h_spar_in, b_spar_in)

        # Connector inertias about their own centroid (one connector)
        I_xx_conn_centroid = self.calc_I_square(b_conn_out, h_conn_out) - self.calc_I_square(b_conn_in, h_conn_in)
        I_yy_conn_centroid = self.calc_I_square(h_conn_out, b_conn_out) - self.calc_I_square(h_conn_in, b_conn_in)

        # Parallel axis for Ixx of top+bottom connectors
        d = h_spar_out / 2 + h_conn_out / 2
        I_xx_conn_one = self.calc_parallel_axis_I(I_xx_conn_centroid, A_conn, d)

        I_xx_conn_total = 2 * I_xx_conn_one
        I_yy_conn_total = 2 * I_yy_conn_centroid

        A_total = A_spar + 2 * A_conn

        return {
            # decision vars (echo)
            "t_spar": t_spar,
            "t_connection": t_connection,
            "f_spar": f_spar,

            # dims
            "b_spar_out": b_spar_out, "h_spar_out": h_spar_out,
            "b_spar_in":  b_spar_in,  "h_spar_in":  h_spar_in,
            "b_conn_out": b_conn_out, "h_conn_out": h_conn_out,
            "b_conn_in":  b_conn_in,  "h_conn_in":  h_conn_in,
            "d_conn": d,

            # areas
            "A_spar": A_spar,
            "A_conn_one": A_conn,
            "A_total": A_total,

            # inertias
            "I_xx_spar": I_xx_spar,
            "I_yy_spar": I_yy_spar,
            "I_xx_conn_total": I_xx_conn_total,
            "I_yy_conn_total": I_yy_conn_total,
        }
        
    def calc_circular_geometry_outputs(self, t_spar, t_bolts, t_clamp):
        #Constraints
        t_rubber = self.t_rubber
        b = self.width
        h = self.airfoil_geometry.h_region
        if h > b:
            raise ValueError("Height of available region must be less than width for circular geometry.")
        #Dims
        r_clamp = h/2 - t_clamp - t_rubber #width cannot be used
        r_spar_out = r_clamp - t_rubber
        r_spar_in = r_spar_out - t_spar
        Clamp_width = 2*(r_clamp + t_rubber + t_bolts)
        Clamp_height = 2*(r_clamp + t_rubber + t_clamp)
        
        #Areas
        A_spar = self.calc_A_circular(2*r_spar_out) - self.calc_A_circular(2*r_spar_in)
        A_clamp = self.calc_A_square(b, h) - 2*t_rubber*(b) - self.calc_A_circular(2*r_clamp)
        A_total_CFRP = A_spar + A_clamp
        A_rubber = 2*t_rubber*(b) + self.calc_A_circular(2*r_clamp) - self.calc_A_circular(2*r_spar_out)
        #Spar Inertias
        I_xx_spar = self.calc_I_circular(2*r_spar_out) - self.calc_I_circular(2*r_spar_in)
        I_yy_spar = self.calc_I_circular(2*r_spar_out) - self.calc_I_circular(2*r_spar_in)

        #Clamp Inertias, for conservativeness assume no rubber (area closer to neutral axis)
        I_xx_clamp = self.calc_I_square(b,h) - self.calc_I_circular(2*r_clamp)
        I_yy_clamp = self.calc_I_square(h,b) - self.calc_I_circular(2*r_clamp)
        #NO PAT!!!!
        return {
            # decision vars (echo)
            "t_spar": t_spar,
            "r_clamp": r_clamp,
            "t_bolts": t_bolts,
            "t_clamp": t_clamp,
            # dims
            "r_clamp": r_clamp,
            "r_spar_out": r_spar_out,
            "r_spar_in": r_spar_in,
            "Clamp_width": Clamp_width,
            "Clamp_height": Clamp_height,
            # areas
            "A_total_CFRP": A_total_CFRP,
            "A_spar": A_spar,
            "A_clamp": A_clamp,
            "A_rubber": A_rubber,
            # inertias
            "I_xx_spar": I_xx_spar,
            "I_yy_spar": I_yy_spar,
            "I_xx_clamp": I_xx_clamp,
            "I_yy_clamp": I_yy_clamp,
            }


    def optimize_geometry_square_with_asb(self):
        import aerosandbox as asb
        import aerosandbox.numpy as np  # ensures any np ops you add later are symbolic-safe

        opti = asb.Opti()

        t_spar = opti.variable(init_guess=float(self.t_spar), lower_bound=self.required_thickness, upper_bound=0.10)
        t_conn = opti.variable(init_guess=float(self.t_connection), lower_bound=self.required_thickness, upper_bound=0.05)
        f_spar = opti.variable(init_guess=float(self.Available_Height_Spar), lower_bound=0.2, upper_bound=0.8)

        g = self.calc_geometry_square_outputs(t_spar, t_conn, f_spar)

        # feasibility
        opti.subject_to(g["b_spar_in"] > 1e-4)
        opti.subject_to(g["h_spar_in"] > 1e-4)
        opti.subject_to(g["b_conn_in"] > 1e-4)
        opti.subject_to(g["h_conn_in"] > 1e-4)
        opti.subject_to(g["h_conn_out"] > 1e-4)

        # requirements
        opti.subject_to(g["I_xx_spar"] >= float(self.I_xx_spar_req))
        opti.subject_to(g["I_yy_spar"] >= float(self.I_yy_spar_req))
        opti.subject_to(g["I_xx_conn_total"] >= float(self.I_xx_connection_req))
        opti.subject_to(g["I_yy_conn_total"] >= float(self.I_yy_connection_req))

        # objective
        opti.minimize(g["A_total"])

        sol = opti.solve()

        result = {}
        for k, v in g.items():
            try:
                result[k] = float(sol(v))
            except (RuntimeError, TypeError):
                # If it fails, it's a constant (self.width, etc.)
                result[k] = float(v)
        return result
    
    def optimize_geometry_circular_with_asb(self):
        import aerosandbox as asb
        import aerosandbox.numpy as np  # ensures any np ops you add later are symbolic-safe

        opti = asb.Opti()

        t_spar = opti.variable(init_guess=float(self.t_spar), lower_bound=self.CFRP_LIMIT_t, upper_bound=0.10)
        t_bolts = opti.variable(init_guess=float(self.t_bolts), lower_bound=self.Bolt_LIMIT_t, upper_bound=0.05)
        t_clamp = opti.variable(init_guess=float(self.t_clamp), lower_bound=self.Clamp_LIMIT_t, upper_bound=0.10)

        g = self.calc_circular_geometry_outputs(t_spar, t_bolts, t_clamp)

        # feasibility
        b = self.width
        h = self.airfoil_geometry.h_region
        opti.subject_to(g["Clamp_width"] <= b)
        opti.subject_to(g["Clamp_height"] <= h)

        # requirements
        opti.subject_to(g["I_xx_spar"] >= float(self.I_xx_spar_req))
        opti.subject_to(g["I_yy_spar"] >= float(self.I_yy_spar_req))
        opti.subject_to(g["I_xx_clamp"] >= float(self.I_xx_clamp_req))
        opti.subject_to(g["I_yy_clamp"] >= float(self.I_yy_clamp_req))

        # objective
        opti.minimize(g["A_total_CFRP"])

        sol = opti.solve()

        result = {}
        for k, v in g.items():
            try:
                result[k] = float(sol(v))
            except (RuntimeError, TypeError):
                # If it fails, it's a constant (self.width, etc.)
                result[k] = float(v)
        return result
    
    def plot_rect_lines(self, xc, yc, b, h, color, linestyle="-", label=None):
        xL, xR = xc - b / 2, xc + b / 2
        yB, yT = yc - h / 2, yc + h / 2

        # bottom + top
        plt.plot([xL, xR], [yB, yB], color=color, linestyle=linestyle, label=label)
        plt.plot([xL, xR], [yT, yT], color=color, linestyle=linestyle)

        # left + right (verticals)
        plt.plot([xL, xL], [yB, yT], color=color, linestyle=linestyle)
        plt.plot([xR, xR], [yB, yT], color=color, linestyle=linestyle)

    def plot_geometry_square(self,b_spar_out, h_spar_out,   b_spar_in,  h_spar_in,
        b_connection_out, h_connection_out,
        b_connection_in,  h_connection_in
        ):

        xc = self.airfoil_geometry.x_centroid
        yc = self.airfoil_geometry.y_centroid

        plt.figure(figsize=(10, 5))

        # Middle spar
        self.plot_rect_lines(xc, yc, b_spar_out, h_spar_out, color="blue", linestyle="-",  label="Spar Outer")
        self.plot_rect_lines(xc, yc, b_spar_in,  h_spar_in,  color="blue", linestyle="--", label="Spar Inner")

        # Top connector (above spar)
        yc_top = yc + h_spar_out / 2 + h_connection_out / 2
        self.plot_rect_lines(xc, yc_top, b_connection_out, h_connection_out, color="red", linestyle="-",  label="Top Connector Outer")
        # inner connector centered in the same connector layer:
        self.plot_rect_lines(xc, yc_top, b_connection_in,  h_connection_in,  color="red", linestyle="--", label="Top Connector Inner")

        # Bottom connector (below spar)
        yc_bot = yc - h_spar_out / 2 - h_connection_out / 2
        self.plot_rect_lines(xc, yc_bot, b_connection_out, h_connection_out, color="red", linestyle="-",  label="Bottom Connector Outer")
        self.plot_rect_lines(xc, yc_bot, b_connection_in,  h_connection_in,  color="red", linestyle="--", label="Bottom Connector Inner")

        plt.title("Optimized Geometry Cross-Section (3 layers: connector / spar / connector)")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.axis("equal")
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_geometry_circular(self, r_clamp, r_spar_out, r_spar_in, t_spar, t_bolts, t_clamp):
        xc = self.airfoil_geometry.x_centroid
        yc = self.airfoil_geometry.y_centroid

        plt.figure(figsize=(10, 5))
        #Available Region
        y_upper = self.airfoil_geometry.y_upper
        y_lower = self.airfoil_geometry.y_lower
        x_min = self.airfoil_geometry.x_min
        x_max = self.airfoil_geometry.x_max
        plt.fill_betweenx([y_lower, y_upper], x_min, x_max, color='gray', alpha=0.3, label='Available Region')

        # Spar (tube)
        circle_outer = plt.Circle((xc, yc), r_spar_out, color="blue", fill=False, linestyle="-", label="Spar Outer")
        circle_inner = plt.Circle((xc, yc), r_spar_in, color="blue", fill=False, linestyle="--", label="Spar Inner")
        plt.gca().add_patch(circle_outer)
        plt.gca().add_patch(circle_inner)
        
        #Rubber btw spar and clamp
        rubber = pth.Wedge(
            (xc, yc),
            r_clamp,
            0, 360,
            width=r_clamp - r_spar_out,
            facecolor="red",
            alpha=0.4,
            label="Rubber"
        )
        plt.gca().add_patch(rubber)
        # Clamp (Offsetted Arc Part))
        circle_clamp = plt.Circle((xc, yc), r_clamp, color="purple", fill=False, linestyle=":", label="Clamp Outer")
        plt.gca().add_patch(circle_clamp)
        half_w = r_clamp + t_bolts
        half_h = r_clamp + t_clamp

        rect_clamp = plt.Rectangle(
            (xc - half_w, yc - half_h),   # bottom-left corner (corrected)
            2 * half_w,
            2 * half_h,
            angle=0,
            edgecolor="purple",
            fill=False,
            linestyle="-",
            label="Clamp Outer"
        )
        plt.gca().add_patch(rect_clamp)
    


        plt.title("Optimized Geometry Cross-Section (circular clamp + circular spar)")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.axis("equal")
        plt.grid(True)
        plt.legend()
        plt.show()
class StructuralWeightCalculator:
    def __init__(self, aircraft_properties, geometry_outputs, airfoil_geometry, material_properties):
        self.aircraft_properties = aircraft_properties
        self.geometry_outputs = geometry_outputs
        self.airfoil_geometry = airfoil_geometry
        self.material_properties = material_properties

    def calculate_weight_supportive_structure(self):
        A_total = self.geometry_outputs["A_total_CFRP"]
        weight_supportive_structure = A_total * (self.aircraft_properties["span"]/np.cos(np.radians(self.aircraft_properties["sweep"]))) * self.material_properties["density_carbon_fiber_composite"]
        return weight_supportive_structure
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
    def calculate_airfoil_skin_weight(self):
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
        
        plt.figure(figsize=(10, 5))
        plt.plot(x_total, y_total, color='blue', label="Outer Airfoil Contour")
        plt.plot(x_offset, y_offset, color='red', label="Inner Airfoil Contour (offset by skin thickness)")
        plt.fill_between(x_total, y_total, color='blue', alpha=0.3, label=f'Outer Airfoil Area: {area_airfoil_total:.4f} m^2')
        plt.fill_between(x_offset, y_offset, color='red', alpha=0.3, label=f'Inner Airfoil Area: {area_airfoil_offset:.4f} m^2')
        plt.title(f"Airfoil Contours and Skin Area, total skin area = {area_airfoil_skin:.4f} m^2")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.legend()
        plt.axis('equal')
        plt.grid()
        plt.show()

        span_effective = (
            self.aircraft_properties["span"]
            / np.cos(np.radians(self.aircraft_properties["sweep"]))
        )

        weight_airfoil_skin = (
            area_airfoil_skin
            * span_effective
            * self.material_properties["LDPET_density"]
        )

        return weight_airfoil_skin

if __name__ == "__main__":
    required_geometry = {
        "I_xx_spar": 6.8e-6, #m^4
        "I_yy_spar": 2.6e-12, #m^4
        "I_xx_connection": 5.2e-7, #m^4
        "I_yy_connection": 2.6e-8, #m^4
        "I_xx_clamp": 5.2e-6, #m^4
        "I_yy_clamp": 2.6e-12, #m^4
        "Minimim Thickness of Spar and Connector": 0.000001, #m
    }
    airfoil_geometry = {
        "airfoil": asb.Airfoil("e344"),  # Use NACA 2412 as an example
        "chord_length": 1.2,  #m
        "Available width": 0.15,  #Fraction of chord length
        "t_skin_airfoil": 0.0002,  #m
    }
    optimize_variables_square = {
        #THESE ARE INITIAL GUESS VALUES, NOT FINAL OPTIMIZED VALUES
        "t_spar": 0.005, #m
        "t_connection": 0.005, #m
        "Available_Height_Spar": 0.6, #m
    }
    optimize_variables_circular = {
        "t_rubber": 0.001, #m
        #THESE ARE INITIAL GUESS VALUES, NOT FINAL OPTIMIZED VALUES
        "t_bolts": 0.025, #m
        "t_spar": 0.001, #m
        "r_clamp": 0.05, #m
        "t_clamp": 0.025, #m
    }
    material_properties = {
        "CFRP_LIMIT_t" : 0.001,
        "Bolt_LIMIT_t": 0.03,
        "Clamp_LIMIT_t": 0.01,
        "density_titanium": 4500.0, #kg/m^3
        "density_carbon_fiber_composite": 1800.0, #kg/m^3 (carbon fiber composite)
        "LDPET_density": 1000.0, #kg/m^3 (for comparison)
    }

    aircraft_properties = {
        "span": 30.1, #m
        "sweep": 15.0, #degrees
    }
    constraint_geometry = AirfoilGeometry(airfoil_geometry, plot=True)
    #Assume Centroid in middle of main spar
    geometry_optimizer = GeometryOptimization(required_geometry, constraint_geometry, optimize_variables_square, optimize_variables_circular, material_properties)
    # optimized_geometry_square = geometry_optimizer.optimize_geometry_square_with_asb()
    # geometry_optimizer.plot_geometry_square(
    #     optimized_geometry_square["b_spar_out"], 
    #     optimized_geometry_square["h_spar_out"],
    #     optimized_geometry_square["b_spar_in"], 
    #     optimized_geometry_square["h_spar_in"],
    #     optimized_geometry_square["b_conn_out"],      
    #     optimized_geometry_square["h_conn_out"],     
    #     optimized_geometry_square["b_conn_in"], 
    #     optimized_geometry_square["h_conn_in"]        
    # )
    optimized_geometry_circular = geometry_optimizer.optimize_geometry_circular_with_asb()
    print(f"Optimized Circular Geometry Outputs: {optimized_geometry_circular}")
    geometry_optimizer.plot_geometry_circular(
        optimized_geometry_circular["r_clamp"],
        optimized_geometry_circular["r_spar_out"],
        optimized_geometry_circular["r_spar_in"],
        optimized_geometry_circular["t_spar"],
        optimized_geometry_circular["t_bolts"],
        optimized_geometry_circular["t_clamp"],
    )
    
    # weight_calculator = StructuralWeightCalculator(aircraft_properties, optimized_geometry_square, constraint_geometry, material_properties)
    # weight_supportive_structure = weight_calculator.calculate_weight_supportive_structure()
    # weight_airfoil_skin = weight_calculator.calculate_airfoil_skin_weight()
    # print(f"Weight of Supportive Structure: {weight_supportive_structure:.2f} kg")
    # print(f"Weight of Airfoil Skin: {weight_airfoil_skin:.2f} kg")
    # print(f"Total Structural Weight: {weight_supportive_structure + weight_airfoil_skin:.2f} kg")