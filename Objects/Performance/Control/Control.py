import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox import Atmosphere


class FlyingWingWithWingletsAeroBuildup:
    def __init__(
        self,
        span=35.479,
        root_chord=1.774,
        tip_chord=1.774,
        sweep_deg=15,
        dihedral_deg=1.5,
        twist_root_deg=0,
        twist_tip_deg=-4.675,
        wing_airfoil="mh91",

        aspect_ratio = 20.0,

        winglet_height=2.0,
        winglet_root_chord=1.774,
        winglet_tip_chord=1.774,
        winglet_sweep_deg=15,
        winglet_cant_deg=0,
        winglet_airfoil="naca0012",

        aileron_tip_gap_frac=0.10,   # gap from tip before aileron starts
        aileron_span_frac=0.20,      # aileron length, going inwards
        elevator_span_frac=0.10,     # elevator immediately inboard of aileron
        control_hinge=0.75,
        rudder_span_frac=0.8,
        rudder_hinge=0.75,

        velocity=32.70,
        altitude=60000 * 0.3048,
        alpha=10.436,
        beta=0,
        x_cg=2.64,
        z_cg=0,

        control_layout="elevator_outboard",  # or "aileron_outboard"
    ):
        self.span = span
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.sweep_deg = sweep_deg
        self.dihedral_deg = dihedral_deg
        self.twist_root_deg = twist_root_deg
        self.twist_tip_deg = twist_tip_deg
        self.wing_airfoil = asb.Airfoil(wing_airfoil)
        self.wing_area = span ** 2 / aspect_ratio

        self.winglet_height = winglet_height
        self.winglet_root_chord = winglet_root_chord
        self.winglet_tip_chord = winglet_tip_chord
        self.winglet_sweep_deg = winglet_sweep_deg
        self.winglet_cant_deg = winglet_cant_deg
        self.winglet_airfoil = asb.Airfoil(winglet_airfoil)

        self.aileron_tip_gap_frac = aileron_tip_gap_frac
        self.aileron_span_frac = aileron_span_frac
        self.elevator_span_frac = elevator_span_frac
        self.control_hinge = control_hinge
        self.rudder_span_frac = rudder_span_frac
        self.rudder_hinge = rudder_hinge

        self.velocity = velocity
        self.altitude = altitude
        self.alpha = alpha
        self.beta = beta

        self.x_cg = x_cg
        self.z_cg = z_cg

        self.control_layout = control_layout

    def chord_at_y(self, y_abs):
        half_span = self.span / 2
        return self.root_chord + (self.tip_chord - self.root_chord) * y_abs / half_span

    def twist_at_y(self, y_abs):
        half_span = self.span / 2
        return self.twist_root_deg + (
            self.twist_tip_deg - self.twist_root_deg
        ) * y_abs / half_span

    def leading_edge_at_y(self, y):
        y_abs = abs(y)
        x = y_abs * np.tand(self.sweep_deg)
        z = y_abs * np.tand(self.dihedral_deg)
        return [x, y, z]

    def make_main_wing(
            self,
            left_aileron=0,
            right_aileron=0,
            left_elevator=0,
            right_elevator=0,
    ):
        half_span = self.span / 2
        y_tip = half_span

        y_control_outer = half_span * (1 - self.aileron_tip_gap_frac)

        if self.control_layout == "elevator_outboard":
            # tip -> elevator -> aileron -> center
            y_elevator_outer = y_control_outer
            y_elevator_inner = y_elevator_outer - half_span * self.elevator_span_frac

            y_aileron_outer = y_elevator_inner
            y_aileron_inner = y_aileron_outer - half_span * self.aileron_span_frac

        elif self.control_layout == "aileron_outboard":
            # tip -> aileron -> elevator -> center
            y_aileron_outer = y_control_outer
            y_aileron_inner = y_aileron_outer - half_span * self.aileron_span_frac

            y_elevator_outer = y_aileron_inner
            y_elevator_inner = y_elevator_outer - half_span * self.elevator_span_frac

        else:
            raise ValueError(
                "control_layout must be either "
                "'elevator_outboard' or 'aileron_outboard'."
            )

        def xsec(y, controls=None):
            if controls is None:
                controls = []

            y_abs = abs(y)

            return asb.WingXSec(
                xyz_le=self.leading_edge_at_y(y),
                chord=self.chord_at_y(y_abs),
                twist=self.twist_at_y(y_abs),
                airfoil=self.wing_airfoil,
                control_surfaces=controls,
            )

        left_aileron_surface = asb.ControlSurface(
            name="left_aileron",
            symmetric=False,
            deflection=left_aileron,
            hinge_point=self.control_hinge,
        )

        right_aileron_surface = asb.ControlSurface(
            name="right_aileron",
            symmetric=False,
            deflection=right_aileron,
            hinge_point=self.control_hinge,
        )

        left_elevator_surface = asb.ControlSurface(
            name="left_elevator",
            symmetric=False,
            deflection=left_elevator,
            hinge_point=self.control_hinge,
        )

        right_elevator_surface = asb.ControlSurface(
            name="right_elevator",
            symmetric=False,
            deflection=right_elevator,
            hinge_point=self.control_hinge,
        )

        if self.control_layout == "elevator_outboard":
            xsecs = [
                xsec(-y_tip),
                xsec(-y_elevator_outer, [left_elevator_surface]),
                xsec(-y_elevator_inner, [left_elevator_surface, left_aileron_surface]),
                xsec(-y_aileron_inner, [left_aileron_surface]),
                xsec(0),
                xsec(y_aileron_inner, [right_aileron_surface]),
                xsec(y_elevator_inner, [right_aileron_surface, right_elevator_surface]),
                xsec(y_elevator_outer, [right_elevator_surface]),
                xsec(y_tip),
            ]

        else:  # "aileron_outboard"
            xsecs = [
                xsec(-y_tip),
                xsec(-y_aileron_outer, [left_aileron_surface]),
                xsec(-y_aileron_inner, [left_aileron_surface, left_elevator_surface]),
                xsec(-y_elevator_inner, [left_elevator_surface]),
                xsec(0),
                xsec(y_elevator_inner, [right_elevator_surface]),
                xsec(y_aileron_inner, [right_elevator_surface, right_aileron_surface]),
                xsec(y_aileron_outer, [right_aileron_surface]),
                xsec(y_tip),
            ]

        return asb.Wing(
            name="Main Flying Wing",
            symmetric=False,
            xsecs=xsecs,
        )

    def make_winglet(self, side="right", rudder_deflection=0):
        sign = 1 if side == "right" else -1
        half_span = self.span / 2

        x_root, y_root, z_root = self.leading_edge_at_y(sign * half_span)

        x_tip = x_root + self.winglet_height * np.tand(self.winglet_sweep_deg)
        y_tip = y_root + sign * self.winglet_height * np.tand(self.winglet_cant_deg)
        z_tip = z_root + self.winglet_height

        rudder = asb.ControlSurface(
            name=f"{side}_rudder",
            symmetric=False,
            deflection=rudder_deflection,
            hinge_point=self.rudder_hinge,
        )

        rudder_frac = self.rudder_span_frac

        x_rudder_end = x_root + rudder_frac * (x_tip - x_root)
        y_rudder_end = y_root + rudder_frac * (y_tip - y_root)
        z_rudder_end = z_root + rudder_frac * (z_tip - z_root)

        return asb.Wing(
            name=f"{side.capitalize()} Winglet",
            symmetric=False,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[x_root, y_root, z_root],
                    chord=self.winglet_root_chord,
                    dihedral=90,
                    airfoil=self.winglet_airfoil,
                    control_surfaces=[rudder],
                ),
                asb.WingXSec(
                    xyz_le=[x_rudder_end, y_rudder_end, z_rudder_end],
                    chord=self.winglet_root_chord + rudder_frac * (
                            self.winglet_tip_chord - self.winglet_root_chord
                    ),
                    dihedral=90,
                    airfoil=self.winglet_airfoil,
                    control_surfaces=[rudder],
                ),
                asb.WingXSec(
                    xyz_le=[x_tip, y_tip, z_tip],
                    chord=self.winglet_tip_chord,
                    dihedral=90,
                    airfoil=self.winglet_airfoil,
                ),
            ],
        )

    def make_airplane(
        self,
        left_aileron=0,
        right_aileron=0,
        left_elevator=0,
        right_elevator=0,
        left_rudder=0,
        right_rudder=0,
    ):
        main_wing = self.make_main_wing(
            left_aileron=left_aileron,
            right_aileron=right_aileron,
            left_elevator=left_elevator,
            right_elevator=right_elevator,
        )

        left_winglet = self.make_winglet(
            side="left",
            rudder_deflection=left_rudder,
        )

        right_winglet = self.make_winglet(
            side="right",
            rudder_deflection=right_rudder,
        )

        s_ref = self.span * (self.root_chord + self.tip_chord) / 2

        taper = self.tip_chord / self.root_chord
        c_ref = (
            2 / 3
            * self.root_chord
            * (1 + taper + taper ** 2)
            / (1 + taper)
        )

        return asb.Airplane(
            name="Flying Wing With Winglets",
            xyz_ref=[self.x_cg, 0, self.z_cg],
            wings=[main_wing, left_winglet, right_winglet],
            s_ref=s_ref,
            b_ref=self.span,
            c_ref=c_ref,
        )

    def make_op_point(self):
        return asb.OperatingPoint(
            velocity=self.velocity,
            alpha=self.alpha,
            beta=self.beta,
            atmosphere=asb.Atmosphere(altitude=self.altitude),
        )

    def run_aerobuildup(self, **control_deflections):
        airplane = self.make_airplane(**control_deflections)

        aero = asb.AeroBuildup(
            airplane=airplane,
            op_point=self.make_op_point(),
        )

        return aero.run_with_stability_derivatives()

    def s_ref_value(self):
        return self.span * (self.root_chord + self.tip_chord) / 2

    def c_ref_value(self):
        taper = self.tip_chord / self.root_chord
        return (
            2 / 3
            * self.root_chord
            * (1 + taper + taper ** 2)
            / (1 + taper)
        )

    def body_axis_coefficients(self, result, alpha_deg):
        alpha = np.radians(alpha_deg)

        CL = result["CL"]
        print(CL)
        CD = result["CD"]
        print(CD)
        CY = result["CY"]
        Cl = result["Cl"]
        Cm = result["Cm"]
        Cn = result["Cn"]

        CX = -CD * np.cos(alpha) + CL * np.sin(alpha)
        CZ = -CD * np.sin(alpha) - CL * np.cos(alpha)

        return {
            "CX": CX,
            "CY": CY,
            "CZ": CZ,
            "Cl": Cl,
            "Cm": Cm,
            "Cn": Cn,
        }

    def run_aerobuildup_custom(
        self,
        velocity=None,
        alpha=None,
        beta=None,
        p=0,
        q=0,
        r=0,
        **control_deflections,
    ):
        if velocity is None:
            velocity = self.velocity
        if alpha is None:
            alpha = self.alpha
        if beta is None:
            beta = self.beta

        airplane = self.make_airplane(**control_deflections)

        op_point = asb.OperatingPoint(
            velocity=velocity,
            alpha=alpha,
            beta=beta,
            p=p,
            q=q,
            r=r,
            atmosphere=asb.Atmosphere(altitude=self.altitude),
        )

        aero = asb.AeroBuildup(
            airplane=airplane,
            op_point=op_point,
        ).run()

        return self.body_axis_coefficients(aero, alpha)

    def compute_all_derivatives(
        self,
        du_frac=0.01,
        dalpha_deg=0.25,
        dbeta_deg=0.25,
        rate_hat=0.01,
        dcontrol_deg=1.0,
    ):
        """
        Returns:
        CXu, CXa, CXq, CXde,
        CZu, CZa, CZq, CZde,
        Cm, Cmu, Cma, Cmq, Cmde,
        CYb, CYp, CYr, CYda, CYdr,
        Clb, Clp, Clr, Clda, Cldr,
        Cnb, Cnp, Cnr, Cnda, Cndr

        Derivative conventions:
            u derivative: per u/V
            alpha, beta: per radian
            p derivative: per p*b/(2V)
            q derivative: per q*c/(2V)
            r derivative: per r*b/(2V)
            controls: per radian
        """

        V = self.velocity
        c_ref = self.c_ref_value()
        b_ref = self.span

        d_alpha_rad = np.radians(dalpha_deg)
        d_beta_rad = np.radians(dbeta_deg)
        d_control_rad = np.radians(dcontrol_deg)

        # Base
        base = self.run_aerobuildup_custom()
        Cm_base = base["Cm"]

        # u derivatives
        plus = self.run_aerobuildup_custom(velocity=V * (1 + du_frac))
        minus = self.run_aerobuildup_custom(velocity=V * (1 - du_frac))

        CXu = (plus["CX"] - minus["CX"]) / (2 * du_frac)
        CZu = (plus["CZ"] - minus["CZ"]) / (2 * du_frac)
        Cmu = (plus["Cm"] - minus["Cm"]) / (2 * du_frac)

        # alpha derivatives
        plus = self.run_aerobuildup_custom(alpha=self.alpha + dalpha_deg)
        minus = self.run_aerobuildup_custom(alpha=self.alpha - dalpha_deg)

        CXa = (plus["CX"] - minus["CX"]) / (2 * d_alpha_rad)
        CZa = (plus["CZ"] - minus["CZ"]) / (2 * d_alpha_rad)
        Cma = (plus["Cm"] - minus["Cm"]) / (2 * d_alpha_rad)

        # beta derivatives
        plus = self.run_aerobuildup_custom(beta=self.beta + dbeta_deg)
        minus = self.run_aerobuildup_custom(beta=self.beta - dbeta_deg)

        CYb = (plus["CY"] - minus["CY"]) / (2 * d_beta_rad)
        Clb = (plus["Cl"] - minus["Cl"]) / (2 * d_beta_rad)
        Cnb = (plus["Cn"] - minus["Cn"]) / (2 * d_beta_rad)

        # pitch-rate derivatives q_hat = q*c/(2V)
        q_rate = rate_hat * V / c_ref # todo removed *2

        plus = self.run_aerobuildup_custom(q=q_rate)
        minus = self.run_aerobuildup_custom(q=-q_rate)

        CXq = (plus["CX"] - minus["CX"]) / (2 * rate_hat)
        CZq = (plus["CZ"] - minus["CZ"]) / (2 * rate_hat)
        Cmq = (plus["Cm"] - minus["Cm"]) / (2 * rate_hat)

        # roll-rate derivatives p_hat = p*b/(2V)
        p_rate = rate_hat * 2 * V / b_ref

        plus = self.run_aerobuildup_custom(p=p_rate)
        minus = self.run_aerobuildup_custom(p=-p_rate)

        CYp = (plus["CY"] - minus["CY"]) / (2 * rate_hat)
        Clp = (plus["Cl"] - minus["Cl"]) / (2 * rate_hat)
        Cnp = (plus["Cn"] - minus["Cn"]) / (2 * rate_hat)

        # yaw-rate derivatives r_hat = r*b/(2V)
        r_rate = rate_hat * 2 * V / b_ref

        plus = self.run_aerobuildup_custom(r=r_rate)
        minus = self.run_aerobuildup_custom(r=-r_rate)

        CYr = (plus["CY"] - minus["CY"]) / (2 * rate_hat)
        Clr = (plus["Cl"] - minus["Cl"]) / (2 * rate_hat)
        Cnr = (plus["Cn"] - minus["Cn"]) / (2 * rate_hat)

        # elevator derivatives
        plus = self.run_aerobuildup_custom(
            left_elevator=dcontrol_deg,
            right_elevator=dcontrol_deg,
        )
        minus = self.run_aerobuildup_custom(
            left_elevator=-dcontrol_deg,
            right_elevator=-dcontrol_deg,
        )

        CXde = (plus["CX"] - minus["CX"]) / (2 * d_control_rad)
        CZde = (plus["CZ"] - minus["CZ"]) / (2 * d_control_rad)
        Cmde = (plus["Cm"] - minus["Cm"]) / (2 * d_control_rad)

        # aileron derivatives
        plus = self.run_aerobuildup_custom(
            left_aileron=-dcontrol_deg,
            right_aileron=dcontrol_deg,
        )
        minus = self.run_aerobuildup_custom(
            left_aileron=dcontrol_deg,
            right_aileron=-dcontrol_deg,
        )

        CYda = (plus["CY"] - minus["CY"]) / (2 * d_control_rad)
        Clda = (plus["Cl"] - minus["Cl"]) / (2 * d_control_rad)
        Cnda = (plus["Cn"] - minus["Cn"]) / (2 * d_control_rad)

        # rudder derivatives
        plus = self.run_aerobuildup_custom(
            left_rudder=-dcontrol_deg,
            right_rudder=-dcontrol_deg,
        )
        minus = self.run_aerobuildup_custom(
            left_rudder=dcontrol_deg,
            right_rudder=dcontrol_deg,
        )

        CYdr = (plus["CY"] - minus["CY"]) / (2 * d_control_rad)
        Cldr = (plus["Cl"] - minus["Cl"]) / (2 * d_control_rad)
        Cndr = (plus["Cn"] - minus["Cn"]) / (2 * d_control_rad)

        return {
            "CXu": CXu,
            "CXa": CXa,
            "CXq": CXq,
            "CXde": CXde,

            "CZu": CZu,
            "CZa": CZa,
            "CZq": CZq,
            "CZde": CZde,

            "Cm": Cm_base,
            "Cmu": Cmu,
            "Cma": Cma,
            "Cmq": Cmq,
            "Cmde": Cmde,

            "CYb": CYb,
            "CYp": CYp,
            "CYr": CYr,
            "CYda": CYda,
            "CYdr": CYdr,

            "Clb": Clb,
            "Clp": Clp,
            "Clr": Clr,
            "Clda": Clda,
            "Cldr": Cldr,

            "Cnb": Cnb,
            "Cnp": Cnp,
            "Cnr": Cnr,
            "Cnda": Cnda,
            "Cndr": Cndr,
        }

    def Coefficients(self):
        """
        Returns coefficients in the exact order expected by:

        CXu, CXa, CXq, CXde,
        CZu, CZa, CZq, CZde,
        Cm, Cmu, Cma, Cmq, Cmde,
        CYb, CYp, CYr, CYda, CYdr,
        Clb, Clp, Clr, Clda, Cldr,
        Cnb, Cnp, Cnr, Cnda, Cndr
        """

        d = self.compute_all_derivatives()

        return (
            d["CXu"],
            d["CXa"],
            d["CXq"],
            d["CXde"],

            d["CZu"],
            d["CZa"],
            d["CZq"],
            d["CZde"],

            d["Cm"],
            d["Cmu"],
            d["Cma"],
            d["Cmq"],
            d["Cmde"],

            d["CYb"],
            d["CYp"],
            d["CYr"],
            d["CYda"],
            d["CYdr"],

            d["Clb"],
            d["Clp"],
            d["Clr"],
            d["Clda"],
            d["Cldr"],

            d["Cnb"],
            d["Cnp"],
            d["Cnr"],
            d["Cnda"],
            d["Cndr"],
        )

    def DynamicAnalysisInputs(self):
        """
        Returns the aircraft quantities typically required by
        the dynamic-analysis classes.
        """

        op_point = self.make_op_point()

        return {
            "op_point": op_point,
            "V0": self.velocity,
            "b": self.span,
            "c": self.c_ref_value(),
            "S": self.s_ref_value(),
            "rho": op_point.atmosphere.density(),
            "alpha_deg": self.alpha,
            "beta_deg": self.beta,
        }

    #test

    # def body_axis_from_wind_axis(self, aero):
    #     alpha = np.radians(self.alpha)
    #
    #     CL = aero["CL"]
    #     CD = aero["CD"]
    #
    #     CX = -CD * np.cos(alpha) + CL * np.sin(alpha)
    #     CZ = -CD * np.sin(alpha) - CL * np.cos(alpha)
    #
    #     # d/dalpha transform, approximate using returned CL/CD alpha derivatives
    #     CXa = (
    #             -aero["CDa"] * np.cos(alpha)
    #             + CD * np.sin(alpha)
    #             + aero["CLa"] * np.sin(alpha)
    #             + CL * np.cos(alpha)
    #     )
    #
    #     CZa = (
    #             -aero["CDa"] * np.sin(alpha)
    #             - CD * np.cos(alpha)
    #             - aero["CLa"] * np.cos(alpha)
    #             + CL * np.sin(alpha)
    #     )
    #
    #     CXq = -aero["CDq"] * np.cos(alpha) + aero["CLq"] * np.sin(alpha)
    #     CZq = -aero["CDq"] * np.sin(alpha) - aero["CLq"] * np.cos(alpha)
    #
    #     return CX, CZ, CXa, CZa, CXq, CZq
    #
    # def compute_all_derivatives_test(self):
    #     aero = self.run_aerobuildup()
    #
    #     CX, CZ, CXa, CZa, CXq, CZq = self.body_axis_from_wind_axis(aero)
    #
    #     derivatives = {
    #         "CXu": 0.0,  # AeroBuildup does not directly return velocity derivatives
    #         "CXa": CXa,
    #         "CXq": CXq,
    #         "CXde": 0.0,  # use manual finite difference if needed
    #
    #         "CZu": 0.0,
    #         "CZa": CZa,
    #         "CZq": CZq,
    #         "CZde": 0.0,
    #
    #         "Cm": aero["Cm"],
    #         "Cmu": 0.0,
    #         "Cma": aero["Cma"],
    #         "Cmq": aero["Cmq"],
    #         "Cmde": 0.0,
    #
    #         "CYb": aero["CYb"],
    #         "CYp": aero["CYp"],
    #         "CYr": aero["CYr"],
    #         "CYda": 0.0,
    #         "CYdr": 0.0,
    #
    #         "Clb": aero["Clb"],
    #         "Clp": aero["Clp"],
    #         "Clr": aero["Clr"],
    #         "Clda": 0.0,
    #         "Cldr": 0.0,
    #
    #         "Cnb": aero["Cnb"],
    #         "Cnp": aero["Cnp"],
    #         "Cnr": aero["Cnr"],
    #         "Cnda": 0.0,
    #         "Cndr": 0.0,
    #     }
    #
    #     return derivatives
    #
    # def scalarize(self, v):
    #     arr = np.asarray(v)
    #
    #     if arr.size == 1:
    #         return float(arr.reshape(-1)[0])
    #
    #     return arr
    #
    # def print_derivatives_as_text(self):
    #     derivatives = self.compute_all_derivatives()
    #
    #     print("\n--- AeroBuildup Stability Derivatives ---")
    #     for k, v in derivatives.items():
    #         v = self.scalarize(v)
    #
    #         if isinstance(v, float):
    #             print(f"{k:>8s}: {v:+.6e}")
    #         else:
    #             print(f"{k:>8s}: {v}")

    def Yaw_Check(self, T_eng=100, fraction_outer_engine=0.67):
        op_point = self.make_op_point()
        derivatives = self.compute_all_derivatives()

        Cndr = derivatives["Cndr"]
        Cnr = derivatives["Cnr"]

        # Engine lateral arm should usually be fraction of half-span, not full span
        y_eng = fraction_outer_engine * (self.span / 2)

        M_engine = T_eng * y_eng

        rho = op_point.atmosphere.density()
        q_dyn = 0.5 * rho * op_point.velocity ** 2

        # Use same S_ref as the aero model, not self.wing_area unless intentional
        Cn_OEI = M_engine / (q_dyn * self.s_ref_value() * self.span)

        deflection_OEI = -Cn_OEI / Cndr

        deflection_max = np.radians(30)

        r_max = (
                (Cndr / Cnr)
                * (deflection_max - deflection_OEI)
                * 2 * self.velocity / self.span
        )

        print("Cndr:", Cndr)
        print("Cnr:", Cnr)
        print("y_eng:", y_eng, "m")
        print("Cn_OEI:", Cn_OEI)
        print("OEI rudder deflection:", np.degrees(deflection_OEI), "deg")
        print("Remaining yaw rate:", np.degrees(r_max), "deg/s")

        return {
            "Cn_OEI": Cn_OEI,
            "deflection_OEI_rad": deflection_OEI,
            "deflection_OEI_deg": np.degrees(deflection_OEI),
            "r_max_rad_s": r_max,
            "r_max_deg_s": np.degrees(r_max),
        }

if __name__ == "__main__":
    aircraft = FlyingWingWithWingletsAeroBuildup()

    aircraft.make_airplane().draw()

    derivatives = aircraft.compute_all_derivatives()

    for k, v in derivatives.items():
        print(
            k,
            type(v),
            np.asarray(v).shape,
            v
        )

    print("\n--- Requested Stability and Control Derivatives ---")
    for k, v in derivatives.items():
        print(f"{k:>8s}: {v}")

    print(aircraft.Yaw_Check())