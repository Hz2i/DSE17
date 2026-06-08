import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib
from fontTools.feaLib import error
from scipy.sparse.linalg import eigen

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class Control_Surface_Sizing():
    def __init__(self):
        self.coeff = None
        self.airplane = None
        self.wing_airfoil = asb.Airfoil("mh60")
        self.tail_airfoil = asb.Airfoil("naca0012")

        self.wing_sweep = 0.2618      # radians
        self.b = 28.80                # full span [m]
        self.c = 1.44                 # chord [m]
        self.S = self.b * self.c      # Wing area [m^2]
        self.dihedral = np.radians(1.0)

        self.x_cg = 2.5

        self.inner_elevon_frac = 0.05
        self.outer_elevon_frac = 0.20
        self.height_winglet = 1.0 # height of winglet above main wing [m]
        self.rudder_frac = 0.8
        self.fraction_outer_engine = None

        self.half_span = self.b / 2
        self.start_inner_elevon = None
        self.elevon_connection   = None
        self.end_outer_elevon    = None

        # Operating point (can be overridden before calling vlm_run)
        self.op_point = asb.OperatingPoint(velocity=27.94, alpha=7)

        self.q_check = True
        self.p_check = True
        self.r_check = True

        self.d_deflect = 5
        self.deflection_points = np.arange(-25, 25 + self.d_deflect, self.d_deflect) # don't use floats with arange
        self.print_plots = False

        self.T_eng = None
        self.y_eng = None

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------
    def Airplane_Geo(self, delta_inner=0, delta_outer=0, delta_rudder=0, outer_symmetric=True):
        """
        Build the airplane geometry.

        Parameters
        ----------
        delta_inner     : inner elevon deflection [deg]
        delta_outer     : outer elevon deflection [deg]
        delta_rudder    : rudder deflection [deg]
        outer_symmetric : True  → both outer panels deflect equally (pitch)
                          False → panels deflect opposite (roll / aileron)
        """


        self.start_inner_elevon = self.half_span * (0.9 - self.inner_elevon_frac - self.outer_elevon_frac)
        self.elevon_connection = self.half_span * (0.9 - self.outer_elevon_frac)
        self.end_outer_elevon = self.half_span * 0.9

        self.airplane = asb.Airplane(
            name="AHAPS",
            xyz_ref=[self.x_cg, 0, 0],   # CG at ~37% chord — adjust to your actual CG
            wings=[
                asb.Wing(
                    name="Main Wing",
                    symmetric=True,
                    xsecs=[
                        # Root
                        asb.WingXSec(
                            xyz_le=[0, 0, 0],
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                        ),

                        # Inner elevon starts here (10 % semi-span)
                        asb.WingXSec(
                            xyz_le=[
                                np.tan(self.wing_sweep) * self.start_inner_elevon,
                                self.start_inner_elevon,
                                np.tan(self.dihedral)*self.start_inner_elevon,
                            ],
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="inner_elevon",
                                    hinge_point=0.75,
                                    deflection=delta_inner,
                                    trailing_edge=True,
                                    symmetric=True,         # always symmetric for pitch
                                ),
                            ],
                        ),

                        # Transition — inner elevon ends / outer elevon starts (50 % semi-span)
                        asb.WingXSec(
                            xyz_le=[
                                np.tan(self.wing_sweep) * self.elevon_connection,
                                self.elevon_connection,
                                np.tan(self.dihedral)*self.elevon_connection,
                            ],
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="outer_elevon",
                                    hinge_point=0.75,
                                    deflection=delta_outer,
                                    trailing_edge=True,
                                    symmetric=outer_symmetric,  # caller decides
                                ),
                            ],
                        ),

                        # Outer elevon ends (90 % semi-span)
                        asb.WingXSec(
                            xyz_le=[
                                np.tan(self.wing_sweep) * self.end_outer_elevon,
                                self.end_outer_elevon,
                                np.tan(self.dihedral)*self.end_outer_elevon,
                            ],
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                        ),

                        # Tip
                        asb.WingXSec(
                            xyz_le=[
                                np.tan(self.wing_sweep) * self.half_span,
                                self.half_span,
                                np.tan(self.dihedral)*self.half_span,
                            ],
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                        ),
                    ],
                ),
                               # ================= LEFT WINGLET =================
                asb.Wing(
                    name="Left Winglet",
                    symmetric=True,
                    xsecs=[
                        asb.WingXSec(
                            xyz_le=[(np.tan(self.wing_sweep) * self.half_span)+(np.tan(self.wing_sweep) * self.height_winglet), self.half_span, self.height_winglet+(np.tan(self.dihedral)*self.half_span)],
                            chord=self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                        ),

                        asb.WingXSec(
                            xyz_le=[(np.tan(self.wing_sweep) * self.half_span)+(np.tan(self.wing_sweep)*(self.rudder_frac * self.height_winglet)), self.half_span, self.rudder_frac * self.height_winglet+(np.tan(self.dihedral)*self.half_span)],
                            chord=self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="rudder_left",
                                    hinge_point=0.8,
                                    deflection=delta_rudder,
                                    trailing_edge=True,
                                    symmetric=False,
                                )
                            ],
                        ),
                        asb.WingXSec(
                            xyz_le=[np.tan(self.wing_sweep) * self.half_span, self.half_span, np.tan(self.dihedral)*self.half_span],
                            chord=self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                        ),
                    ],
                ),
            ],
        )
        return self.airplane

    # ------------------------------------------------------------------
    # VLM / AeroBuildup runner
    # ------------------------------------------------------------------
    def vlm_run(self, delta_inner=0, delta_outer=0, delta_rudder=0, outer_symmetric=True, verbose=False):
        """
        Rebuild the geometry and run AeroBuildup.

        Returns the aero dict (also stored in self.coeff).
        """
        self.Airplane_Geo(delta_inner, delta_outer, delta_rudder, outer_symmetric)

        if verbose:
            print(f"\ndeflections: inner={delta_inner}°  outer={delta_outer}°  rudder={delta_rudder}"
                  f"outer_symmetric={outer_symmetric}")

        deflected = self.airplane.with_control_deflections({
            "inner_elevon": delta_inner,
            "outer_elevon": delta_outer,
            "rudder_left": delta_rudder,
        })

        try:
            ab = asb.AeroBuildup(
                airplane=deflected,
                op_point=self.op_point,
            )
            aero = ab.run_with_stability_derivatives()
        except Exception:
            import traceback
            print("AeroBuildup failed:")
            traceback.print_exc()
            aero = None

        self.coeff = aero
        return aero

    # ------------------------------------------------------------------
    # Sweep helper
    # ------------------------------------------------------------------
    def _sweep(self, deflection_points, delta_inner_fn, delta_outer_fn, delta_rudder_fn,
               outer_symmetric, coeff_key):
        """
        Generic deflection sweep.  Returns a list of scalar coefficient values.
        """
        results = []
        for i in deflection_points:
            aero = self.vlm_run(
                delta_inner=delta_inner_fn(i),
                delta_outer=delta_outer_fn(i),
                delta_rudder=delta_rudder_fn(i),
                outer_symmetric=outer_symmetric,
            )
            val = aero.get(coeff_key, None) if aero is not None else None
            # Flatten to scalar
            if val is not None:
                try:
                    val = float(np.asarray(val).flat[0])
                except Exception:
                    val = None
            results.append(val)
        return results

    # ------------------------------------------------------------------
    # Sweep helper (single run) | Only for damping coefficients!!! |
    # ------------------------------------------------------------------
    def _sweep_single(self, deflection_points, delta_inner_fn, delta_outer_fn, delta_rudder_fn,
               outer_symmetric, coeff_key):
        """
        Generic deflection sweep.  Returns a list of scalar coefficient values.
        """
        aero = self.vlm_run(
            delta_inner=delta_inner_fn(0),
            delta_outer=delta_outer_fn(0),
            delta_rudder=delta_rudder_fn(0),
            outer_symmetric=outer_symmetric
        )
        val = aero.get(coeff_key, None) if aero is not None else None
        # Flatten to scalar
        if val is not None:
            try:
                val = float(np.asarray(val).flat[0])
            except Exception:
                val = None
        return val

    # ------------------------------------------------------------------
    # Main analysis
    # ------------------------------------------------------------------
    def Pitching_Coefficients(self, print_plots=False):
        """
        Three sweeps:
          1. Inner elevon (symmetric) → Cm — pitch authority
          2. Outer elevon (antisymmetric) → Cl — roll authority
          3. Rudder → Cn — yaw authority
        """
        print("-------------------------------")
        print("Running inner elevon Cm sweep …")
        Cm_inner = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: i,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=True,
            coeff_key="Cm",
        )
        x_np_inner = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: i,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=True,
            coeff_key="x_np",
        )

        if print_plots:
            # ── Plot ──────────────────────────────────────────────────────
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle("Elevon Control Effectiveness  (V=" + str(self.op_point.velocity) + "m/s, α=" + str(
                self.op_point.alpha) + "°)")

            # Left: pitch sweeps
            axes[0].plot(self.deflection_points, Cm_inner, color="red", label="Cm — inner elevon (sym)")
            axes[0].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[0].set_xlabel("Elevon deflection [deg]")
            axes[0].set_ylabel("Cm")
            axes[0].set_title("Pitch authority")
            # axes[0, 0].legend()
            axes[0].grid(True)

            axes[1].plot(self.deflection_points, x_np_inner)
            axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[1].set_xlabel("Elevon deflection [deg]")
            axes[1].set_ylabel("x_np")
            axes[1].set_title("x_np for pitch")
            # axes[1, 0].legend()
            axes[1].grid(True)

            plt.tight_layout()
            plt.show()

        # ── Control Coefficient Linearisation ─────────────────────────
        Cmde = (Cm_inner[np.size(self.deflection_points) - 1] - Cm_inner[0]) / np.radians(self.deflection_points[np.size(self.deflection_points) - 1] - self.deflection_points[0])

        print("Running Cmq sweep …")
        Cmq = self._sweep_single(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cmq",
        )
        print("Cmq:", Cmq, "/rad")
        print("Cmde:", Cmde, "/rad")

        return Cmde, Cmq

    def Rolling_Coefficients(self):
        print("-------------------------------")
        print("Running outer elevon Cl sweep …")
        Cl_outer = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: i,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cl",
        )
        x_np_outer = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: i,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=True,
            coeff_key="x_np",
        )

        if self.print_plots:
            # ── Plot ──────────────────────────────────────────────────────
            fig, axes = plt.subplots(2, 3, figsize=(12, 5))
            fig.suptitle("Elevon Control Effectiveness  (V=" + str(self.op_point.velocity) + "m/s, α=" + str(self.op_point.alpha) + "°)")

            # Right: roll sweep
            axes[0].plot(self.deflection_points, Cl_outer, color="blue", label="Cl — outer elevon (antisym)")
            axes[0].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[0].set_xlabel("Elevon deflection [deg]")
            axes[0].set_ylabel("Cl")
            axes[0].set_title("Roll authority")
            # axes[0, 1].legend()
            axes[0].grid(True)

            axes[1].plot(self.deflection_points, x_np_outer)
            axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[1].set_xlabel("Elevon deflection [deg]")
            axes[1].set_ylabel("x_np")
            axes[1].set_title("x_np for roll")
            # axes[1, 1].legend()
            axes[1].grid(True)

            plt.tight_layout()
            plt.show()

        # ── Control Coefficient Linearisation ─────────────────────────
        Clda = (Cl_outer[np.size(self.deflection_points) - 1] - Cl_outer[0]) / np.radians(self.deflection_points[np.size(self.deflection_points) - 1] - self.deflection_points[0])

        print("Running Clp sweep …")
        Clp = self._sweep_single(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Clp",
        )
        print("Clp:", Clp, "/rad")
        print("Clda:", Clda, "/rad")

        return Clda, Clp

    def Yawing_Coefficients(self):
        print("-------------------------------")
        print("Running rudder Cn sweep (antisymmetric / yaw) …")
        Cn_rudder = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: i,
            outer_symmetric=False,
            coeff_key="Cn",
        )
        x_np_rudder = self._sweep(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: i,
            outer_symmetric=True,
            coeff_key="x_np",
        )
        # print(min(x_np_inner), min(x_np_outer), min(x_np_rudder))

        if self.print_plots:
            # ── Plot ──────────────────────────────────────────────────────
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle("Elevon Control Effectiveness  (V=" + str(self.op_point.velocity) + "m/s, α=" + str(self.op_point.alpha) + "°)")

            # Right: yaw sweep
            axes[0].plot(self.deflection_points, Cn_rudder, color="green", label="Cn — rudder (antisym)")
            axes[0].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[0].set_xlabel("Elevon deflection [deg]")
            axes[0].set_ylabel("Cn")
            axes[0].set_title("Yaw authority")
            # axes[0, 2].legend()
            axes[0].grid(True)

            axes[1].plot(self.deflection_points, x_np_rudder)
            axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
            axes[1].set_xlabel("Elevon deflection [deg]")
            axes[1].set_ylabel("x_np")
            axes[1].set_title("x_np for yaw")
            # axes[1, 2].legend()
            axes[1].grid(True)

            plt.tight_layout()
            plt.show()

        # ── Control Coefficient Linearisation ─────────────────────────
        Cndr = (Cn_rudder[np.size(self.deflection_points)-1] - Cn_rudder[0])/np.radians(self.deflection_points[np.size(self.deflection_points)-1] - self.deflection_points[0])

        print("Running Cmq sweep …")
        Cnr = self._sweep_single(
            self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cnr",
        )
        print("Cnr:", Cnr, "/rad")
        print("Cndr:", Cndr, "/rad")

        return Cndr, Cnr

    # ------------------------------------------------------------------
    # Control requirements check (placeholder)
    # ------------------------------------------------------------------
    def Pitch_Check(self):
        Cmde, Cmq = self.Pitching_Coefficients()
        q_req = np.radians(7)  # pitch rate [rad/s]

        q = Cmde/Cmq*(np.radians(self.deflection_points[np.size(self.deflection_points)-1]))*self.op_point.velocity/self.c
        print("Q:", np.degrees(q), "deg/s")

        return q, q_req

    def Roll_Check(self):
        Clda, Clp = self.Rolling_Coefficients()
        p_req = np.radians(10)  # roll  rate [rad/s]

        p = Clda / Clp * (np.radians(
            self.deflection_points[np.size(self.deflection_points) - 1])) * 2 * self.op_point.velocity / self.b
        print("P:", np.degrees(p), "deg/s")

        return p, p_req

    def Yaw_Check(self, T_eng=17, fraction_outer_engine=0.6):
        Cndr, Cnr = self.Yawing_Coefficients()
        r_req = np.radians(5)   # yaw  rate, now defined nose to left [rad/s]

        #OEI
        y_eng=fraction_outer_engine*self.half_span
        M_engine = T_eng*y_eng
        rho_cruise = 0.116
        k = 2
        Cn_OEI_counter = k*M_engine/(0.5*rho_cruise*self.op_point.velocity**2*self.S*self.b) # Required Cn to counteract OEI yawing moment
        deflection_OEI = -Cn_OEI_counter/Cndr
        print("OEI Deflection", np.degrees(deflection_OEI), "deg/s")
        deflection_max = np.radians(30)

        r = Cndr / Cnr * (deflection_max-deflection_OEI) * 2 * self.op_point.velocity / self.b
        # r = Cndr / Cnr * (np.radians(self.deflection_points[np.size(self.deflection_points) - 1])) * 2 * self.op_point.velocity / self.b
        print("R:", np.degrees(r), "deg/s")

        return r, r_req

    def Spiral_Check(self):
        aero = self.vlm_run(
            delta_inner=0, delta_outer=0, delta_rudder=0,
            outer_symmetric=False,
        )
        Clb = float(np.asarray(aero["Clb"]).flat[0])
        Cnb = float(np.asarray(aero["Cnb"]).flat[0])
        Clr = float(np.asarray(aero["Clr"]).flat[0])
        Cnr = float(np.asarray(aero["Cnr"]).flat[0])

        criterion = Clb * Cnr - Cnb * Clr

        print(f"  Cl_beta = {Clb:.5f}")
        print(f"  Cn_beta = {Cnb:.5f}")
        print(f"  Cl_r    = {Clr:.5f}")
        print(f"  Cn_r    = {Cnr:.5f}")
        print(f"  Spiral criterion (Cl_beta*Cn_r - Cn_beta*Cl_r) = {criterion:.6f}")

        if criterion > 0:
            print("Spiral mode is stable")
        else:
            print("Spiral mode is unstable")

        return criterion

    def Control_Sizing(self):
        p, p_req = self.Roll_Check()
        d_size_aileron = 0.1
        if p > p_req:
            while self.p_check:
                print("Aileron fraction:", self.outer_elevon_frac)
                self.outer_elevon_frac -= d_size_aileron

                if self.outer_elevon_frac < 0.0001:
                    self.outer_elevon_frac += d_size_aileron/2
                    d_size_aileron = d_size_aileron/2

                p, p_req = self.Roll_Check()

                if p < p_req:
                    if d_size_aileron > 0.001:
                        print("Aileron fraction failed guess:", self.outer_elevon_frac)
                        self.outer_elevon_frac += d_size_aileron #!!!!!
                        d_size_aileron = d_size_aileron/2
                    else:
                        self.p_check = False
                        self.outer_elevon_frac += d_size_aileron
                        p, p_req = self.Roll_Check()
                        print("----Final aileron fraction:", self.outer_elevon_frac)
                        print("----Final roll rate:", np.degrees(p), "deg/s")
                        # cs.airplane.draw()
        elif p < p_req:
            while self.p_check:
                print("Aileron fraction:", self.outer_elevon_frac)
                self.outer_elevon_frac += d_size_aileron
                p, p_req = self.Roll_Check()
                if p > p_req:
                    if d_size_aileron > 0.001:
                        print("Aileron fraction failed guess:", self.outer_elevon_frac)
                        self.outer_elevon_frac -= d_size_aileron #!!!!!
                        d_size_aileron = d_size_aileron / 2
                    else:
                        self.p_check = False
                        # self.outer_elevon_frac -= d_size_aileron
                        p, p_req = self.Roll_Check()
                        print("----Final aileron fraction:", self.outer_elevon_frac)
                        print("----Final roll rate:", np.degrees(p), "deg/s")
                        # cs.airplane.draw()


        q, q_req = self.Pitch_Check()
        d_size_elevator = 0.1
        if q > q_req:
            while self.q_check:
                print("Elevator fraction", self.inner_elevon_frac)
                self.inner_elevon_frac -= d_size_elevator

                if self.inner_elevon_frac < 0.0001:
                    self.inner_elevon_frac += d_size_elevator/2
                    d_size_elevator = d_size_elevator/2

                q, q_req = self.Pitch_Check()
                if q < q_req:
                    if d_size_elevator > 0.001:
                        print("Elevator fraction failed guess:", self.inner_elevon_frac)
                        self.inner_elevon_frac += d_size_elevator #!!!!!
                        d_size_elevator = d_size_elevator/ 2
                    else:
                        self.q_check = False
                        self.inner_elevon_frac += d_size_elevator
                        q, q_req = self.Pitch_Check()
                        print("----Final elevator fraction:", self.inner_elevon_frac)
                        print("----Final pitch rate:", np.degrees(q), "deg/s")
                        #cs.airplane.draw()

        elif q < q_req:
            while self.q_check:
                print("Elevator fraction", self.inner_elevon_frac)
                self.inner_elevon_frac += d_size_elevator
                q, q_req = self.Pitch_Check()
                if q > q_req:
                    if d_size_elevator > 0.001:
                        print("Elevator fraction failed guess:", self.inner_elevon_frac)
                        self.inner_elevon_frac -= d_size_elevator #!!!!!
                        d_size_elevator = d_size_elevator/ 2
                    else:
                        self.q_check = False
                        # self.inner_elevon_frac -= d_size_elevator
                        q, q_req = self.Pitch_Check()
                        print("----Final elevator fraction:", self.inner_elevon_frac)
                        print("----Final pitch rate:", np.degrees(q), "deg/s")
                        #cs.airplane.draw()

        r, r_req = self.Yaw_Check()
        d_size_rudder = 0.1
        if r > r_req:
            while self.r_check:
                print("Rudder fraction:", self.rudder_frac)
                self.rudder_frac -= d_size_rudder

                if self.rudder_frac < 0.0001:
                    self.rudder_frac += d_size_rudder/2
                    d_size_rudder = d_size_rudder/2

                r, r_req = self.Yaw_Check()
                if r < r_req:
                    if d_size_rudder > 0.001:
                        print("Rudder fraction failed guess:", self.rudder_frac)
                        self.rudder_frac += d_size_rudder #!!!!!
                        d_size_rudder = d_size_rudder/ 2
                    else:
                        self.r_check = False
                        self.rudder_frac += d_size_rudder
                        r, r_req = self.Yaw_Check()
                        print("----Final rudder fraction:", self.rudder_frac)
                        print("----Final yaw rate:", np.degrees(r), "deg/s")
                        cs.airplane.draw()
        elif r < r_req:
            while self.r_check:
                print("Rudder fraction:", self.rudder_frac)
                self.rudder_frac += d_size_rudder
                r, r_req = self.Yaw_Check()
                if r > r_req:
                    if d_size_rudder > 0.001:
                        print("Rudder fraction failed guess:", self.rudder_frac)
                        self.rudder_frac -= d_size_rudder #!!!!!
                        d_size_rudder = d_size_rudder / 2
                    else:
                        self.r_check = False
                        # self.rudder_frac -= d_size_rudder
                        r, r_req = self.Yaw_Check()
                        print("----Final rudder fraction:", self.rudder_frac)
                        print("----Final yaw rate:", np.degrees(r), "deg/s")
                        cs.airplane.draw()

        print("------------Final Values--------------")
        print("Final aileron fraction:", self.outer_elevon_frac)
        print("----Final roll rate:", np.degrees(p), "deg/s")
        print("Final elevator fraction:", self.inner_elevon_frac)
        print("----Final pitch rate:", np.degrees(q), "deg/s")
        print("Final rudder fraction:", self.rudder_frac)
        print("----Final yaw rate:", np.degrees(r), "deg/s")

    def Cm_check(self):
        cg_range = np.linspace(0., 3., 30)
        cm_list = []
        cma_list =[]

        for i in cg_range:
            self.x_cg = i
            cs.Airplane_Geo()
            cm = cs._sweep_single(self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cm"
            )
            cm_list.append(cm)

            cma = cs._sweep_single(self.deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cma"
            )
            cma_list.append(cma)

        plt.plot(cg_range, cm_list)
        plt.plot(cg_range, cma_list)
        plt.show()





        # todo OEI
        # todo controllability at forward cg
        # todo dutch roll/spiral stability





if __name__ == "__main__":
    print("Starting simulation")
    cs = Control_Surface_Sizing()
    cs.Airplane_Geo()
    # cs.airplane.draw()
    # cs.Control_Check()
    # cs.Control_Sizing()
    # cs.Pitching_Coefficients(print_plots=True)
    # cs.Spiral_Check()
    cs.Cm_check()