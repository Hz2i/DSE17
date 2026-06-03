import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class Control_Surface_Sizing():
    def __init__(self):
        self.coeff = None
        self.airplane = None
        self.wing_airfoil = asb.Airfoil("mh91")
        self.tail_airfoil = asb.Airfoil("naca0012")

        self.wing_sweep = 0.2618      # radians
        self.b = 30.08                # full span [m]
        self.c = 1.203                # chord [m]
        self.dihedral = 0.0
        self.height_winglet = 4      # height of winglet above main wing [m]

        self.half_span = self.b / 2
        self.start_inner_elevon = self.half_span * 0.1
        self.elevon_connection   = self.half_span * 0.5
        self.end_outer_elevon    = self.half_span * 0.9

        # Operating point (can be overridden before calling vlm_run)
        self.op_point = asb.OperatingPoint(velocity=27.94, alpha=7)

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
        outer_symmetric : True  → both outer panels deflect equally (pitch)
                          False → panels deflect opposite (roll / aileron)
        """
        self.airplane = asb.Airplane(
            name="AHAPS",
            xyz_ref=[0.45, 0, 0],   # CG at ~37% chord — adjust to your actual CG
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
                                self.dihedral,
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
                                self.dihedral,
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
                                self.dihedral,
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
                                self.dihedral,
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
                    symmetric=False,
                    xsecs=[
                        asb.WingXSec(
                            xyz_le=[np.tan(self.wing_sweep) * self.half_span, self.half_span, 0],
                            chord=1.5*self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="rudder_left",
                                    hinge_point=0.6,
                                    deflection=0.0,
                                    trailing_edge=True,
                                    symmetric=False,
                                )
                            ],
                        ),
                        asb.WingXSec(
                            xyz_le=[(np.tan(self.wing_sweep) * self.half_span)+(np.tan(self.wing_sweep) * self.height_winglet), self.half_span, self.height_winglet],
                            chord=1.5*self.c,
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
    # Main analysis
    # ------------------------------------------------------------------
    def Control_Coefficients(self):
        """
        Three sweeps:
          1. Inner elevon (symmetric) → Cm   — pitch authority
          2. Outer elevon (antisymmetric) → Cl — roll authority
        """
        d_deflect = 1
        deflection_points = np.arange(-20, 20 + d_deflect, d_deflect)

        self.op_point = asb.OperatingPoint(velocity=27.94, alpha=7)

        print("Running inner elevon Cm sweep …")
        Cm_inner = self._sweep(
            deflection_points,
            delta_inner_fn=lambda i: i,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=True,
            coeff_key="Cm",
        )

        print("Running outer elevon Cl sweep (antisymmetric / roll) …")
        Cl_outer = self._sweep(
            deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: i,
            delta_rudder_fn=lambda i: 0,
            outer_symmetric=False,
            coeff_key="Cl",
        )

        print("Running rudder Cn sweep (antisymmetric / yaw) …")
        Cn_rudder = self._sweep(
            deflection_points,
            delta_inner_fn=lambda i: 0,
            delta_outer_fn=lambda i: 0,
            delta_rudder_fn=lambda i: i,
            outer_symmetric=False,
            coeff_key="Cn",
        )

        # ── Plot ──────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle("Elevon Control Effectiveness  (V=27.94 m/s, α=7°)")

        # Left: pitch sweeps
        axes[0].plot(deflection_points, Cm_inner, color="red",   label="Cm — inner elevon (sym)")
        axes[0].axhline(0, color="black", linewidth=0.8, linestyle="--")
        axes[0].set_xlabel("Elevon deflection [deg]")
        axes[0].set_ylabel("Cm")
        axes[0].set_title("Pitch authority")
        axes[0].legend()
        axes[0].grid(True)

        # Right: roll sweep
        axes[1].plot(deflection_points, Cl_outer, color="blue", label="Cl — outer elevon (antisym)")
        axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
        axes[1].set_xlabel("Elevon deflection [deg]")
        axes[1].set_ylabel("Cl")
        axes[1].set_title("Roll authority")
        axes[1].legend()
        axes[1].grid(True)

        # Right: yaw sweep
        axes[1].plot(deflection_points, Cl_outer, color="blue", label="Cn — rudder (antisym)")
        axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
        axes[1].set_xlabel("Elevon deflection [deg]")
        axes[1].set_ylabel("Cn")
        axes[1].set_title("Roll authority")
        axes[1].legend()
        axes[1].grid(True)

        plt.tight_layout()
        plt.show()

    # ------------------------------------------------------------------
    # Control requirements check (placeholder)
    # ------------------------------------------------------------------
    def Control_Check(self):
        q_req = np.radians(3)    # pitch rate [rad/s]
        p_req = np.radians(10)   # roll  rate [rad/s]
        r_req = np.radians(5)    # yaw   rate [rad/s]
        # todo


if __name__ == "__main__":
    print("Starting simulation")
    cs = Control_Surface_Sizing()
    cs.Airplane_Geo()
    cs.airplane_draw()
    # cs.Control_Coefficients()
