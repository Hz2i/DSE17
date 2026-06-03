import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox import Airplane
from aerosandbox.aerodynamics.aero_3D.test_aero_3D.geometries.conventional import wing_airfoil, tail_airfoil
from matplotlib import pyplot as plt
from numpy.ma.core import shape


class Control_Surface_Sizing():
    def __init__(self):
        self.coeff = None
        self.airplane = None
        self.wing_airfoil = asb.Airfoil("e344")
        self.tail_airfoil = asb.Airfoil("naca0010")
        self.op_point = asb.OperatingPoint(
            velocity=27.94,  # m/s
            alpha=7,  # degree
        )
        self.wing_sweep = 0.2618
        self.b = 30.08
        self.c = 1.203
        self.dihedral = 0.0
        self.start_inner_elevon = (self.b/2)*0.1
        self.elevon_connection = (self.b/2)*0.5
        self.end_outer_elevon = (self.b/2)*0.9
        self.half_span = self.b/2
        self.height_winglet = 1.5
        #self.control_requirements = control_requirements

    def Airplane_Geo(self, delta_inner=0, delta_outer=0):
        ### Define the 3D geometry you want to analyze/optimize.
        # Here, all distances are in meters and all angles are in degrees.
        self.airplane = asb.Airplane(
            name="AHAPS",
            xyz_ref=[1.8, 0, 0],  # CG location
            wings=[
                asb.Wing(
                    name="Main Wing",
                    symmetric=True,  # Should this wing be mirrored across the XZ plane?
                    xsecs=[
                        # The wing's cross ("X") sections
                        asb.WingXSec(  # Root
                            xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                            chord=self.c,  # Chord length at the XSec
                            twist=0,  # degrees
                            airfoil=self.wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[np.tan(self.wing_sweep) * self.start_inner_elevon, self.start_inner_elevon, self.dihedral], #[sweep, , dihedral]
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="inner_elevon",
                                    hinge_point=0.75,
                                    deflection=0.0,
                                    trailing_edge=True,
                                    symmetric=True,
                                ),
                            ],
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[np.tan(self.wing_sweep) * self.elevon_connection, self.elevon_connection, self.dihedral], #[sweep, , dihedral]
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="outer_elevon",
                                    hinge_point=0.75,
                                    deflection=0.0,
                                    trailing_edge=True,
                                    symmetric=False,
                                ),
                            ],
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[np.tan(self.wing_sweep) * self.end_outer_elevon, self.end_outer_elevon, self.dihedral], #[sweep, , dihedral]
                            chord=self.c,
                            twist=0,
                            airfoil=self.wing_airfoil,
                        ),

                        asb.WingXSec(  # Tip
                            xyz_le=[np.tan(self.wing_sweep) * self.half_span, self.half_span, self.dihedral],
                            chord=self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="rudder",
                                    hinge_point=0.75,
                                    deflection=0.0,
                                    trailing_edge=True,
                                    symmetric=True,
                                ),
                            ],
                        ),

                        asb.WingXSec(  # Tip
                            xyz_le=[(np.tan(self.wing_sweep) * self.half_span)+(np.tan(self.wing_sweep) * self.height_winglet), self.half_span, self.height_winglet],
                            chord=self.c,
                            twist=0,
                            airfoil=self.tail_airfoil,
                        ),
                    ],
                ),
            
                    # control_surfaces=[
                    #     asb.ControlSurface(
                    #         name='Inboard Elevon',
                    #         xsec_start=0.1,  # Start at 10% of the wing span
                    #         xsec_end=0.25,     # End at 40% of the wing span
                    #         hinge_line_xsec=0.75,  # Hinge line at 75% chord
                    #         deflection=10,
                    #     ),
                    #     asb.ControlSurface(
                    #         name='Outboard Elevon',
                    #         xsec_start=0.5,  # Start at 50% of the wing span
                    #         xsec_end=0.9,     # End at 90% of the wing span
                    #         hinge_line_xsec=0.75,  # Hinge line at 75% chord
                    #         deflection=10,
                    #     ),
                    # ],

                # asb.Wing(
                #     name="Horizontal Stabilizer",
                #     symmetric=True,
                #     xsecs=[
                #         asb.WingXSec(  # root
                #             xyz_le=[0, 0, 0],
                #             chord=0.1,
                #             twist=-10,
                #             airfoil=tail_airfoil,
                #         ),
                #         asb.WingXSec(  # tip
                #             xyz_le=[0.02, 0.17, 0], chord=0.08, twist=-10, airfoil=tail_airfoil
                #         ),
                #     ],
                # ).translate([0.6, 0, 0.06]),

            #     asb.Wing(
            #         name="Vertical Stabilizer",
            #         symmetric=False,
            #         xsecs=[
            #             asb.WingXSec(
            #                 xyz_le=[0, 0, 0],
            #                 chord=0.1,
            #                 twist=0,
            #                 airfoil=tail_airfoil,
            #             ),
            #             asb.WingXSec(
            #                 xyz_le=[0.04, 0, 0.15], chord=0.06, twist=0, airfoil=tail_airfoil
            #             ),
            #         ],
            #     ).translate([0.6, 0, 0.07]),
            # ],

        #     fuselages=[
        #         asb.Fuselage(
        #             name="Fuselage",
        #             xsecs=[
        #                 asb.FuselageXSec(
        #                     xyz_c=[0.8 * xi - 0.1, 0, 0.1 * xi - 0.03],
        #                     radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi),
        #                 )
        #                 for xi in np.cosspace(0, 1, 30)
        #             ],
        #         )
            ],
        )
        # self.airplane.with_control_deflections({"inner_elevon": 15.0})
        return self.airplane

    def vlm_run(self, delta_inner=0, delta_outer=0):
        
        # Rebuild airplane with new deflections
        self.Airplane_Geo(delta_inner, delta_outer)

        print(f"\nRequested deflections: inner={delta_inner}°, outer={delta_outer}°")

        def format_aero_value(v):
            try:
                arr = np.asarray(v)
            except Exception:
                return v
            if arr.shape == ():
                return float(arr)
            if arr.size == 1:
                return float(arr.flat[0])
            return arr.tolist() if arr.ndim <= 2 else arr

        # ===== AeroBuildup Analysis (DOES model control-surface effects) =====
        deflected_airplane = self.airplane.with_control_deflections({
            "inner_elevon": delta_inner,
            "outer_elevon": delta_outer,
        })
        
        try:
            ab = asb.AeroBuildup(
                airplane=deflected_airplane,
                op_point=self.op_point,
            )
            aero_ab = ab.run_with_stability_derivatives()
            cm_ab = aero_ab.get("Cm", None)
            if cm_ab is not None:
                cm_ab_scalar = format_aero_value(cm_ab)
                print(f"AeroBuildup Cm = {cm_ab_scalar:.6f}")
            else:
                print("AeroBuildup Cm = N/A")

            print("\nAeroBuildup results:")
            for k, v in aero_ab.items():
                try:
                    print(f"  {k} : {format_aero_value(v)}")
                except Exception:
                    print(f"  {k} : {v}")
        except Exception:
            import traceback
            print(f"AeroBuildup failed with error:")
            traceback.print_exc()
            aero_ab = None

        self.coeff = aero_ab
        return aero_ab

    def Control_Coefficients(self):
        d_deflect = 1
        deflection_points = np.arange(-20, 20 + d_deflect, d_deflect)
        Cm_list = []
        Cl_list = []
        Cn_list = []

        for i in deflection_points:
            self.op_point = asb.OperatingPoint(velocity=27.94, alpha=7)
            self.vlm_run(delta_inner=i, delta_outer=0)
            if self.coeff is not None:
                Cm_list.append(self.coeff.get("Cm", None))
            else:
                Cm_list.append(None)

        for i in deflection_points:
            self.op_point = asb.OperatingPoint(velocity=10, alpha=7)
            self.vlm_run(delta_inner=0, delta_outer=i)
            if self.coeff is not None:
                Cl_list.append(self.coeff.get("Cl", None))
            else:
                Cl_list.append(None)

        # print(Cm_list)
        plt.plot(deflection_points, Cm_list, color="red", label="Cm")
        plt.plot(deflection_points, Cl_list, color="blue", label="Cl")
        plt.xlabel("Elevon deflection (deg)")
        plt.ylabel("Cm")
        plt.title("Control coefficient sweep")
        plt.grid(True)
        plt.show()



#    def Control_Check(self):

if __name__ == "__main__":
    print("Starting simulation")
    cs = Control_Surface_Sizing()
    cs.Airplane_Geo(0, 0)
    cs.airplane.draw()
    # cs.vlm_run(-20, 0)
    # cs.Airplane_Geo()
    # cs.vlm_run()
    # cs.Control_Coefficients()
