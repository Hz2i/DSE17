import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.aerodynamics.aero_3D.test_aero_3D.geometries.conventional import wing_airfoil, tail_airfoil


class Control_Surface_Sizing():
    def __init__(self):
        self.wing_airfoil = asb.Airfoil("sd7037")
        self.tail_airfoil = asb.Airfoil("naca0010")

    def Airplane_Geo(self):
        ### Define the 3D geometry you want to analyze/optimize.
        # Here, all distances are in meters and all angles are in degrees.
        airplane = asb.Airplane(
            name="Peter's Glider",
            xyz_ref=[0, 0, 0],  # CG location
            wings=[
                asb.Wing(
                    name="Main Wing",
                    symmetric=True,  # Should this wing be mirrored across the XZ plane?
                    xsecs=[
                        # The wing's cross ("X") sections
                        asb.WingXSec(  # Root
                            xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                            chord=0.12,
                            twist=0,  # degrees
                            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[0.00, 0.15, 0], #[sweep, , dihedral]
                            chord=0.12,
                            twist=0,
                            airfoil=wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="inner elevon",
                                    hinge_point=0.75,
                                    deflection=0.0,
                                    symmetric=True,
                                ),
                            ],
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[0.00, 0.35, 0], #[sweep, , dihedral]
                            chord=0.12,
                            twist=0,
                            airfoil=wing_airfoil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="outer elevon",
                                    hinge_point=0.75,
                                    deflection=0.0,
                                    symmetric=True,
                                ),
                            ],
                        ),

                        asb.WingXSec(  # Mid
                            xyz_le=[0.00, 0.55, 0], #[sweep, , dihedral]
                            chord=0.12,
                            twist=0,
                            airfoil=wing_airfoil,
                        ),

                        asb.WingXSec(  # Tip
                            xyz_le=[0.00, 0.6, 0],
                            chord=0.12,
                            twist=0,
                            airfoil=wing_airfoil,
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
        return self.airplane

    def vlm_run(self):
        vlm = asb.VortexLatticeMethod(
            airplane=self.airplane,
            op_point=asb.OperatingPoint(
                velocity=28,  # m/s
                alpha=0,  # degree
            ),
        )
        aero = vlm.run_with_stability_derivatives()  # Returns a dictionary
        for k, v in aero.items():
            print(f"{k.rjust(4)} : {v}")

        # vlm.draw(show_kwargs=dict(jupyter_backend="static"))

control_surface = Control_Surface_Sizing()
control_surface.Airplane_Geo()
control_surface.vlm_run()
