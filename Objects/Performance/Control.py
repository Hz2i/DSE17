import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.aerodynamics.aero_3D.test_aero_3D.geometries.conventional import wing_airfoil, tail_airfoil


class Control_Surface_Sizing():
    def __init__(self):
        self.airplane = None
        self.wing_airfoil = asb.Airfoil("sd7037")
        self.tail_airfoil = asb.Airfoil("naca0010")
        self.op_point=asb.OperatingPoint(
                velocity=27.94,  # m/s
                alpha=10,  # degree
            )
        self.wing_sweep = 0.2618
        self.b = 30.08
        self.c = 1.203
        self.dihedral = 0.0
        self.start_inner_elevon = (self.b/2)*0.1
        self.elevon_connection = (self.b/2)*0.5
        self.end_outer_elevon = (self.b/2)*0.9
        self.half_span = self.b/2
        #self.control_requirements = control_requirements

    def Airplane_Geo(self, delta_inner=0, delta_outer=0):
        ### Define the 3D geometry you want to analyze/optimize.
        # Here, all distances are in meters and all angles are in degrees.
        self.airplane = asb.Airplane(
            name="AHAPS",
            xyz_ref=[0, 0, 0],  # CG location
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
                                    deflection=delta_inner,
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
                                    deflection=delta_outer,
                                    trailing_edge=True,
                                    symmetric=True,
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
                            airfoil=self.wing_airfoil,
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
        
        # ===== VLM Analysis (doesn't model control-surface effects) =====
        vlm = asb.VortexLatticeMethod(
            airplane=self.airplane,
            op_point=self.op_point,
        )
        aero_vlm = vlm.run_with_stability_derivatives()
        print("\nVLM results:")
        for k, v in aero_vlm.items():
            try:
                # Try printing scalars nicely
                if hasattr(v, "shape") and getattr(v, "shape") != ():
                    print(f"  {k} : array shape {np.asarray(v).shape}")
                else:
                    print(f"  {k} : {float(v):.6f}")
            except Exception:
                print(f"  {k} : {v}")
        
        # ===== AeroBuildup Analysis (DOES model control-surface effects) =====
        # Create deflected airplane copy for AeroBuildup
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
            cm_ab = aero_ab.get('Cm', None)
            if cm_ab is not None:
                # Handle numpy array: flatten and get first scalar
                cm_ab_scalar = float(np.asarray(cm_ab).flat[0])
                print(f"AeroBuildup Cm = {cm_ab_scalar:.6f}")
            else:
                print("AeroBuildup Cm = N/A")
            print("\nAeroBuildup results:")
            for k, v in aero_ab.items():
                try:
                    if hasattr(v, "shape") and getattr(v, "shape") != ():
                        print(f"  {k} : array shape {np.asarray(v).shape}")
                    else:
                        print(f"  {k} : {float(v):.6f}")
                except Exception:
                    print(f"  {k} : {v}")
        except Exception as e:
            import traceback
            print(f"AeroBuildup failed with error:")
            traceback.print_exc()
            aero_ab = None
    
        # Return both results for programmatic consumption
        return {"vlm": aero_vlm, "aerobuildup": aero_ab}

#    def Control_Check(self):
#       control_surface = Control_Surface_Sizing()
#       control_surface.Airplane_Geo()
#       control_surface.vlm_run()

if __name__ == "__main__":
    print("Starting simulation")
    cs = Control_Surface_Sizing()
    cs.Airplane_Geo(0, 0)
    cs.vlm_run(-20, -20)