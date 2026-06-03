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

    def Airplane_Geo(self, delta_inner=0, delta_outer=0, rudder=0):
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
            ],
        )
        return self.airplane

    def vlm_run(self, delta_inner=0, delta_outer=0, rudder=0):
        
        # Rebuild airplane with new deflections
        self.Airplane_Geo(delta_inner, delta_outer, rudder)

        print(f"\nRequested deflections: inner={delta_inner}°, outer={delta_outer}°, rudder={rudder}°")

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
            "rudder": rudder
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
            self.op_point = asb.OperatingPoint(velocity=10, alpha=7)
            self.vlm_run(delta_inner=i, delta_outer=0, rudder=0)
            if self.coeff is not None:
                Cm_list.append(self.coeff.get("Cm", None))
            else:
                Cm_list.append(None)

        for i in deflection_points:
            self.op_point = asb.OperatingPoint(velocity=10, alpha=7)
            self.vlm_run(delta_inner=0, delta_outer=i, rudder=0)
            if self.coeff is not None:
                Cl_list.append(self.coeff.get("Cl", None))
            else:
                Cl_list.append(None)

        for i in deflection_points:
            self.op_point = asb.OperatingPoint(velocity=10, alpha=7)
            self.vlm_run(delta_inner=0, delta_outer=0, rudder=i)
            if self.coeff is not None:
                Cn_list.append(self.coeff.get("Cn", None))
            else:
                Cn_list.append(None)

        # print(Cm_list)
        plt.plot(deflection_points, Cm_list, color="red", label="Cm")
        plt.plot(deflection_points, Cl_list, color="blue", label="Cl")
        plt.plot(deflection_points, Cn_list, color="green", label="Cn")
        plt.xlabel("Elevon deflection (deg)")
        plt.ylabel("Cm")
        plt.title("Control coefficient sweep")
        plt.grid(True)
        plt.show()



#    def Control_Check(self):

if __name__ == "__main__":
    print("Starting simulation")
    cs = Control_Surface_Sizing()
    # cs.Airplane_Geo(0, 0, 0)
    # cs.airplane.draw()
    # cs.vlm_run(-20, 0)
    # cs.Airplane_Geo()
    # cs.vlm_run()
    cs.Control_Coefficients()
