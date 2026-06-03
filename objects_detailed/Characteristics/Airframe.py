import numpy as np
import aerosandbox as asb
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from objects_detailed.Characteristics.ReferenceGeometries import *
import aerosandbox.tools.pretty_plots as p

class fuselage:
    def __init__(self,D=0.3,L1=0.5,L2=3,L3=1.5):
        self.D = D
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3

        self.m = None           # Currently initialised with None; Class 2 estimation methods required!
        self.v_total = None     # Currently initialised with None



class nacelles:
    def __init__(self, nr_of_engines = 4, diameter = 0.5, length = 1.0, pos = [-0.6, -0.3, 0.3, 0.6]):
        self.nr_of_engines = nr_of_engines
        self.positions = pos
        self.D = diameter
        self.L = length

        self.m = None

class airframe:
    def __init__(self, S=36.0, A=25.0, qc_sweep=0.0, taper=1.0, dihedral=0.0 , airfoil=asb.Airfoil("e335"), fus = fuselage(), nac = nacelles(), display=False):
        self.foil = airfoil
        self.AR = A
        self.taper = taper
        self.qc_sweep = qc_sweep        # RADIANS!!!
        self.le_sweep = np.arctan(np.tan(self.qc_sweep) + (1 - self.taper)/((1+self.taper)*self.AR))
        self.dihedral = dihedral        # RADIANS!!!
        self.S = S                
        self.CL_grad = None          # Currently initialised with None; Add method to compute
        self.CL_max = None           # Currently initialised with None; Add method to compute
        self.e = None
        self.x_ac = None

        self.m_wing = None                # Currently initialised with None; Class 2 estimation methods required!
        self.v_int_wing = None            # Currently initialised with None; Stress calculations and internal design required!
        self.v_tot_wing = None            # Currently initialised with None

        self.fuselage = fus             # Fuselage geometry and parameters
        self.nacelles = nac             # Fuselage Geometry and parameters

        self.define_geometry()
        self.display_plane(display)


    def define_geometry(self):
        self.le_sweep = np.arctan(np.tan(self.qc_sweep) + (1 - self.taper)/((1+self.taper)*self.AR))
        self.c_r = 2.0 * (self.S/self.AR) ** 0.5 / (1 + self.taper)
        self.b = (self.S * self.AR) ** 0.5
        self.c_t = self.c_r * self.taper

        self.geometry_asb = asb.Airplane(
            name="AHAPS",
            xyz_ref = [0.0, 0.0, 0.0],
            wings = [
                asb.Wing(
                    name="Main Planform",
                    symmetric = True,
                    xsecs=[
                        asb.WingXSec(
                            xyz_le=[0.0,
                                    0.0,
                                    0.0],
                            chord=self.c_r,
                            twist=0,
                            airfoil=self.foil),
                        asb.WingXSec(
                            xyz_le=[(self.b/2) * np.tan(self.le_sweep),
                                    self.b/2,
                                    (self.b/2) * np.tan(self.dihedral)],
                            chord=self.c_t,
                            twist=0,
                            airfoil=self.foil)
                            ]
                        )
                    ],
            fuselages=[
                asb.Fuselage(
                    name="Nacelle_"+str(pos),
                    xsecs=[
                        asb.FuselageXSec(
                            xyz_c=[self.c_t * xi + abs(pos)*self.b/2*np.tan(self.le_sweep), self.b/2*pos, pos*self.b/2*np.tan(self.dihedral)],
                            radius=self.c_t * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
                            )
                        for xi in np.linspace(0.0, self.c_t, 30)
                        ]
                    )
                for pos in self.nacelles.positions
                ]
            )

    def display_plane(self, display):
        if display:
            drawn_airplane = self.geometry_asb.deepcopy()
            drawn_airplane.wings = [
                w.subdivide_sections(15, np.linspace) for w in drawn_airplane.wings
            ]
            drawn_airplane.fuselages = [f.subdivide_sections(2) for f in drawn_airplane.fuselages]
            drawn_airplane.draw_three_view()
            # self.geometry_asb.draw_three_view()
            p.show_plot(dpi=600)


    def vlm_analysis(self):         # Call to AeroSandbox for VLM implementation
        pass


    def llt_analysis(self, series=False, alpha=5.0):
        op_point = asb.OperatingPoint(
            atmosphere=asb.Atmosphere(altitude=0),
            velocity=15.0
        )
        if series:
            N_runs = len(alpha)
            llt_op_pt = op_point.copy()
            llt_op_pt.alpha = alpha

            llt_batch = [
                asb.LiftingLine(
                    airplane=self.geometry_asb,
                    op_point=op,
                    xyz_ref=self.geometry_asb.xyz_ref
                    ).run()
                for op in llt_op_pt
                ]

            llt_results = {}
            for param in llt_batch[0].keys():
                llt_results[param] = np.array([result[param] for result in llt_batch])
            llt_results["alpha"] = alpha


        else:
            op = op_point.copy()
            op.alpha = alpha
            llt_results = asb.LiftingLine(
                            airplane=self.geometry_asb,
                            op_point=op,
                            xys_ref=self.geometry_asb.xyz_ref
                            ).run()

        return llt_results

