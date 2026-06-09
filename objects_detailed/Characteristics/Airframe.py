import numpy as np
import aerosandbox as asb
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from objects_detailed.Characteristics.ReferenceGeometries import *
import aerosandbox.tools.pretty_plots as p
from objects_detailed.ModifiedLibraries.lifting_line_ADJUSTED import LiftingLine as LLT_Adjusted

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
    def __init__(self, S=36.0, A=20.0, qc_sweep=15.0/180*np.pi, taper=1.0, dihedral=0.0 , airfoil=asb.Airfoil("mh91"), fus = fuselage(), nac = nacelles(), display=False, init_polar=True):
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
        if init_polar:
            self.compute_polar()
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
                    symmetric = False,
                    xsecs=[
                        # Left Half:
                        asb.WingXSec(
                            xyz_le=[(self.b/2) * np.tan(self.le_sweep),
                                    -self.b/2,
                                    (self.b/2) * np.tan(self.dihedral)],
                            chord=self.c_t,
                            twist=0,
                            airfoil=self.foil),
                        # Root Section::
                        asb.WingXSec(
                            xyz_le=[0.0,
                                    0.0,
                                    0.0],
                            chord=self.c_r,
                            twist=0,
                            airfoil=self.foil),
                        # Right Half:
                        asb.WingXSec(
                            xyz_le=[(self.b/2) * np.tan(self.le_sweep),
                                    self.b/2,
                                    (self.b/2) * np.tan(self.dihedral)],
                            chord=self.c_t,
                            twist=0,
                            airfoil=self.foil),  
                        ]
                        ),
                    asb.Wing(
                    name="Rudder",
                    symmetric = True,
                    xsecs=[
                        asb.WingXSec(
                            xyz_le=[(self.b/2) * np.tan(self.le_sweep),
                                    self.b/2,
                                    (self.b/2) * np.tan(self.dihedral)],
                            chord=self.c_t,
                            twist=0,
                            airfoil=self.foil),
                        asb.WingXSec(
                            xyz_le=[(self.b/2) * np.tan(self.le_sweep),
                                    self.b/2,
                                    (self.b/2) * np.tan(self.dihedral)+1.5],
                            chord=self.c_t,
                            twist=0,
                            airfoil=self.foil),
                        ]
                        ),
                    ],
            # fuselages=[
            #     asb.Fuselage(
            #         name="Nacelle_"+str(pos),
            #         xsecs=[
            #             asb.FuselageXSec(
            #                 xyz_c=[self.c_t * xi + abs(pos)*self.b/2*np.tan(self.le_sweep), self.b/2*pos, pos*self.b/2*np.tan(self.dihedral)],
            #                 radius=self.c_t * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
            #                 )
            #             for xi in np.linspace(0.0, self.c_t, 30)
            #             ]
            #         )
            #     for pos in self.nacelles.positions
            #     ]
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


    def llt_analysis(self, series=False, alpha=5.0, alt=18500.0, TAS=25.0, resolution=5):
        op_point = asb.OperatingPoint(
            atmosphere=asb.Atmosphere(altitude=alt),
            velocity=TAS
        )
        if series:
            N_runs = len(alpha)
            llt_op_pt = op_point.copy()
            llt_op_pt.alpha = alpha

            llt_batch = [
                LLT_Adjusted(
                    airplane=self.geometry_asb,
                    op_point=op,
                    xyz_ref=self.geometry_asb.xyz_ref,
                    spanwise_resolution=resolution
                    ).run()
                for op in llt_op_pt
                ]

            llt_results = {}
            for param in llt_batch[0].keys():
                llt_results[param] = np.array([result[param] for result in llt_batch])
            llt_results["alpha"] = alpha

            return llt_results


        else:
            op = op_point.copy()
            op.alpha = alpha
            llt_an = LLT_Adjusted(
                            airplane=self.geometry_asb,
                            op_point=op,
                            xyz_ref=self.geometry_asb.xyz_ref,
                            spanwise_resolution=resolution
                            )
            llt_results = llt_an.run()

            return llt_results, llt_an


    def compute_polar(self, alpha_range=np.linspace(-10.0, 20.0, 30), alt=18500.0, TAS=25.0, res=5):
        llt_data = self.llt_analysis(series=True, alpha=alpha_range, alt=alt, TAS=TAS, resolution=res)
        CL_data = llt_data["CL"]
        CD_data = llt_data["CD"]

        i_stall = np.argmax(CL_data)
        self.CL_max = CL_data[i_stall]

        CL_min = CL_data[np.argmin(CL_data)]
        self.CL_min = CL_min

        i_min = np.argmin(abs(CL_data - 0.2*CL_min))
        i_max = np.argmin(abs(CL_data - 0.95*self.CL_max))
        i_max_lift = np.argmin(abs(CL_data - 0.6*self.CL_max))

        lift_curve_coeff = np.polynomial.polynomial.polyfit(alpha_range[i_min:i_max_lift+1], CL_data[i_min:i_max_lift+1], 1)
        drag_polar_coeff = np.polynomial.polynomial.polyfit(CL_data[i_min:i_max+1], CD_data[i_min:i_max+1], 2)

        self.CL0 = lift_curve_coeff[0]
        self.CL_alpha = lift_curve_coeff[1]

        self.CD0 = drag_polar_coeff[0]
        self.K1 = drag_polar_coeff[1]
        self.K2 = drag_polar_coeff[2]

        CL_CD_data = CL_data/CD_data
        self.CL_CD_max = CL_CD_data[np.argmax(CL_CD_data)]


    def compute_load_distribution(self, alpha=5.0, alt=18500.0, TAS=25.0, res = 20):
        llt_data, llt_an = self.llt_analysis(series=False, alpha=alpha, alt=alt, TAS=TAS, resolution=res)

        total_F = llt_an.forces_inviscid_geometry + llt_an.forces_profile_geometry
        total_M = llt_an.moments_inviscid_geometry + llt_an.moments_profile_geometry

        self.vortex_coords = llt_an.vortex_centers

        panel_widths = abs(llt_an.front_right_vertices[:,1] - llt_an.front_left_vertices[:,1])

        self.panel_widths =panel_widths

        self.dFz_dy_current = total_F[:,2]/panel_widths
        self.dFx_dy_current = total_F[:,0]/panel_widths
        self.dFy_dy_current = total_F[:,1]/panel_widths

        self.dMx_dy_current = total_M[:,0]/panel_widths
        self.dMy_dy_current = total_M[:,1]/panel_widths
        self.dMz_dy_current = total_M[:,2]/panel_widths

        # return dFx_dy, dFy_dy, dFz_dy, dMx_dy, dMy_dy, dMz_dy, panel_widths, vortex_coords, llt_an
