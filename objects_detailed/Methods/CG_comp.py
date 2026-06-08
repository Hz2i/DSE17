import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Characteristics')))

import Airframe
import Components_Materials

import aerosandbox.numpy as np
import plotly.graph_objects as go
import aerosandbox as asb
import StructuralAnalysis as sa

import matplotlib.pyplot as plt

def Calculate_Airfoil_Centroid(airfoil, chord_length):

    # Datum: leading edge, x along chord, z upward

    upper_coords = airfoil.upper_coordinates()  # (x, z)
    lower_coords = airfoil.lower_coordinates()  # (x, z)

    # nondimensional x grid
    x = np.linspace(0, 1, 200)

    # interpolate surfaces (make sure ordering is correct)
    z_upper = np.interp(x, upper_coords[:, 0][::-1], upper_coords[:, 1][::-1])
    z_lower = np.interp(x, lower_coords[:, 0], lower_coords[:, 1])

    thickness = z_upper - z_lower

    # differential area
    dA = thickness

    # total area (nondimensional)
    A = np.trapezoid(dA, x)

    # centroid x (nondimensional)
    x_centroid = np.trapezoid(x * dA, x) / A

    # centroid z (nondimensional)
    z_mid = 0.5 * (z_upper + z_lower)
    z_centroid = np.trapezoid(z_mid * dA, x) / A

    # scale to physical coordinates
    x_centroid *= chord_length
    z_centroid *= chord_length   # assumes geometric similarity scaling

    return x_centroid, z_centroid


class CG_comp:
    def __init__(self, x_cg_goal=2.4, airframe=Airframe.airframe(), t_skin_airfoil=0.0002):
        self.airframe = airframe
        self.airfoil = airframe.foil
        self.x_cg_goal = x_cg_goal
        self.panel_size = 0.01
        self.t_skin_airfoil = t_skin_airfoil
        # self.plot_mesh_plotly()
        self.x_cg_skin()

    def AC_layout(self):
        b = self.airframe.b
        c_r = self.airframe.c_r
        c_t = self.airframe.c_t
        sweep = self.airframe.le_sweep

        edge_margin = 0.005
        span = b / 2

        # -------------------------
        # Panel grid (centers)
        # -------------------------
        y_edges = np.arange(-span + edge_margin,
                             span - edge_margin + self.panel_size,
                             self.panel_size)

        x_edges = np.arange(edge_margin,
                            c_r + edge_margin,
                            self.panel_size)

        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        x_centers = (x_edges[:-1] + x_edges[1:]) / 2

        Xc, Yc = np.meshgrid(x_centers, y_centers)
        Xc = Xc + np.abs(Yc) * np.tan(sweep)  # Shift x-coordinates based on sweep

        # -------------------------
        # Wing geometry (FIXED SWEEP)
        # -------------------------
        chord = c_r - (np.abs(Yc) / span) * (c_r - c_t)
        chord = np.maximum(chord, 1e-6)

        x_le = np.abs(Yc) * np.tan(sweep)  

        # -------------------------
        # Airfoil thickness
        # -------------------------
        x_over_c = (Xc - x_le) / chord
        x_over_c = np.clip(x_over_c, 0.0, 1.0)

        t_over_c = self.airfoil.local_thickness(x_over_c)

        Z = t_over_c * chord  # meters

        print("Xc shape: ", Xc.shape)
        print("Yc shape: ", Yc.shape)
        print("Z shape: ", Z.shape)
        return Xc, Yc, Z

    def x_cg_skin(self):
        area_out, perimeter_out = sa.airfoil_properties(self.airfoil, chord_length=self.airframe.c_r)
        area_skin = perimeter_out * self.t_skin_airfoil  # Approximate skin area as perimeter times thickness, which is reasonable for thin skins around complex shapes like airfoils
        x_centroid, z_centroid = Calculate_Airfoil_Centroid(self.airfoil, self.airframe.c_r)
        print(f"Z_centroid of skin: {z_centroid:.4f} m from leading edge, x_centroid: {x_centroid:.4f} m from leading edge")



    def plot_mesh_plotly(self):

        X, Y, Z = self.AC_layout()

        b = self.airframe.b
        c_r = self.airframe.c_r
        c_t = self.airframe.c_t
        sweep = self.airframe.le_sweep

        span = b / 2

        # -------------------------
        # Wing outline (2D reference)
        # -------------------------
        y = np.linspace(-b/2, b/2, 300)

        chord = c_r - (np.abs(y) / span) * (c_r - c_t)
        x_le = np.abs(y) * np.tan(sweep)   # FIXED
        x_te = x_le + chord
        print(X)
        # -------------------------
        # Create figure
        # -------------------------
        fig = go.Figure()

        # -------------------------
        # 1. Wing surface
        # -------------------------
        fig.add_trace(go.Surface(
            x=X,
            y=Y,
            z=Z,
            colorscale='Viridis',
            showscale=True,
            name="Wing surface",
            opacity=1.0
        ))

        # -------------------------
        # 2. Leading edge
        # -------------------------
        fig.add_trace(go.Scatter3d(
            x=x_le,
            y=y,
            z=np.zeros_like(y),
            mode='lines',
            line=dict(color='black', width=5),
            name="Leading edge"
        ))

        # -------------------------
        # 3. Trailing edge
        # -------------------------
        fig.add_trace(go.Scatter3d(
            x=x_te,
            y=y,
            z=np.zeros_like(y),
            mode='lines',
            line=dict(color='black', width=5),
            name="Trailing edge"
        ))

        # -------------------------
        # 4. Chordwise guide lines (IMPORTANT FOR SWEEP VISIBILITY)
        # -------------------------
        for yi in np.linspace(-b/2, b/2, 10):

            chord_i = c_r - (np.abs(yi) / span) * (c_r - c_t)
            x_le_i = np.abs(yi) * np.tan(sweep)

            fig.add_trace(go.Scatter3d(
                x=[x_le_i, x_le_i + chord_i],
                y=[yi, yi],
                z=[0, 0],
                mode='lines',
                line=dict(color='gray', width=2),
                showlegend=False
            ))

        # -------------------------
        # Layout
        # -------------------------
        fig.update_layout(
            title=f"HAPS Wing Panel Mesh ({self.panel_size} m × {self.panel_size} m) with Sweep + Thickness",
            scene=dict(
                xaxis_title="x (m)",
                yaxis_title="y (m)",
                zaxis_title="thickness (m)",
                aspectmode="data"
            ),
            margin=dict(l=0, r=0, t=40, b=0)
        )

        fig.show()


if __name__ == "__main__":
    airframe = Airframe.airframe(qc_sweep=np.radians(15), init_polar=False, display=False)
    cg_calculator = CG_comp(x_cg_goal=2.4, airframe=airframe)