import numpy as np
import aerosandbox as asb
def distributed_engine_weight(array_half_span, y_positions, mass_eng, width_eng, g=9.81, shape="rect"):
    """
    Returns a 1D array w(y) [N/m] defined on array_half_span (meters from root to tip)
    representing distributed engine weight along the half-span.

    - array_half_span: 1D monotonic increasing grid (m)
    - y_positions: engine center locations on half-span (m)
    - mass_eng: mass per engine (kg)
    - width_eng: distribution width along span (m)
    - shape: "rect" (uniform over width) or "tri" (triangular over width)
    """
    y = np.asarray(array_half_span)
    dy = np.diff(y)
    # local grid spacing per node for integrating (node-centered control volumes)
    dy_node = np.empty_like(y)
    dy_node[1:-1] = 0.5 * (dy[:-1] + dy[1:])
    dy_node[0] = dy[0]
    dy_node[-1] = dy[-1]

    w = np.zeros_like(y, dtype=float)  # N/m

    for yc in np.atleast_1d(y_positions):
        if shape == "rect":
            # uniform load over [yc - width/2, yc + width/2]
            mask = (y >= yc - width_eng / 2) & (y <= yc + width_eng / 2)
            if not np.any(mask):
                # if the width is smaller than grid resolution, fall back to nearest node
                j = np.argmin(np.abs(y - yc))
                w[j] += (mass_eng * g) / dy_node[j]
            else:
                # distribute so that integral equals mass_eng*g
                w[mask] += (mass_eng * g) / (width_eng)

        elif shape == "tri":
            # triangular distribution: peak at center, zero at edges, normalized to mass_eng*g
            # base width = width_eng, so area = 0.5 * base * peak => peak = 2*W/base
            W = mass_eng * g
            peak = 2 * W / width_eng  # N/m
            dist = np.abs(y - yc)
            mask = dist <= width_eng / 2
            w[mask] += peak * (1 - (2 * dist[mask] / width_eng))
        else:
            raise ValueError("shape must be 'rect' or 'tri'")

    return w

class VM:
    def __init__(self, span, dy, lift_distribution=None, weight_wing_distribution=None, weight_battery_distribution=None, weight_propulsion_distribution=None):
        self.half_span = span / 2
        self.dy = dy
        self.y_int = np.arange(0, self.half_span + dy, dy)  # y positions along the half-span for integration
        self.LiftDistribution = lift_distribution
        self.Weight_Wing_distribution = weight_wing_distribution
        self.Weight_Battery_distribution = weight_battery_distribution
        self.Weight_Propulsion_distribution = weight_propulsion_distribution

        compute_forces = self.compute_force_distributions()
        print("Force distributions computed successfully.")
    #Convention: normal force out is positive, normal force in is negative. 
    #Moment wing up is positive, wing down is negative. 
    #Shear force up is positive, shear force down is negative.

    def compute_force_distributions(self):
        q = (
            np.asarray(self.LiftDistribution)
            - np.asarray(self.Weight_Wing_distribution)
            - np.asarray(self.Weight_Battery_distribution)
            - np.asarray(self.Weight_Propulsion_distribution)
        )

        if len(q) != len(self.y_int):
            raise ValueError("Force distributions must have the same length as y_int.")

        y = self.y_int

        # Work in reversed order so that boundary conditions at the tip are enforced:
        # V(tip)=0, M(tip)=0
        y_r = y[::-1]
        q_r = q[::-1]

        dy_r = np.diff(y_r)  # negative values (since y_r decreases)
        # Trapezoidal integration
        dV_r = -0.5 * (q_r[1:] + q_r[:-1]) * dy_r
        V_r = np.concatenate(([0.0], np.cumsum(dV_r)))

        dM_r = 0.5 * (V_r[1:] + V_r[:-1]) * dy_r
        M_r = np.concatenate(([0.0], np.cumsum(dM_r)))

        # Flip back to increasing y for plotting
        self.total_normal_force_distribution = q
        self.V = V_r[::-1]
        self.M = M_r[::-1]

        return self.total_normal_force_distribution, self.V, self.M

    def Plot(self):
        import matplotlib.pyplot as plt

        plt.figure(figsize=(12, 8))
        plt.subplot(3, 1, 1)
        plt.plot(self.y_int, self.total_normal_force_distribution, label="Total Force Distribution")
        plt.xlabel("Spanwise Position (m)")
        plt.ylabel("Total Force (N/m)")
        plt.title("Total Force Distribution Along the Span")
        plt.grid()
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(self.y_int, self.V, label="Shear Force Distribution", color='orange')
        plt.xlabel("Spanwise Position (m)")
        plt.ylabel("Shear Force (N)")
        plt.title("Shear Force Distribution Along the Span")
        plt.grid()
        plt.legend()

        plt.subplot(3, 1, 3)
        plt.plot(self.y_int, self.M, label="Bending Moment Distribution", color='green')
        plt.xlabel("Spanwise Position (m)")
        plt.ylabel("Bending Moment (Nm)")
        plt.title("Bending Moment Distribution Along the Span")
        plt.grid()
        plt.legend()

        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    span = 30.1
    dy = 0.1
    array_half_span = np.arange(0, span/2 + dy, dy)
    #Reasonable AHAPS
    span = 30.1
    half_span = span / 2
    y = array_half_span
    g = 9.81

    # Aircraft mass + required lift (1g)
    m_total = 200.0
    W_total = m_total * g
    L_total_required = W_total

    # Elliptical lift distribution L'(y) [N/m] on the half span
    b_semi = half_span
    phi = np.sqrt(np.clip(1 - (y / b_semi) ** 2, 0, None))

    # Scale ellipse so that ∫_0^{b_semi} L'(y) dy = L_total_required / 2
    # since ∫_0^{b_semi} sqrt(1-(y/b_semi)^2) dy = b_semi * pi/4
    L0 = (L_total_required / 2) / (b_semi * (np.pi / 4))
    lift_distribution = L0 * phi

    # Example mass splits (adjust if you know better)
    m_wing_total = 70.0        # kg for both wings
    m_battery_total = 50.0     # kg total batteries (distributed in wings)

    weight_wing_distribution = ((m_wing_total * g) / 2 / b_semi) * np.ones_like(y)
    weight_battery_distribution = ((m_battery_total * g) / 2 / b_semi) * np.ones_like(y)

    # Engines
    n_eng = 4
    n_half = n_eng // 2

    # exactly 2 engines on the half span, not at root/tip
    y_positions = b_semi * (np.arange(1, n_half + 1) / (n_half + 1))

    mass_eng = 5.0             # kg per engine
    width_eng = 0.4            # m footprint along span
    weight_propulsion_distribution = distributed_engine_weight(
        array_half_span=array_half_span,
        y_positions=y_positions,
        mass_eng=mass_eng,
        width_eng=width_eng,
        g=9.81,
        shape="tri",     # or "tri"
    )
    
    VM = VM(span, dy, lift_distribution=lift_distribution, weight_wing_distribution=weight_wing_distribution, weight_battery_distribution=weight_battery_distribution, weight_propulsion_distribution=weight_propulsion_distribution)
    VM.Plot()
    