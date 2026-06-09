from aerosandbox import AeroSandboxObject
import numpy as np

class Mass_moments():
    def __init__(self):
        self.spar_major_axis = 0.5
        self.spar_minor_axis = 0.4

        self.b = 28.80
        self.c = 1.44
        self.sweep = np.radians(- 15)
        self.dihedral = np.radians(- 2)
        self.twist = np.radians(4.675)

    import numpy as np

    def spar_inertia_fd(
            self,
            mass,
            length,
            sweep_deg,
            dihedral_deg,
            twist_deg=4.675,
            cg=(-2.19, 0.0, 0.0),
            root=(0.0, 0.0, 0.0)
    ):
        """
        Centroid + CG-shift inertia tensor of a swept/dihedral spar
        in standard flight dynamics body axes:

            x: forward
            y: right wing
            z: down

        Parameters
        ----------
        mass : float
            Spar mass [kg]
        length : float
            Spar length [m]
        sweep_deg : float
            Sweep angle Λ [deg]
        dihedral_deg : float
            Dihedral angle Γ [deg]
        twist_deg : float
            Twist angle (not used in centroid inertia, included for completeness)
        cg : tuple
            Aircraft CG location in body axes
        root : tuple
            Spar root location in body axes
        """

        Λ = np.radians(sweep_deg)
        Γ = np.radians(dihedral_deg)
        twist = np.radians(twist_deg)  # stored for extension

        # Unit direction vector of spar centerline (FD axes)
        u = np.array([
            np.cos(Γ) * np.sin(Λ),  # x (forward)
            np.cos(Γ) * np.cos(Λ),  # y (right)
            -np.sin(Γ)  # z (down positive)
        ])

        # centroid position
        rc = np.array(root) + 0.5 * length * u

        # shift to reference point (usually CG)
        d = rc - np.array(cg)
        dx, dy, dz = d

        # centroidal inertia tensor (slender rod assumption)
        Ic = (mass * length ** 2 / 12.0) * (np.eye(3) - np.outer(u, u))

        # parallel axis theorem
        Ipa = mass * np.array([
            [dy ** 2 + dz ** 2, -dx * dy, -dx * dz],
            [-dx * dy, dx ** 2 + dz ** 2, -dy * dz],
            [-dx * dz, -dy * dz, dx ** 2 + dy ** 2]
        ])

        return Ic + Ipa


cs = Mass_moments()
print(cs.spar_inertia_fd(mass=20,length=28.8/2,sweep_deg=15,dihedral_deg=2))
