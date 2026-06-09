import numpy as np

MH91_COORDS = np.array([
    # upper surface (LE→TE)
    [0.0000,  0.0000],
    [0.0125,  0.0236],
    [0.0250,  0.0334],
    [0.0500,  0.0468],
    [0.0750,  0.0570],
    [0.1000,  0.0651],
    [0.1500,  0.0780],
    [0.2000,  0.0880],
    [0.2500,  0.0956],
    [0.3000,  0.1007],
    [0.3500,  0.1033],
    [0.4000,  0.1034],
    [0.4500,  0.1011],
    [0.5000,  0.0968],
    [0.5500,  0.0908],
    [0.6000,  0.0833],
    [0.6500,  0.0746],
    [0.7000,  0.0648],
    [0.7500,  0.0541],
    [0.8000,  0.0428],
    [0.8500,  0.0311],
    [0.9000,  0.0192],
    [0.9500,  0.0075],
    [1.0000,  0.0000],
    # lower surface (TE→LE)
    [0.9500, -0.0107],
    [0.9000, -0.0192],
    [0.8500, -0.0249],
    [0.8000, -0.0281],
    [0.7500, -0.0292],
    [0.7000, -0.0286],
    [0.6500, -0.0267],
    [0.6000, -0.0237],
    [0.5500, -0.0200],
    [0.5000, -0.0160],
    [0.4500, -0.0119],
    [0.4000, -0.0080],
    [0.3500, -0.0044],
    [0.3000, -0.0012],
    [0.2500,  0.0015],
    [0.2000,  0.0035],
    [0.1500,  0.0048],
    [0.1000,  0.0050],
    [0.0750,  0.0043],
    [0.0500,  0.0027],
    [0.0250,  0.0002],
    [0.0125, -0.0022],
    [0.0000,  0.0000],
])


class Mass_moments:

    def __init__(self):
        self.spar_major_axis = 0.5
        self.spar_minor_axis = 0.4

        self.cg = (-2.19, 0.0, 0.0)
        self.root = (0.0, 0.0, 0.0) # body axis central point

        self.b        = 28.80
        self.half_b = self.b / 2
        self.c        = 1.44
        self.sweep    = np.radians(-15)
        self.dihedral = np.radians(-2)
        self.twist    = np.radians(4.675)
        self.twist_rate = np.radians(self.twist/self.half_b)

        # MH91 leading-20% battery box geometry
        self.battery_chord_fraction     = 0.20
        self.battery_thickness_fraction = 0.10
        self.battery_density            = 2180.0

    # ------------------------------------------------------------------ #
    #  internal helpers
    # ------------------------------------------------------------------ #

    def _spar_direction(self):
        """Unit vector along the spar in FD body axes (x fwd, y right, z down)."""
        Λ = self.sweep
        Γ = self.dihedral
        return np.array([
            np.cos(Γ) * np.sin(Λ),
            np.cos(Γ) * np.cos(Λ),
            -np.sin(Γ),
        ])

    def _parallel_axis(self, mass, d):
        """Parallel-axis inertia contribution for a point mass at offset d."""
        dx, dy, dz = d
        return mass * np.array([
            [dy**2 + dz**2, -dx*dy,        -dx*dz       ],
            [-dx*dy,         dx**2 + dz**2, -dy*dz       ],
            [-dx*dz,        -dy*dz,          dx**2 + dy**2],
        ])

    def _mh91_perimeter(self):
        """Arc length of full MH91 airfoil perimeter [m]."""
        xy = MH91_COORDS * self.c
        diffs = np.diff(xy, axis=0)
        return np.hypot(diffs[:, 0], diffs[:, 1]).sum()

    def _mh91_skin_centroid(self):
        """
        Arc-length-weighted centroid of the MH91 skin in body axes.
        Returns (x_offset, z_offset) relative to the leading edge.
        x: chordwise (positive aft), z: positive down (negated from airfoil y).
        """
        xy  = MH91_COORDS * self.c
        dls = np.hypot(*np.diff(xy, axis=0).T)
        mids = 0.5 * (xy[:-1] + xy[1:])
        x_c =  (mids[:, 0] * dls).sum() / dls.sum()
        z_c = -(mids[:, 1] * dls).sum() / dls.sum()   # flip sign: airfoil y-up → FD z-down
        return x_c, z_c

    # ------------------------------------------------------------------ #
    #  spar
    # ------------------------------------------------------------------ #

    def spar_inertia_fd(
            self,
            mass,
            length,
    ):
        """
        Inertia tensor of a swept/dihedral spar (slender rod) about the CG.

        Parameters
        ----------
        mass       : float   spar mass [kg]
        length     : float   spar length [m]
        sweep_deg  : float   sweep angle [deg]
        dihedral_deg: float  dihedral angle [deg]
        twist_deg  : float   (stored for completeness, not used)
        cg         : tuple   aircraft CG in body axes
        root       : tuple   spar root in body axes
        """
        u  = self._spar_direction()
        rc = np.array(self.root) + 0.5 * length * u
        d  = rc - np.array(self.cg)

        Ic  = (mass * length**2 / 12.0) * (np.eye(3) - np.outer(u, u))
        return Ic + self._parallel_axis(mass, d)

    # ------------------------------------------------------------------ #
    #  batteries
    # ------------------------------------------------------------------ #

    def batteries_inertia_fd(
            self,
            total_mass,
            span_sections,
    ):
        """
        Inertia tensor of batteries in the leading 20% of the MH91 airfoil.

        Total mass + density back-calculates exact volume; that volume is
        distributed across span_sections proportionally to section length.

        Parameters
        ----------
        total_mass    : float              known battery mass [kg]
        span_sections : list of (float, float)
                        (y_inner, y_outer) pairs [m] along the half-span
        sweep_deg     : float
        dihedral_deg  : float
        chord         : float, optional    defaults to self.c
        cg            : tuple
        root          : tuple

        Returns
        -------
        I_total      : np.ndarray (3×3)  [kg·m²]
        section_info : list of dict
        """

        u = self._spar_direction()

        batt_width  = self.battery_chord_fraction     * self.c
        batt_height = self.battery_thickness_fraction * self.c
        A_cross     = batt_width * batt_height

        total_volume = total_mass / self.battery_density

        section_lengths = np.array([y2 - y1 for y1, y2 in span_sections])
        if np.any(section_lengths <= 0):
            raise ValueError("All span sections must have y_outer > y_inner.")

        volume_per_section = total_volume * (section_lengths / section_lengths.sum())
        fill_lengths       = volume_per_section / A_cross

        for i, (fl, sl) in enumerate(zip(fill_lengths, section_lengths)):
            if fl > sl + 1e-9:
                raise ValueError(
                    f"Section {i} ({span_sections[i]}): fill length {fl:.3f} m "
                    f"exceeds available span {sl:.3f} m. Add more sections."
                )

        I_total      = np.zeros((3, 3))
        section_info = []

        for (y_inner, _), fill_length in zip(span_sections, fill_lengths):
            mass    = self.battery_density * A_cross * fill_length
            p_inner = np.array(self.root) + y_inner * u
            rc      = p_inner + 0.5 * fill_length * u
            d       = rc - np.array(self.cg)

            L2  = fill_length**2
            wh2 = batt_width**2 + batt_height**2
            Ic  = (mass / 12.0) * (L2 * np.eye(3) + (wh2 - L2) * np.outer(u, u))

            I_total += Ic + self._parallel_axis(mass, d)
            section_info.append({
                "section":     (y_inner, y_inner + fill_length),
                "fill_length": fill_length,
                "mass":        mass,
                "volume":      A_cross * fill_length,
            })

        return I_total, section_info

    # ------------------------------------------------------------------ #
    #  skin
    # ------------------------------------------------------------------ #

    def skin_inertia_fd(
            self,
            skin_density,
            skin_thickness,
            span_sections,
    ):
        """
        Inertia tensor of a constant-thickness shell skin following the MH91
        profile, over user-specified spanwise sections.

        Parameters
        ----------
        skin_density    : float   material density [kg/m³]
        skin_thickness  : float   shell thickness [m]
        span_sections   : list of (float, float)
        sweep_deg       : float
        dihedral_deg    : float
        chord           : float, optional
        cg              : tuple
        root            : tuple

        Returns
        -------
        I_total      : np.ndarray (3×3)  [kg·m²]
        section_info : list of dict
        """

        u         = self._spar_direction()
        perimeter = self._mh91_perimeter()
        x_skin, z_skin = self._mh91_skin_centroid()

        # Gyration radius of the skin cross-section about the spar axis
        # (thin-ring equivalent: r² = P² / 4π²)
        r_gyr2 = perimeter**2 / (4 * np.pi**2)

        I_total      = np.zeros((3, 3))
        section_info = []

        for (y_inner, y_outer) in span_sections:
            seg_len = y_outer - y_inner
            if seg_len <= 0:
                raise ValueError(
                    f"Invalid section ({y_inner}, {y_outer}): y_outer must be > y_inner."
                )

            mass = skin_density * skin_thickness * perimeter * seg_len

            # Centroid: spar midpoint + cross-section offset
            p_inner = np.array(self.root) + y_inner * u
            rc_spar = p_inner + 0.5 * seg_len * u
            rc      = rc_spar + np.array([x_skin, 0.0, z_skin])
            d       = rc - np.array(self.cg)

            L2 = seg_len**2
            Ic = mass * (
                (L2 / 12.0) * (np.eye(3) - np.outer(u, u))
                + r_gyr2    *  np.outer(u, u)
            )

            I_total += Ic + self._parallel_axis(mass, d)
            section_info.append({
                "section":   (y_inner, y_outer),
                "span_len":  seg_len,
                "perimeter": perimeter,
                "mass":      mass,
            })

        return I_total, section_info

    def rib_point_mass_inertia_full(
            self,
            mass,
            y_span,
    ):
        """
        Rib point-mass inertia about aircraft CG
        with sweep, dihedral, and twist included.

        Flight dynamics axes:
            x forward
            y right wing
            z down
        """

        Λ = self.sweep
        Γ = self.dihedral

        θ0 = self.twist
        θt = self.twist_rate

        x0, y0, z0 = self.root
        xcg, ycg, zcg = self.cg

        # spanwise twist (for completeness only)
        twist = θ0 + θt * y_span

        # position in aircraft body axes (root frame)
        x_rib = x0 + y_span * np.tan(Λ)
        y_rib = y0 + y_span
        z_rib = z0 - y_span * np.tan(Γ)

        # shift to CG frame
        x = x_rib - xcg
        y = y_rib - ycg
        z = z_rib - zcg

        # point-mass inertia about CG
        I = mass * np.array([
            [y ** 2 + z ** 2, -x * y, -x * z],
            [-x * y, x ** 2 + z ** 2, -y * z],
            [-x * z, -y * z, x ** 2 + y ** 2]
        ])

        return I, twist

    def wing_rib_inertia_full(
            self,
            ribs,
    ):

        I_total = np.zeros((3, 3))
        M_total = 0.0

        for mass, y in ribs:
            I, _ = self.rib_point_mass_inertia_full(
                mass=mass,
                y_span=y,
            )

            I_total += I
            M_total += mass

        return M_total, I_total


# ======================================================================== #
#  Example
# ======================================================================== #

cs        = Mass_moments()
half_span = cs.b / 2   # 14.4 m

# --- spar (one half-wing) ---
I_spar = cs.spar_inertia_fd(
    mass=20,
    length=half_span,
)
print("Spar inertia tensor [kg·m²]:")
print(np.round(I_spar, 4))

# --- batteries: inner and outer third of half-span ---
battery_sections = [
    (0.0,               half_span / 3),   # 0.0 – 4.8 m
    (2 * half_span / 3, half_span),        # 9.6 – 14.4 m
]

I_batt, batt_info = cs.batteries_inertia_fd(
    total_mass=45.0,
    span_sections=battery_sections,
)
print("\nBattery sections:")
for s in batt_info:
    print(f"  {s['section'][0]:.2f}–{s['section'][1]:.2f} m | "
          f"fill {s['fill_length']:.3f} m | mass {s['mass']:.3f} kg")
print(f"Battery inertia tensor [kg·m²]:")
print(np.round(I_batt, 4))

# --- skin: full half-span ---
I_skin, skin_info = cs.skin_inertia_fd(
    skin_density=1390.0,
    skin_thickness=0.002,
    span_sections=[(0.0, half_span)],
)
print("\nSkin sections:")
for s in skin_info:
    print(f"  {s['section'][0]:.2f}–{s['section'][1]:.2f} m | "
          f"perimeter {s['perimeter']:.4f} m | mass {s['mass']:.3f} kg")
print(f"Skin inertia tensor [kg·m²]:")
print(np.round(I_skin, 4))

# --- ribs: full half span ---
ribs = [
    (0.6, 0),
    (0.6, 0.5*cs.half_b),
    (0.6, cs.half_b)
]

M_ribs, I_ribs = cs.wing_rib_inertia_full(ribs)
print("\nWing ribs mass:")
print(np.round(M_ribs, 4))
print(f"Rib inertia tensor [kg·m²]:")
print(np.round(I_ribs, 4))

# --- combined ---
I_total = I_spar + I_batt + I_skin + I_ribs
print("\nCombined inertia tensor [kg·m²]:")
print(np.round(I_total, 4))