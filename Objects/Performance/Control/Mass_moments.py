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
        # self.spar_major_axis = 0.5
        # self.spar_minor_axis = 0.4

        self.cg = (-2.19, 0.0, 0.0)
        self.root = (-0.7, 0.0, 0.0)
        # self.body_centre = (0.0, 0.0, 0.0) # body axis central point

        self.b        = 28.80
        self.half_b   = self.b / 2
        self.c        = 1.44
        self.S        = self.b * self.c
        self.AR       = 20
        self.sweep    = np.radians(15)
        self.dihedral = np.radians(2)
        self.twist    = np.radians(4.675) # twist not used - assumption that effect is minimal
        self.twist_rate = self.twist/self.half_b

        # MH91 leading-20% battery box geometry
        self.battery_chord_fraction     = 0.20
        self.battery_thickness_fraction = 0.10
        self.battery_density            = 2180.0

        self.m_total = 170.0  # MTOM [kg]

        I_spar = self.spar_inertia_fd(mass=20, length=self.half_b)

        I_batt, _ = self.batteries_inertia_fd(
            total_mass=40.5,
            span_sections=[(0.0, self.half_b / 3), (2 * self.half_b / 3, self.half_b)],
        )
        I_skin, _ = self.skin_inertia_fd(
            skin_density=1390.0,
            skin_thickness=0.0002,
            span_sections=[(0.0, self.half_b)],
        )
        ribs = [(0.6, 0), (0.6, 0.5 * self.half_b), (0.6, self.half_b)]
        _, I_ribs = self.wing_rib_inertia_full(ribs)

        I_solar, _ = self.solar_panel_inertia_fd(
            surface_density=0.665,
            span_limits=(self.half_b * 0.05, self.half_b * 0.95),
            chord_limits=(0.05, 0.95),
        )
        motors = [(2.0, 0.25 * self.half_b, 0.5), (2.0, 0.5 * self.half_b, 0.5)]
        _, I_motors = self.motor_inertia_full(motors)

        I_payload = self.payload_inertia_fd(mass=20.0, x_pos=0.75)

        I_combined = I_spar + I_batt + I_skin + I_ribs + I_solar + I_motors + I_payload

        # KX2, KZ2, KXZ — normalised by b
        (self.k_x_nd,
         self.k_y_nd,
         self.k_z_nd,
         self.k_xz_nd) = self.non_dimensional_radius_of_gyration(I_combined, m_total, self.b)


    # ------------------------------------------------------------------ #
    #  internal helpers
    # ------------------------------------------------------------------ #

    def _spar_direction(self):
        """Unit vector along the spar in FD body axes (x fwd, y right, z down)."""
        Λ = self.sweep
        Γ = self.dihedral
        return np.array([
            -np.cos(Γ) * np.sin(Λ),
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

    def _include_second_half(self, I_left):
        T = np.diag([1, -1, 1])

        I_right = T @ I_left @ T.T

        I_total = I_left + I_right
        return I_total

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

        I_left = Ic + self._parallel_axis(mass, d)

        I_total = self._include_second_half(I_left)

        return I_total

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

        I_left      = np.zeros((3, 3))
        section_info = []

        for (y_inner, _), fill_length in zip(span_sections, fill_lengths):
            mass    = self.battery_density * A_cross * fill_length
            p_inner = np.array(self.root) + y_inner * u # todo check if root use is correct and not body centre
            rc      = p_inner + 0.5 * fill_length * u
            d       = rc - np.array(self.cg)

            L2  = fill_length**2
            wh2 = batt_width**2 + batt_height**2
            Ic  = (mass / 12.0) * (L2 * np.eye(3) + (wh2 - L2) * np.outer(u, u))

            I_left += Ic + self._parallel_axis(mass, d)

            section_info.append({
                "section":     (y_inner, y_inner + fill_length),
                "fill_length": fill_length,
                "mass":        mass,
                "volume":      A_cross * fill_length,
            })

        I_total = self._include_second_half(I_left)

        return I_total, section_info # doubled due to wing symmetry

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

        I_left      = np.zeros((3, 3))
        section_info = []

        for (y_inner, y_outer) in span_sections:
            seg_len = y_outer - y_inner
            if seg_len <= 0:
                raise ValueError(
                    f"Invalid section ({y_inner}, {y_outer}): y_outer must be > y_inner."
                )

            mass = skin_density * skin_thickness * perimeter * seg_len

            # Centroid: spar midpoint + cross-section offset
            p_inner = np.array(self.root) + y_inner * u # todo
            rc_spar = p_inner + 0.5 * seg_len * u
            rc      = rc_spar + np.array([-x_skin, 0.0, z_skin])
            d       = rc - np.array(self.cg)

            L2 = seg_len**2
            Ic = mass * (
                (L2 / 12.0) * (np.eye(3) - np.outer(u, u))
                + r_gyr2    *  np.outer(u, u)
            )

            I_left += Ic + self._parallel_axis(mass, d)
            section_info.append({
                "section":   (y_inner, y_outer),
                "span_len":  seg_len,
                "perimeter": perimeter,
                "mass":      mass,
            })

            T = np.diag([1, -1, 1])

            I_right = T @ I_left @ T.T

            I_total = I_left + I_right

        return I_total, section_info # doubled due to wing symmetry

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
        x_rib = x0 - y_span * np.tan(Λ)
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

        I_left = np.zeros((3, 3))
        M_total = 0.0

        for mass, y in ribs:
            I, _ = self.rib_point_mass_inertia_full(
                mass=mass,
                y_span=y,
            )

            I_left += I
            M_total += mass

        T = np.diag([1, -1, 1])

        I_right = T @ I_left @ T.T

        I_total = I_left + I_right

        return M_total, I_total

    def solar_panel_inertia_fd(
            self,
            surface_density,  # kg/m²
            span_limits,  # (y_start, y_end) in meters
            chord_limits,  # (x_start_frac, x_end_frac) as fraction of chord,
    ):
        """
        Inertia tensor of a solar panel modeled as a thin plate with surface mass density.

        Parameters
        ----------
        surface_density : float   surface mass density of the solar panel [kg/m²]
        span_limits     : tuple   (y_start, y_end) spanwise limits of the solar panel [m]
        chord_limits    : tuple   (x_start_frac, x_end_frac) chordwise limits (fraction of chord)
        sweep_deg       : float   sweep angle of the wing [deg]
        dihedral_deg    : float   dihedral angle of the wing [deg]
        cg              : tuple   aircraft CG in body axes
        root            : tuple   wing root in body axes
        chord           : float   chord length of the wing [m] (defaults to self.c)

        Returns
        -------
        I_total         : np.ndarray (3×3)  inertia tensor about the aircraft CG [kg·m²]
        section_info    : dict                 details about the solar panel
        """

        y_start, y_end = span_limits
        x_start_frac, x_end_frac = chord_limits

        # Unit vector along the spar (spanwise direction)
        u_spar = self._spar_direction()

        # Chordwise limits in meters
        x_start = x_start_frac * self.c
        x_end = x_end_frac * self.c

        # Area and mass of the solar panel
        area = (y_end - y_start) * (x_end - x_start)
        mass = surface_density * area

        # Centroid of the solar panel in the local wing coordinate system
        x_centroid_local = 0.5 * (x_start + x_end)
        y_centroid_local = 0.5 * (y_start + y_end)
        z_centroid_local = 0.0  # Assuming the solar panel lies on the wing surface

        # Centroid in the global body axes
        centroid_global = (
                np.array(self.root) +
                y_centroid_local * u_spar +
                np.array([-x_centroid_local, 0.0, 0.0])  # Assuming chord is aligned with x-axis
        )

        # Offset from the aircraft CG
        d = centroid_global - np.array(self.cg)

        # Inertia tensor of a thin plate about its centroid in the local coordinate system
        L = y_end - y_start  # Spanwise length
        W = x_end - x_start  # Chordwise length
        I_local = (mass / 12) * np.array([
            [W ** 2, 0, 0],
            [0, L ** 2 + W ** 2, 0],
            [0, 0, L ** 2]
        ])

        # Rotation matrix from local to global body axes
        e_span = u_spar

        e_chord = np.array([1.0, 0.0, 0.0])
        e_chord -= np.dot(e_chord, e_span) * e_span
        e_chord /= np.linalg.norm(e_chord)

        e_normal = np.cross(e_span, e_chord)

        R = np.column_stack((e_chord, e_span, e_normal))

        # Rotate the local inertia tensor to the global body axes
        I_rotated = R @ I_local @ R.T

        # Apply the parallel axis theorem
        I_left = I_rotated + self._parallel_axis(mass, d)

        T = np.diag([1, -1, 1])

        I_right = T @ I_left @ T.T

        I_total = I_left + I_right

        # Section info
        section_info = {
            "mass": mass,
            "area": area,
            "centroid_global": centroid_global,
            "span_limits": span_limits,
            "chord_limits": chord_limits,
        }

        return I_total, section_info

    def payload_inertia_fd(
            self,
            mass,
            x_pos,
    ):
        x0, y0, z0 = self.root
        xcg, ycg, zcg = self.cg

        x_payload = x0 + x_pos
        y_payload = y0
        z_payload = z0

        x = x_payload - xcg
        y = y_payload - ycg
        z = z_payload - zcg

        I = mass * np.array([
            [y ** 2 + z ** 2, -x * y, -x * z],
            [-x * y, x ** 2 + z ** 2, -y * z],
            [-x * z, -y * z, x ** 2 + y ** 2]
        ])

        return I

    def motor_inertia_fd(
            self,
            mass,
            y_span,
            x_offset=0.5,
    ):

        """
        Motor point-mass inertia about aircraft CG
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

        # position in aircraft body axes (root frame)
        x_motor = x0 - y_span * np.tan(Λ) + x_offset
        y_motor = y0 + y_span
        z_motor = z0 - y_span * np.tan(Γ)

        # shift to CG frame
        x = x_motor - xcg
        # print("Motor x", x)
        y = y_motor - ycg
        z = z_motor - zcg

        # point-mass inertia about CG
        I = mass * np.array([
            [y ** 2 + z ** 2, -x * y, -x * z],
            [-x * y, x ** 2 + z ** 2, -y * z],
            [-x * z, -y * z, x ** 2 + y ** 2]
        ])

        return I

    def motor_inertia_full(
            self,
            motors,
    ):

        I_left = np.zeros((3, 3))
        M_total = 0.0

        for mass, y, x_offset in motors:
            I = self.motor_inertia_fd(
                mass=mass,
                y_span=y,
                x_offset=x_offset,
            )

            I_left += I
            M_total += mass

            T = np.diag([1, -1, 1])

            I_right = T @ I_left @ T.T

            I_total = I_left + I_right

        return M_total, I_total # doubled due to wing symmetry

    def radius_of_gyration(self, I_total, M_total):
        """
        Calculate the dimensional radius of gyration about the x, y, and z axes.

        Parameters
        ----------
        I_total : np.ndarray (3x3)
            Combined inertia tensor of the aircraft [kg·m²].
        M_total : float
            Total mass of the aircraft [kg].

        Returns
        -------
        k_x, k_y, k_z : float
            Radius of gyration about the x, y, and z axes [m].
        """
        I_xx = I_total[0, 0]
        I_yy = I_total[1, 1]
        I_zz = I_total[2, 2]
        I_xz = I_total[0,2]

        k_x = np.sqrt(I_xx / M_total)
        k_y = np.sqrt(I_yy / M_total)
        k_z = np.sqrt(I_zz / M_total)
        k_xz = -I_xz / M_total

        return k_x, k_y, k_z, k_xz

    def non_dimensional_radius_of_gyration(self, I_total, M_total, reference_length=None):
        if reference_length is None:
            reference_length = self.b

        k_x, k_y, k_z, k_xz = self.radius_of_gyration(I_total, M_total)

        k_x_nd = k_x / reference_length
        k_y_nd = k_y / reference_length  # ← use reference_length consistently, not hardcoded c
        k_z_nd = k_z / reference_length
        k_xz_nd = k_xz / (reference_length ** 2)

        return k_x_nd, k_y_nd, k_z_nd, k_xz_nd

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
    (0.0, half_span / 3),   # 0.0 – 4.8 m
    (2 * half_span / 3, half_span),        # 9.6 – 14.4 m
]

I_batt, batt_info = cs.batteries_inertia_fd(
    total_mass=40.5,
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
    skin_thickness=0.0002,
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

# Example: Solar panel covering the middle 50% of the half-span and 80% of the chord
solar_span_limits = (half_span * 0.05, half_span*0.95)  # 3.6 m to 10.8 m
solar_chord_limits = (0.05, 0.95)  # 10% to 90% of the chord

I_solar, solar_info = cs.solar_panel_inertia_fd(
    surface_density=0.665,  # kg/m²
    span_limits=solar_span_limits,
    chord_limits=solar_chord_limits,
)

print("\nSolar panel inertia tensor [kg·m²]:")
print(np.round(I_solar, 4))
print("\nSolar panel info:")
print(f"Mass: {solar_info['mass']:.3f} kg")
print(f"Area: {solar_info['area']:.3f} m²")
print(f"Centroid (global): {solar_info['centroid_global']}")

motors = [
    (2.0, 0.25*cs.half_b, 0.5),
    (2.0, 0.5*cs.half_b, 0.5)
]

M_motors, I_motors = cs.motor_inertia_full(motors)
print("\nmotors mass:")
print(np.round(M_motors, 4))
print(f"motor inertia tensor [kg·m²]:")
print(np.round(I_motors, 4))

# Calculate payload inertia tensor
I_payload = cs.payload_inertia_fd(mass=20.0, x_pos=0.75)
print("\nPayload inertia tensor [kg·m²]:")
print(np.round(I_payload, 4))
# --- combined ---

I_total = I_spar + I_batt + I_skin +I_solar + I_ribs + I_motors + I_payload
print("\nCombined inertia tensor [kg·m²]:")
print(np.round(I_total, 4))

# --- Calculate dimensional radius of gyration ---
k_x, k_y, k_z, k_xz = cs.radius_of_gyration(I_total, M_total=200.0)
print("\nDimensional Radius of Gyration:")
print(f"k_x (roll axis): {k_x:.4f} m")
print(f"k_y (pitch axis): {k_y:.4f} m")
print(f"k_z (yaw axis): {k_z:.4f} m")
print(f"k_xz: {k_xz:.4f} m")

# --- Calculate non-dimensional radius of gyration (normalized by wingspan) ---
k_x_nd, k_y_nd, k_z_nd, k_xz_nd = cs.non_dimensional_radius_of_gyration(I_total, M_total=200.0)
print("\nNon-Dimensional Radius of Gyration (normalized by wingspan):")
print(f"k_x/b: {k_x_nd:.4f}")
print(f"k_y/b: {k_y_nd:.4f}")
print(f"k_z/b: {k_z_nd:.4f}")
print(f"k_xz/b^2: {k_xz_nd:.4f} m")