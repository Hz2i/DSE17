import aerosandbox as asb
import aerosandbox.numpy as np


def make_wing(
    span: float,
    ys_over_half_span: np.ndarray,
    chords: np.ndarray,
    twists: np.ndarray,
    offsets: np.ndarray = None,
    heave_displacements: np.ndarray = None,
    twist_displacements: np.ndarray = None,
    x_ref_over_chord: float = 0.33,
    airfoil: asb.Airfoil = asb.Airfoil("dae11"),
    color="black",
) -> asb.Wing:
    """
    Generates a wing based on a given set of per-cross-section characteristics.

    Args:

        span: Span of the wing [meters].

        ys_over_half_span: Array of the y-locations of each cross-section, normalized by half-span. Should be between 0 and 1.

        chords: Array of the chord lengths of each cross-section [meters].

        twists: Array of the twist angles of each cross-section [degrees].

        offsets: Array of the x-offsets of the leading edge of each cross-section [meters]. Defaults to -chords / 4, yielding an unswept quarter-chord.

        heave_displacements: Array of the vertical displacements of the shear center of each cross-section [meters]. Defaults to zero.

        twist_displacements: Array of the twist displacements of each cross-section [degrees], as measured about the shear center. Defaults to zero.

        x_ref_over_chord: The x-location of the shear center (i.e., torsion axis), normalized by the chord. Defaults to 0.33.

        airfoil: The airfoil to use for all cross-sections. Defaults to the DAE11.

    Returns:

        A Wing object.
    """
    if offsets is None:
        offsets = -chords / 4
    if heave_displacements is None:
        heave_displacements = np.zeros_like(ys_over_half_span)
    if twist_displacements is None:
        twist_displacements = np.zeros_like(ys_over_half_span)

    xsecs = []

    for i in range(len(ys_over_half_span)):
        xyz_le = np.array(
            [-chords[i] * x_ref_over_chord, ys_over_half_span[i] * (span / 2), 0]
        )
        xyz_le = (
            np.rotation_matrix_3D(
                angle=np.radians(twists[i] + twist_displacements[i]), axis="y"
            )
            @ xyz_le
        )
        xyz_le += np.array(
            [offsets[i] + chords[i] * x_ref_over_chord, 0, heave_displacements[i]]
        )

        xsecs.append(
            asb.WingXSec(
                xyz_le=xyz_le,
                chord=chords[i],
                twist=twists[i] + twist_displacements[i],
                airfoil=airfoil,
            )
        )

    return asb.Wing(
        symmetric=True,
        xsecs=xsecs,
        color=color,
    )

span = 30.1
ys_over_half_span = np.linspace(0, 1)
taper_ratio = 1
chord =1.2
chords = np.linspace(chord, chord * taper_ratio)
twists = np.linspace(0, 0)
sweep = 15 #deg
offsets = (ys_over_half_span*span/2) * np.tan(np.radians(sweep))  # x pos from LE
print(offsets)
airfoil = asb.Airfoil("e387")

wing = make_wing(
    span=span,
    ys_over_half_span=ys_over_half_span,
    chords=chords,
    twists=twists,
    offsets=offsets,
    airfoil=airfoil,
    color="black",
)

wing_deformed = make_wing(
    span=span,
    ys_over_half_span=ys_over_half_span,
    chords=chords,
    twists=twists,
    offsets=offsets,
    heave_displacements=0.5 * np.cos(np.linspace(0, 2 * np.pi)),
    twist_displacements=25 * np.cos(np.linspace(0, np.pi)),
    color="teal",
)

import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

preset_view_angles = np.array([["XZ", "-YZ"], ["XY", "left_isometric"]], dtype="O")

fig, axs = p.figure3d(
    nrows=preset_view_angles.shape[0],
    ncols=preset_view_angles.shape[1],
    figsize=(8, 8),
)

for i in range(axs.shape[0]):
    for j in range(axs.shape[1]):
        ax = axs[i, j]
        preset_view = preset_view_angles[i, j]

        for w in [wing, wing_deformed]:
            w.draw_wireframe(
                ax=ax,
                set_axis_visibility=False if "isometric" in preset_view else None,
                show=False,
            )

        p.set_preset_3d_view_angle(preset_view)

        xres = np.diff(ax.get_xticks())[0]
        yres = np.diff(ax.get_yticks())[0]
        zres = np.diff(ax.get_zticks())[0]

        p.set_ticks(
            xres,
            xres / 4,
            yres,
            yres / 4,
            zres,
            zres / 4,
        )

        ax.xaxis.set_tick_params(color="white", which="minor")
        ax.yaxis.set_tick_params(color="white", which="minor")
        ax.zaxis.set_tick_params(color="white", which="minor")

        if preset_view == "XY" or preset_view == "-XY":
            ax.set_zticks([])
        if preset_view == "XZ" or preset_view == "-XZ":
            ax.set_yticks([])
        if preset_view == "YZ" or preset_view == "-YZ":
            ax.set_xticks([])

axs[1, 0].set_xlabel("$x_g$ [m]")
axs[1, 0].set_ylabel("$y_g$ [m]")
axs[0, 0].set_zlabel("$z_g$ [m]")
axs[0, 0].set_xticklabels([])
axs[0, 1].set_yticklabels([])
axs[0, 1].set_zticklabels([])

plt.subplots_adjust(
    left=-0.08,
    right=1.08,
    bottom=-0.08,
    top=1.08,
    wspace=-0.38,
    hspace=-0.38,
)
plt.suptitle("Illustration of an Example Prescribed Wing Displacement (Heave + Twist)")

p.show_plot(tight_layout=False)