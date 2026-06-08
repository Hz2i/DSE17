import pytest
import numpy as np
from Objects.Performance.Control import Control_Surface_Sizing

@pytest.fixture(scope="module")
def cs():
    """Runs AeroBuildup once at zero deflection."""
    instance = Control_Surface_Sizing()
    instance.Airplane_Geo()
    return instance

@pytest.fixture(scope="module")
def aero_zero(cs):
    """Single AeroBuildup evaluation at zero deflection."""
    return cs.vlm_run(delta_inner=0, delta_outer=0, delta_rudder=0, outer_symmetric=False)

# CLASS UNIT TESTS — initialisation and geometry

class TestInitialisation:
    """Check that default parameters are set correctly on construction."""

    def test_inner_elevon_frac_positive(self, cs):
        assert cs.inner_elevon_frac > 0.0

    def test_outer_elevon_frac_positive(self, cs):
        assert cs.outer_elevon_frac > 0.0

    def test_height_winglet_positive(self, cs):
        assert cs.height_winglet > 0.0

    def test_rudder_frac_in_range(self, cs):
        assert 0.0 < cs.rudder_frac <= 1.0

class TestGeometry:
    """Check that a valid airplane object is made."""

    def test_airplane_object_created(self, cs):
        assert cs.airplane is not None

    def test_airplane_has_two_wings(self, cs):
        assert len(cs.airplane.wings) == 2

    def test_main_wing_symmetric(self, cs):
        assert cs.airplane.wings[0].symmetric is True

    def test_winglet_symmetric(self, cs):
        assert cs.airplane.wings[1].symmetric is True

    def test_main_wing_has_five_sections(self, cs):
        assert len(cs.airplane.wings[0].xsecs) == 5

    def test_winglet_has_three_sections(self, cs):
        # root + rudder transition + tip
        assert len(cs.airplane.wings[1].xsecs) == 3

    def test_elevon_spanwise_ordering(self, cs):
        """Inner elevon start must be inboard of the connection, which is inboard of the outer end."""
        assert cs.start_inner_elevon < cs.elevon_connection < cs.end_outer_elevon

    def test_elevon_stations_within_half_span(self, cs):
        assert cs.end_outer_elevon <= cs.half_span

# ZERO-DEFLECTION TESTS

class TestNullDeflection:
    """At zero deflection the antisymmetric coefficients should be zero."""

    def test_aero_result_not_none(self, aero_zero):
        assert aero_zero is not None

    def test_roll_moment_near_zero(self, aero_zero):
        """Cl should be zero for a symmetric configuration and no deflections."""
        Cl = float(np.asarray(aero_zero["Cl"]).flat[0])
        assert abs(Cl) == 0.0

    def test_yaw_moment_near_zero(self, aero_zero):
        """Cn should be zero for a symmetric configuration and no deflections."""
        Cn = float(np.asarray(aero_zero["Cn"]).flat[0])
        assert abs(Cn) == 0.0

    def test_lift_coefficient_positive(self, aero_zero):
        """At defined 7 AoA the wing should produce positive lift."""
        CL = float(np.asarray(aero_zero["CL"]).flat[0])
        assert CL > 0.0

    def test_drag_coefficient_positive(self, aero_zero):
        CD = float(np.asarray(aero_zero["CD"]).flat[0])
        assert CD > 0.0

# ORDER OF MAGNITUDE TESTS

class TestOrderOfMagnitude:
    """Checks outputs are physically plausible."""

    def test_Cmq_order_of_magnitude(self, aero_zero):
        Cmq = float(np.asarray(aero_zero["Cmq"]).flat[0])
        assert -50.0 < Cmq < -0.1

    def test_Clp_order_of_magnitude(self, aero_zero):
        Clp = float(np.asarray(aero_zero["Clp"]).flat[0])
        assert -10.0 < Clp < -0.01

    def test_Cnr_order_of_magnitude(self, aero_zero): # FAILED but because Cnr is -0.0035
        Cnr = float(np.asarray(aero_zero["Cnr"]).flat[0])
        assert -10.0 < Cnr < -0.01

    def test_x_np_order_of_magnitude(self, aero_zero):
        x_np = float(np.asarray(aero_zero["x_np"]).flat[0])
        assert 0.0 < x_np


# VALUE LIMIT TESTS

class TestValueLimits:
    """Check physical bounds and sign conventions across the deflection range."""

    def test_inner_elevon_produces_negative_Cm(self): #FAILED
        """Positive (trailing-edge-down) inner elevon should pitch nose down (Cm < 0)."""
        cs = Control_Surface_Sizing()
        aero = cs.vlm_run(delta_inner=15, delta_outer=0, delta_rudder=0, outer_symmetric=True)
        Cm = float(np.asarray(aero["Cm"]).flat[0])
        assert Cm < 0.0

    def test_outer_elevon_antisym_produces_nonzero_Cl(self):
        """Antisymmetric outer elevon should produce a non-zero roll moment."""
        cs = Control_Surface_Sizing()
        aero = cs.vlm_run(delta_inner=0, delta_outer=15, delta_rudder=0, outer_symmetric=False)
        Cl = float(np.asarray(aero["Cl"]).flat[0])
        assert abs(Cl) > 0.0

    def test_rudder_produces_nonzero_Cn(self):
        """Rudder deflection should produce a non-zero yaw moment."""
        cs = Control_Surface_Sizing()
        aero = cs.vlm_run(delta_inner=0, delta_outer=0, delta_rudder=15, outer_symmetric=False)
        Cn = float(np.asarray(aero["Cn"]).flat[0])
        assert abs(Cn) > 0.0

    def test_Cmde_sign(self):
        """Cmde (pitch control derivative) should be negative."""
        cs = Control_Surface_Sizing()
        Cmde, _ = cs.Pitching_Coefficients()
        assert Cmde < 0.0

    def test_Clda_sign(self):
        """Clda (roll control derivative) should be non-zero."""
        cs = Control_Surface_Sizing()
        Clda, _ = cs.Rolling_Coefficients()
        assert abs(Clda) > 0.0

    def test_Cndr_sign(self):
        """Cndr (yaw control derivative) should be non-zero."""
        cs = Control_Surface_Sizing()
        Cndr, _ = cs.Yawing_Coefficients()
        assert abs(Cndr) > 0.0

    def test_OEI_deflection_within_limits(self):
        """The deflection required to counteract OEI should not exceed 30 degrees."""
        cs = Control_Surface_Sizing()
        Cndr, Cnr = cs.Yawing_Coefficients()
        y_eng = 0.6 * cs.half_span
        M_engine = 17 * y_eng
        rho_cruise = 0.116
        k = 2
        Cn_OEI = k * M_engine / (0.5 * rho_cruise * cs.op_point.velocity ** 2 * cs.S * cs.b)
        deflection_OEI = abs(np.degrees(-Cn_OEI / Cndr))
        assert deflection_OEI < 30.0, "OEI deflection requirement exceeds physical limit"

