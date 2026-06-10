import pytest
import numpy as np

from Objects.Performance.Control.Mass_moments import Mass_moments


@pytest.fixture(scope="module")
def mm():
    """Runs AeroBuildup once at zero deflection."""
    instance = Mass_moments()
    return instance

#   Class Unit Tests

class TestInitialisation:

    def test_cg_position(self, mm):
        assert mm.cg == (-2.19,0,0)

    def test_root_position(self, mm):
        assert mm.root == (-0.5,0,0)

    def test_body_centre(self, mm):
        assert mm.body_centre == (0,0,0)

    def test_wing_geometry(self, mm):
        assert mm.AR == mm.b / mm.c
        assert mm.S == mm.b * mm.c
        assert mm.AR == mm.b ** 2 / mm.S

    def test_wing_sweep(self, mm):
        assert mm.sweep == np.radians(15)

    def test_wing_dihedral(self, mm):
        assert mm.dihedral == np.radians(2)

    def test_wing_twist(self, mm):
        assert mm.twist == np.radians(4.675)

    def test_battery_data(self, mm):
        assert mm.battery_chord_fraction == 0.20
        assert mm.battery_thickness_fraction == 0.10
        assert mm.battery_density == 2180.0

class TestInternalHelpers:

    def test_spar_direction(self, mm):
        unit_vector = mm._spar_direction()
        assert np.linalg.norm(unit_vector) == 1
        assert unit_vector[0] < 0
        assert unit_vector[1] > 0
        assert unit_vector[2] < 0

    def test_parallel_axis(self, mm):
        test1 = mm._parallel_axis(1,[1,1,1])
        test2 = mm._parallel_axis(2,[1,1,1])
        factor = test2[0,0] / test1[0,0]
        assert factor == 2

        test3 = mm._parallel_axis(1,[1,2,1])
        assert test3[0,1] == 2 * test1[0,1]

    def test_include_second_half(self, mm):
        I_test = [[1,1,1],
                  [1,1,1],
                  [1,1,1]]
        I_test_total = mm._include_second_half(I_test)
        assert I_test_total[0,1] == I_test_total[1,0] == I_test_total[1,2] == I_test_total[2,1] == 0

class TestSpar:

    def test_spar_inertia(self, mm):
        mass1 = 2
        length1 = mm.half_b
        I1 = mm.spar_inertia_fd(mass1, length1)
        I2 = mm.spar_inertia_fd(mass1 * 2, length1)
        I3 = mm.spar_inertia_fd(mass1, length1 - 1)

        # Order of magnitude
        assert I1[0,0] > 0
        assert I1[0,1] == 0
        assert I1[0,2] < 0
        assert I1[1,0] == 0
        assert I1[1,1] > 0
        assert I1[1,2] == 0
        assert I1[2,0] < 0
        assert I1[2,1] == 0
        assert I1[2,2] > 0

        # Sensitivity
        assert I2[0, 0] == I1[0, 0] * 2
        assert I3[0,0] < I1[0,0]
        assert I3[1,1] < I1[1,1]
        assert I3[2,2] < I1[2,2]
        assert I1[0,2] < I3[0,2] < 0

    # def test_battery_inertia(self, mm):
    #     mass1 =