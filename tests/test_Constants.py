import pytest
from Objects.Constants import Constants


def test_constants_values():
    """Verify key constants and container dim/mass"""

    c = Constants()

    assert c.g == pytest.approx(9.81)
    assert c.R == pytest.approx(287.05)
    assert c.gamma == pytest.approx(1.4)
    assert c.solar_constant == 1378

    assert c.container_inner_length == pytest.approx(5.867)
    assert c.container_inner_width == pytest.approx(2.330)
    assert c.container_inner_height == pytest.approx(2.350)
    assert c.container_mass_capacity == 28280
    assert c.container_packing_efficiency == pytest.approx(0.6)