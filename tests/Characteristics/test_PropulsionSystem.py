import pytest
import numpy as np
import ambiance as am
from Objects.Characteristics.PropulsionSystem import PropulsionSystem

# ==============================
# CODE VERIFICATION TESTS
# ==============================

# NULL VALUE TESTS
def test_propulsion_system_null_values():
    """Test propulsion system calculations with mostly zero inputs."""
    propulsion = PropulsionSystem(
        plotdata=False,
        T=0,
        velocity=25,
        alt=0,
        rpm=1000,
        torque=4,
        motor_temp=-40,
        propeller_diameter=2.5,
    )

    assert propulsion.T == 0
    assert propulsion.alt == 0
    assert propulsion.motor_eff == pytest.approx(0.94)
    # assert propulsion.calc_gearbox_eff() == pytest.approx(0.95)
    assert np.isclose(propulsion.propeller_eff, 0.94, atol=0.01)
    # assert propulsion.overall_eff == pytest.approx(0.8387566462415859)
    assert propulsion.power_required == 0.0
    # assert propulsion.lambda_adv == pytest.approx(0.19098593171027442)


# ORDER OF MAGNITUDE TESTS
def test_propulsion_system_order_of_magnitude():
    """Test order-of-magnitude for key outputs"""
    propulsion = PropulsionSystem(
        plotdata=False,
        T=100,
        velocity=25,
        alt=5000,
        rpm=1000,
        torque=4,
        motor_temp=20,
        propeller_diameter=2.5,
    )

    assert 0.8 < propulsion.motor_eff < 1.0
    # assert 0.9 < propulsion.calc_gearbox_eff() < 1.0
    assert 0.8 < propulsion.propeller_eff < 1.0
    # assert 0.7 < propulsion.overall_eff < 1.0
    assert 100.0 < propulsion.power_required < 4000.0
    # assert 0.1 < propulsion.lambda_adv < 1.0
    
# ASSUMED CONSTANTS CROSS-VERIFICATION TESTS
def test_propulsion_system_constants():
    """Test that constants used in calculations are as expected"""
    propulsion = PropulsionSystem(
        plotdata=False,
        T=100,
        velocity=25,
        alt=5000,
        rpm=1000,
        torque=4,
        motor_temp=20,
        propeller_diameter=2.5,
    )

    assert am.Atmosphere(propulsion.alt).density[0] == pytest.approx(0.736, abs=0.01)  # Air density at 5000 m
    assert am.Atmosphere(propulsion.alt).speed_of_sound[0] == pytest.approx(320.5, abs=0.1)  # Speed of sound at 5000 m


# UNIT TEST 
# TODO, unit tests when finalized