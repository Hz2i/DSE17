import pytest
from Objects.Characteristics.GeneralSubsystems import FlightConditionsSystem

def test_flight_conditions_system():
    fcs = FlightConditionsSystem()

    mass, volume, power = fcs.compute_flight_conditions_system()

    expected_mass = 0.0 + 0.055 + 0.0068 + 0.160
    expected_volume = 0.0 + 3.5 * 10**(-7) + 0.031 * 0.031 * 0.004 + 0.19 * 0.019 * 0.019
    expected_power = 0.0 + 1.5 + 0.6 + 30

    assert mass == pytest.approx(expected_mass)
    assert volume == pytest.approx(expected_volume)
    assert power == pytest.approx(expected_power)
