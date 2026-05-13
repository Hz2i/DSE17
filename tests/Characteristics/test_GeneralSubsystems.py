import pytest
from Objects.Characteristics.GeneralSubsystems import FlightConditionsSystem
from Objects.Characteristics.GeneralSubsystems import PayloadSystem

def test_flight_conditions_system():
    fcs = FlightConditionsSystem()

    mass, volume, power = fcs.compute_flight_conditions_system()

    expected_mass = fcs.IMU_mass + fcs.GNSS_mass + fcs.pitot_mass
    expected_volume = fcs.IMU_volume + fcs.GNSS_volume + fcs.pitot_volume
    expected_power = fcs.IMU_power_required + fcs.GNSS_power_required + fcs.pitot_power_required

    assert mass == pytest.approx(expected_mass)
    assert volume == pytest.approx(expected_volume)
    assert power == pytest.approx(expected_power)

def test_payload_system():
    payload = PayloadSystem()

    expected_mass = payload.mass_connector + payload.mass_cables + payload.mass_mounting
    expected_volume = payload.volume
    expected_power = payload.power_loss

    assert payload.PS_mass == pytest.approx(expected_mass)
    assert payload.PS_volume == pytest.approx(expected_volume)
    assert payload.PS_power == pytest.approx(expected_power)