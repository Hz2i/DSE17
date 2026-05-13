import pytest
from Objects.Characteristics.GeneralSubsystems import FlightConditionsSystem
from Objects.Characteristics.GeneralSubsystems import PayloadSystem
from Objects.Characteristics.GeneralSubsystems import ControlSystem

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

def test_control_system():
    control = ControlSystem()

    expected_mass = (control.actuator_mass + control.mass_cables + control.mass_pushrod + control.joints_mass_ratio * control.mass_pushrod) * 4
    expected_volume = (control.actuator_volume + control.volume_cable + control.volume_pushrod) * 4
    expected_power = (control.actuator_power + control.power_loss_cable) * 4

    assert control.CS_mass == pytest.approx(expected_mass)
    assert control.CS_volume == pytest.approx(expected_volume)
    assert control.CS_power_required == pytest.approx(expected_power)