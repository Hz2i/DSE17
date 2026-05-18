import pytest
import numpy as np
from Objects.Characteristics.PowerSystem_sizing import power_storage, power_generation, solar_incidence, power_required
from Objects.Characteristics.Components_Materials import solar_panel, battery, fuel_cell

# ==============================
# CODE VERIFICATION TESTS
# ==============================

# NULL VALUE TESTS
def test_power_system_null_values():
    """Test power_system calculations with zero values"""

    storage = power_storage(power_req=0, latitude=0, days_from_solstice=0, DOD=1)
    storage.compute_weight_volume()
    assert storage.mass == 0.0
    assert storage.volume == 0.0

    generation = power_generation(power_req=0, latitude=0, days_from_solstice=0)
    generation.compute_weight_surface()
    assert generation.mass == 0.0
    assert generation.area == 0.0

    solar = solar_incidence(latitude=0, days_from_solstice=0)
    solar.daylight_cycle()
    assert solar.daylight_time > 0.0

    required = power_required(mass=0, LD=1, V_cruise=0, payload=0, payload_peak=0, payload_frac=0, margin=0)
    assert required == 0.0

# ORDER OF MAGNITUDE TESTS
def test_power_system_order_of_magnitude():
    """Test order-of-magnitude for key outputs"""

    storage = power_storage(power_req=1000, latitude=40, days_from_solstice=0, DOD=0.5)
    storage.compute_weight_volume()
    assert 0.1 < storage.mass < 100.0
    assert 0.01 < storage.volume < 10.0

    generation = power_generation(power_req=1000, latitude=40, days_from_solstice=0)
    generation.compute_weight_surface()
    assert 1.0 < generation.mass < 100.0
    assert 1.0 < generation.area < 100.0

    solar = solar_incidence(latitude=40, days_from_solstice=0)
    solar.daylight_cycle()
    assert 3600 < solar.daylight_time < 86400

    required = power_required(mass=120, LD=40, V_cruise=25, payload=100, payload_peak=150, payload_frac=0.1, margin=300)
    assert 500.0 < required < 2000.0

# ASSUMED CONSTANTS CROSS-VERIFICATION TESTS  - TODO: constant test better in components_materials test file? Or separate constants test file?
def test_power_system_constants():
    """Test that constants used in calculations are as expected"""

    battery_data = battery()
    assert battery_data.massEnergy == 144e4
    assert battery_data.volumeEnergy == 31392e5

    fuel_cell_data = fuel_cell()
    assert fuel_cell_data.massEnergy == 576e4
    assert fuel_cell_data.volumeEnergy == 972e6

    solar_panel_data = solar_panel()
    assert solar_panel_data.efficiency == 0.20
    assert solar_panel_data.powLimS == 200.0
    assert solar_panel_data.powLimM == 1000.0

    # axial tilt is 0.409 rads with error
    solar = solar_incidence(latitude=0, days_from_solstice=0)
    assert np.isclose(solar.axial_tilt, 0.409, atol=0.01)


# UNIT TESTS
def test_power_system_class_initialization():
    """Unit test for class initialization"""
    storage = power_storage(power_req=1000, latitude=40, days_from_solstice=0, DOD=0.5)
    assert storage.power_req == 1000
    assert storage.DOD == 0.475
    assert isinstance(storage.daylight_time, float)
    assert storage.batteries_used is True
    assert storage.bat is not None
    assert storage.f_c is not None
    assert storage.f_c_efficiency == 0.75

    generation = power_generation(power_req=1000, latitude=40, days_from_solstice=0)
    assert generation.power_req == 1000
    assert isinstance(generation.daylight_time, float)
    assert isinstance(generation.avg_incidence, float)
    assert generation.solar is not None

    solar = solar_incidence(latitude=40, days_from_solstice=0)
    assert isinstance(solar.axial_tilt, float)
    assert isinstance(solar.lat, float)
    assert solar.day == 0

def test_power_system_return_types():
    """Unit test for function return types on valid inputs"""

    storage = power_storage(power_req=1000, latitude=40, days_from_solstice=0, DOD=0.5)
    storage.compute_weight_volume()
    assert isinstance(storage.mass, float)
    assert isinstance(storage.volume, float)

    generation = power_generation(power_req=1000, latitude=40, days_from_solstice=0)
    generation.compute_weight_surface()
    assert isinstance(generation.mass, float)
    assert isinstance(generation.area, float)

    solar = solar_incidence(latitude=40, days_from_solstice=0)
    solar.daylight_cycle()
    assert isinstance(solar.daylight_time, float)

    required = power_required(mass=120, LD=40, V_cruise=25, payload=100, payload_peak=150, payload_frac=0.1, margin=300)
    assert isinstance(required, float)


