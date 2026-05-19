import pytest

from Objects.Characteristics.CO2 import (
	calculate_co2_emissions_MaterialWeight,
	calculate_co2_emissions_Solarpanels,
	calculate_co2_emissions_power,
	calculate_co2_emissions_computer,
	calculate_total_co2_emissions,
)


# ==============================
# CODE VERIFICATION TESTS
# ==============================


# NULL VALUE TESTS
def test_co2_material_weight_null_values():
	assert calculate_co2_emissions_MaterialWeight([], []) == 0.0

def test_co2_solarpanels_null_values():
	assert calculate_co2_emissions_Solarpanels(0.0, "GaAs") == 0.0

def test_co2_power_null_values():
	assert calculate_co2_emissions_power("battery", 0.0) == 0.0
	assert calculate_co2_emissions_power("fuelcell", 0.0) == 0.0

def test_co2_computer_null_values():
	assert calculate_co2_emissions_computer(0.0) == 0.0



# ASSUMED CONSTANTS CROSS-VERIFICATION TESTS
def test_co2_constants():
	assert calculate_co2_emissions_MaterialWeight(["Polymer"], [2.0]) == pytest.approx(10.6)
	assert calculate_co2_emissions_Solarpanels(50.0, "GaAs") == pytest.approx(3.5)
	assert calculate_co2_emissions_power("battery", 2.0) == pytest.approx(240.0)
	assert calculate_co2_emissions_power("fuelcell", 2.0) == pytest.approx(44.0)
	assert calculate_co2_emissions_computer(2.0) == pytest.approx(60.0)


# UNIT TESTS
def test_co2_total_emissions():
	total = calculate_total_co2_emissions(
		material_lst=["Polymer"],
		weight_lst=[1.0],
		energy_consumption=100.0,
		solar_panel_type="GaAs",
		battery_fuelcell="battery",
		consumption=1.0,
		weight_computer_system=1.0,
	)

	expected_total = 5.3 + 7.0 + 120.0 + 30.0
	assert total == pytest.approx(expected_total)
