import pytest
from Deployability import SolarPower, MissionProfile
from objects_detailed.AircraftGeneral.Aircraft import Aircraft
from Objects.Characteristics.PropulsionSystem import PropulsionSystem


# ==============================
# CODE VERIFICATION TESTS
# ==============================

# NULL VALUE TESTS
def test_deployability_null_values():
    """Test deployability calculations with zero values"""

    endurance = MissionProfile(self,latitude=30, cruise_power_total=0, Propulsion=PropulsionSystem(v_inf_cruise=0,required_thrust_cruise=0,m_TO=0,S=0,CL_max=0),D=0,p_battery_per_motor=0, solarpower = SolarPower(), Aircraft=Aircraft())

    assert endurance.climb_profile_init(plot=False,h_cloud=0,cloud_cover=0,day_of_year=0,start_time=0,time_step=0) is False
    assert endurance.climb_profile(plot=False,h_cloud=0,cloud_cover=0,day_of_year=0,start_time=0,time_step=0) is False
'''
def test_endurance_power_calc_null_values():
    """Test the internal power calculation with zero collecting area"""

    endurance = Endurance(power_consumption=0, init_bat_capacity=0)

    assert endurance.P(A=0, h=18000, lat=40, init_days_from_solstice=0, time_passed=0) == 0.0




# ORDER OF MAGNITUDE TESTS
def test_endurance_power_calc_order_of_magnitude():
    """Test order-of-magnitude for the solar power calculation"""

    endurance = Endurance(power_consumption=1000, init_bat_capacity=10000)

    power = endurance.P(A=35.0, h=18000, lat=40, init_days_from_solstice=0, time_passed=12 * 3600)

    assert 1000.0 < power < 10000.0


def test_endurance_capacity_reduction_order_of_magnitude():
    """Test order-of-magnitude for battery capacity degradation"""

    endurance = Endurance(power_consumption=1000, init_bat_capacity=10000)

    capacity_frac = endurance.reduced_capacity_frac(200)

    assert 0.85 < capacity_frac < 0.95


def test_endurance_compute_order_of_magnitude():
    """Test order-of-magnitude for endurance computation"""

    endurance = Endurance(power_consumption=0.1, init_bat_capacity=1e9)

    assert endurance.compute_endurance() is True


# ASSUMED CONSTANTS CROSS-VERIFICATION TESTS
def test_endurance_constants():
    """Test that constants used in calculations are as expected"""

    endurance = Endurance(power_consumption=1000, init_bat_capacity=10000)
    # Solar panel efficiency
    assert endurance.solar.efficiency == 0.2
    # Battery cycle limits
    assert endurance.cycle_limit_nr == 400
    assert endurance.cycle_limit_degradation == 0.2


# UNIT TESTS
def test_endurance_initialization():
    """Unit test for class initialization"""
    endurance = Endurance(power_consumption=1000, init_bat_capacity=10000, init_bat_charge=80, S=35.0, latitude=40, height=18000, startingtimeofday=0, days_from_solstice_start=0)

    assert endurance.power_consumption == 1000
    assert endurance.init_bat_capacity == 10000
    assert endurance.init_bat_charge == 80
    assert endurance.S == 35.0
    assert endurance.lat == 40
    assert endurance.h == 18000
    assert endurance.starting_timeofday == 0
    assert endurance.days_from_solstice_start == 0

def test_endurance_function_return_types():
    """Unit test for the solar power calculation method"""

    endurance = Endurance(power_consumption=1000, init_bat_capacity=10000)

    power = endurance.P(A=35.0, h=18000, lat=40, init_days_from_solstice=0, time_passed=12 * 3600)

    assert isinstance(power, float)

'''