import pytest
from Deployability import SolarPower, MissionProfile
from objects_detailed.AircraftGeneral.Aircraft import Aircraft
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.Prop_TO_CLIMB import *
from objects_detailed.Characteristics.Airframe import nacelles, fuselage, airframe

# ==============================
# Initialisation
# ==============================

MTOW = 200
TAS_initial = 25
h_cruise = 60000*0.3048
lat = 30
DoD = 0.8
S = 30

fus_geo = fuselage(D=0.0, L1=0.0, L2=0.0, L3=0.0)
nac_geo = nacelles(nr_of_engines=0, pos=[])
planform = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0, twist=-4.675, winglet_h=2.1, fus=fus_geo, nac=nac_geo, display=False, init_polar=True)

aircraft_class = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=0, lat=30, day_margin=0, DoD=DoD, airframe=planform, use_batt=True, energy_delta=0.0)

D = aircraft_class.prop.D # m (already optimized for cruise)
v_inf_cruise = aircraft_class.TAS_cruise  # m/s
required_thrust_cruise = aircraft_class.T_req  # N
m_TO = MTOW + 10.0  # kg 10 for landing gear!!
CL_max = aircraft_class.airframe.CL_max # -

propulsion = PropulsionSystem(
        v_inf_cruise=v_inf_cruise,
        required_thrust_cruise=required_thrust_cruise,
        m_TO=m_TO,
        S=S,
        CL_max=CL_max,
    )

cruise_power_total, _ = propulsion.run_full_analysis()
propulsion.D = D

cl_interp_to, cd_interp_to = build_takeoff_airfoil_interpolants(propulsion)
takeoff_rpm = solve_power_limited_takeoff_rpm(propulsion, D, cl_interp_to, cd_interp_to)
result = simulate_takeoff_roll(propulsion, D, takeoff_rpm, cl_interp_to, cd_interp_to)

TO_BATTERY_PER_MOTOR = result['power_battery_total'] / 4
CLIMB_BATTERY_PER_MOTOR = TO_BATTERY_PER_MOTOR * 0.70

mission_profile = MissionProfile(latitude=lat, cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=lat),Aircraft=aircraft_class)

# ==============================
# CODE VERIFICATION TESTS
# ==============================

# NULL VALUE TESTS
def test_deployability_null_values():
    """Test deployability class initialisation with zero values"""
    with pytest.raises(ValueError):
        missionprofile = MissionProfile(latitude=0, cruise_power_total=0, Propulsion=PropulsionSystem(v_inf_cruise=0,required_thrust_cruise=0,m_TO=0,S=0,CL_max=0),D=0,p_battery_per_motor=0, solarpower = SolarPower(), Aircraft=Aircraft())
def test_deployability_climb_init_timestep_null_values():
    """Test deployability climb init function with zero timestep"""
    with pytest.raises(ValueError):
        _, _, _ = mission_profile.climb_profile_init(plot=False,extra_power=0,h_cloud=18500,cloud_cover = 4,day_of_year = 0, start_time = 0, time_step = 0)

def test_deployability_climb_timestep_null_values():
    """Test deployability climb function with zero timestep"""
    with pytest.raises(ValueError):
        _, _, _ = mission_profile.climb_profile_init(plot=False,extra_power=0,h_cloud=18500,cloud_cover = 4,day_of_year = 0, start_time = 0, time_step = 0)

# ORDER OF MAGNITUDE TESTS
def test_deployability_CD_order_of_magnitude():
    """Test order-of-magnitude for CD calculation"""
    assert 0 < mission_profile.Calc_CD_total(0.5) < 1.0

def test_deployability_CLopt_order_of_magnitude():
    """Test order-of-magnitude for CLopt calculation"""
    assert 0 < mission_profile.Calc_Cl_opt_climb() < 2.0

def test_deployability_Pa_order_of_magnitude():
    """Test order-of-magnitude for power available calculation"""
    assert 0 < mission_profile.Calc_Pa(0,0)[0] < CLIMB_BATTERY_PER_MOTOR * 4

def test_deployability_Pprop_cruise_order_of_magnitude():
    """Test order-of-magnitude for propulsion power cruise calculation"""
    assert mission_profile.Calc_Pprop_cruise() == cruise_power_total

def test_deployability_V_order_of_magnitude():
    """Test order-of-magnitude for climb speed calculation"""
    assert aircraft_class.TAS_cruise*0.9 < mission_profile.Calc_V_Pr_climb(18288)[0] < aircraft_class.TAS_cruise*1.1

def test_deployability_Pr_climb_order_of_magnitude():
    """Test order-of-magnitude for climb power required calculation"""
    assert aircraft_class.Pow_motor*0.9 < mission_profile.Calc_V_Pr_climb(18288)[1] < aircraft_class.Pow_motor*1.1


# UNIT TESTS
def test_deployability_battery_size():
    big_battery_aircraft = aircraft_class
    big_battery_aircraft.pow_store.mass *= 10
    small_battery_aircraft = aircraft_class
    small_battery_aircraft.pow_store.mass *= 0.1

    mission_profile_big = MissionProfile(latitude=lat, cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=lat),Aircraft=big_battery_aircraft)
    mission_profile_small = MissionProfile(latitude=lat, cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=lat),Aircraft=small_battery_aircraft)
    
    assert mission_profile_big.climb_profile_init(time_step=3600)[0] > mission_profile_small.climb_profile_init(False,time_step=3600)[0] 

def test_deployability_latitudes():
    mission_profile_small_lat = MissionProfile(latitude=20, cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=20),Aircraft=aircraft_class)
    mission_profile_big_lat = MissionProfile(latitude=80, cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=80),Aircraft=aircraft_class)

    assert mission_profile_small_lat.climb_profile_init(day_of_year=0, time_step=3600)[0] > mission_profile_big_lat.climb_profile_init(day_of_year=0, time_step=3600)[0]

def test_deployability_solar_panel_size():
    mission_profile_small_solar = mission_profile
    mission_profile_small_solar.S_solar *= 0.1
    mission_profile_big_solar = mission_profile
    mission_profile_big_solar.S_solar *= 10

    assert mission_profile_small_solar.climb_profile_init(day_of_year=0, time_step=3600)[0] < mission_profile_big_solar.climb_profile_init(day_of_year=0, time_step=3600)[0]



'''
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