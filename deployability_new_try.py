import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import ambiance as am
from Objects.Characteristics.PowerSystem_sizing import solar_incidence
from Objects.Constants import Constants
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.Prop_TO_CLIMB import *
from objects_detailed.Characteristics.ReferenceGeometries import airfoil_e334, airfoil_e387
from objects_detailed.AircraftGeneral.Aircraft import Aircraft
from objects_detailed.Characteristics.Airframe import airframe, fuselage, nacelles
from objects_detailed.Characteristics.Components_Materials import solar_panel
from Objects.Performance.Endurance import *

solar = solar_panel(0.3*0.97**2*0.95)


# Solar Power
class SolarPower:
    def __init__(self, latitude_deg=30.0, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)

        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)
        #self.solar_area_m2 = float(solar_area_m2)

        self.efficiency = 0.3*0.97**2*0.95
        self.I0 = 1378.0
        self.powLimS = 270.0


    def calc_power_per_m2(self,h,h_cloud=8000,cloud_cover=4, day_of_year=0,time_passed=0,starting_timeofday=0):
        lat = self.latitude_rad
        days_passed = (time_passed+starting_timeofday) // (24*60*60)
        days_from_solstice = day_of_year + 10 + days_passed
        time_of_day = time_passed + starting_timeofday - days_passed * 24*60*60

        daylight_analysis = solar_incidence(self.latitude_deg,days_from_solstice)
        daylight_analysis.daylight_cycle()
        daylight_time = daylight_analysis.daylight_time

        sunrise_time = (24*60*60-daylight_time)/2
        sunset_time = sunrise_time + daylight_time

        if sunrise_time < time_of_day and time_of_day < sunset_time:
            current_incidence = np.arccos(np.cos(lat)*np.cos(15*np.pi/180*(daylight_time/2-(time_of_day-sunrise_time))/3600)*np.cos(daylight_analysis.eq_inclination)+np.sin(lat)*np.sin(daylight_analysis.eq_inclination))
        else:
            current_incidence = np.pi/2

        power = np.minimum(self.powLimS,self.efficiency*Constants().solar_constant*np.cos(current_incidence))

        if h < h_cloud:
            power *= (1-0.75*(cloud_cover/8)**3.4)

        return power, sunrise_time, sunset_time, time_of_day
    
    # -----------------------------
# Mission profile
# -----------------------------
class MissionProfile:
    def __init__(self,latitude, cruise_power_total, Propulsion,D,p_battery_per_motor, solarpower = SolarPower(), Aircraft=Aircraft()):
        self.m_battery_guess = Aircraft.pow_store.mass
        self.gamma_guess = 2
        self.S_guess = Aircraft.airframe.S
        self.solarpower = solarpower
        self.propulsion = Propulsion
        self.D = D
        self.climb_battery_per_motor = p_battery_per_motor
        self.cruise_power_total = cruise_power_total
        self.lat = latitude

        self.E_battery_guess = self.m_battery_guess * 500 * 3600 * 0.96  # J

        print(self.S_guess)

        # given
        self.g = 9.81
        self.alt = 60000 * 0.3048
        self.CD0 = Aircraft.airframe.CD0
        self.AR = Aircraft.airframe.AR
        self.K1 = Aircraft.airframe.K1
        self.K2 = Aircraft.airframe.K2
        self.Pavg_climb_subsys = Aircraft.Pow_req - Aircraft.Pow_motor - 100
        self.Pavg_cruise_subsys = Aircraft.Pow_req - Aircraft.Pow_motor
        self.LD = Aircraft.airframe.CL_CD_max
        self.V_cruise = 25.0
        self.CL_max = 0.8*Aircraft.airframe.CL_max
        self.nr_of_engines = 4

        # derived
        self.m_total_guess = Aircraft.MTOW
        self.Cl_opt_climb = self.Calc_Cl_opt_climb()
        self.CD_total_climb = self.Calc_CD_total(self.Cl_opt_climb)

        self.Pprop_cruise = Aircraft.Pow_motor
        self.Pavg_cruise = Aircraft.Pow_req

        print(self.Pavg_cruise)

    def Calc_Cl_opt_climb(self):
        return (self.K1 + np.sqrt(self.K1**2 + 12 * self.CD0 * self.K2))/(2*self.K2)

    def Calc_CD_total(self, CL):
        return self.CD0 + self.K1*CL + self.K2*CL**2

    def Calc_V_Pr_climb(self,h):
        V_opt =  np.sqrt(self.m_total_guess*9.81/self.S_guess * 2/am.Atmosphere(h).density * 1/self.Cl_opt_climb)
        if self.Cl_opt_climb > self.CL_max:
            CL = self.CL_max
            V = np.sqrt(self.m_total_guess*9.81/self.S_guess * 2/am.Atmosphere(h).density * 1/CL)
            CD = self.Calc_CD_total(CL)

            if V < self.V_cruise:
                Pr = self.m_total_guess*9.81*np.sqrt(9.81*self.m_total_guess / self.S_guess * 2/am.Atmosphere(h).density * CD**2/CL**3)
            else:
                V = self.V_cruise
                CL = 2*self.m_total_guess*9.81/(am.Atmosphere(h).density * (V)**2 * self.S_guess)
                CD = self.Calc_CD_total(CL)
                Pr = 1/2 * am.Atmosphere(h).density * V**2 * self.S_guess * CD * V
        elif V_opt < self.V_cruise:
            V = V_opt
            Pr = self.m_total_guess*9.81*np.sqrt(9.81*self.m_total_guess / self.S_guess * 2/am.Atmosphere(h).density * self.CD_total_climb**2/self.Cl_opt_climb**3)
        else:
            V = self.V_cruise
            CL = 2*self.m_total_guess*9.81/(am.Atmosphere(h).density * (V)**2 * self.S_guess)
            CD = self.Calc_CD_total(CL)
            Pr = 1/2 * am.Atmosphere(h).density * V**2 * self.S_guess * CD * V

        return V, Pr

    def Calc_Pprop_cruise(self):
        return self.cruise_power_total

    def Calc_Pa(self,h,extra_power):
        V = self.Calc_V_Pr_climb(h)[0]
        Pr = self.Calc_V_Pr_climb(h)[1]
        propsystem = evaluate_climb_state(self.propulsion,self.D,V,h,self.climb_battery_per_motor)
        return np.asarray(propsystem["power_available"]).item(), self.climb_battery_per_motor * self.nr_of_engines, (np.asarray(propsystem["power_available"]).item()-np.asarray(Pr).item())/(self.m_total_guess*9.81)

    
    def climb_profile_init(self,plot=False,extra_power = 1000,h_cloud=8000,cloud_cover = 4, day_of_year=0, start_time=0,time_step=3600):

        self.power_available_used_ROC_array = []

        self.time_cruise_start = 0

        distance_travelled = 0.0

        cruise = False
        less_than_90percent = False
        n = 0
        newtimestart=0

        t = np.array([start_time])
        h = np.array([0.0])
        Energy = np.array([self.E_battery_guess])

        i = 0
        sunset_passed = False
        correctday = True
        
        power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[0],h_cloud,cloud_cover,day_of_year,t[0],0)

        sunrise_time = sunrise
        sunset_time = sunset

        while timeofday < sunset or cruise == False or correctday == False:

            power_available, power_required, Rate_of_climb = self.Calc_Pa(h[i],extra_power)
            
            if h[i] < self.alt:
                self.power_available_used_ROC_array.append([power_available,power_required,Rate_of_climb])
            else:
                self.power_available_used_ROC_array.append([self.Pprop_cruise,self.Pprop_cruise,0])
                
            power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[i],h_cloud,cloud_cover,day_of_year,t[i],0)
            power_gen *= self.S_guess

            if cruise == False:
                power_req = power_required
                power_req_subsys = self.Pavg_climb_subsys

                Energy = np.append(Energy,Energy[i] + (-power_req-power_req_subsys+power_gen)*time_step)

                distance_travelled += time_step * self.Calc_V_Pr_climb(h[i])[0]
            if cruise == True:
                power_req = self.Pavg_cruise
                Energy = np.append(Energy,Energy[i] + (-power_req+power_gen)*time_step)

            h = np.append(h,h[i] + Rate_of_climb*time_step)

            if not less_than_90percent and Energy[i+1] >= self.E_battery_guess:
                Energy[i+1] = self.E_battery_guess
            
            if Energy[i+1] <= 0:
                timeofday += sunset
                cruise = True
                correctday = True
            
            #if less_than_90percent == False and Energy[i+1] <= self.E_battery_guess*0.95:
                #less_than_90percent = True

            #if Energy[i+1] >= self.E_battery_guess*0.95 and less_than_90percent == True:
                #Energy[i+1] = self.E_battery_guess*0.95
            
            t = np.append(t,t[i]+time_step)

            if h[i+1] >= self.alt:
                h[i+1] = self.alt
                if self.time_cruise_start == 0:
                    self.time_cruise_start = t[i+1]
                    cruise = True

            if timeofday > sunset and cruise == False and correctday == True:
                correctday = False

            if correctday == False and timeofday < sunset:
                correctday = True

            i += 1
            

        self.power_available_used_ROC_array = np.array(self.power_available_used_ROC_array)
        self.h = h

        if Energy[-1] > 0.9 * self.E_battery_guess or Energy[-2] > 0.9 * self.E_battery_guess:
            Profile_passed = 1
        else:
            Profile_passed = 0

        Energy_capacity = np.ones_like(Energy)*self.E_battery_guess*0.9
        Energy_capacity_min = np.ones_like(Energy)*self.E_battery_guess*0.1

        print("Distance travelled during climb:", distance_travelled)

        if plot:
            plt.plot(t/3600,Energy/self.E_battery_guess*100)
            plt.plot(t/3600,Energy_capacity/self.E_battery_guess*100,color='r',linestyle='dashed')
            plt.plot(t/3600,Energy_capacity_min/self.E_battery_guess*100,color='b',linestyle='dashed')
            plt.fill_between(t/3600,np.ones_like(t)*100,Energy_capacity/self.E_battery_guess*100,color='mistyrose')
            plt.fill_between(t/3600,np.zeros_like(t),Energy_capacity_min/self.E_battery_guess*100,color='mistyrose')
            plt.xscale
            plt.xlabel(f'Time [hours]')
            plt.ylabel(f'Battery level [%]')
            plt.hlines(100,t[0]/3600,t[-1]/3600,colors=['r'],linestyles=['dashed'])
            plt.hlines(0,t[0]/3600,t[-1]/3600,colors=['black'],linestyles=['dashed'])
            plt.ylim(-10,110)
            plt.show()

            plt.plot(t/3600,h)
            plt.xlabel(f'Time [hours]')
            plt.ylabel(f'Altitude [m]')
            plt.vlines(self.time_cruise_start/3600,0,h[-1],colors=["r"],linestyles=['dashed'])
            plt.show()
        

        return Profile_passed, sunrise_time, sunset_time
    

    def climb_profile(self,plot=False,extra_power = 1000,h_cloud=8000,cloud_cover = 4, day_of_year=0, start_time=0,time_step=3600):

        self.time_cruise_start = 0

        cruise = False
        less_than_90percent = False
        n = 0
        newtimestart=0

        t = np.array([start_time])
        h = np.array([0.0])
        Energy = np.array([self.E_battery_guess])

        i = 0
        sunset_passed = False
        correctday = True
        
        power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[0],h_cloud,cloud_cover,day_of_year,t[0],0)

        sunrise_time = sunrise
        sunset_time = sunset

        while timeofday < sunset or cruise == False or correctday == False:


            power_available = np.interp(h[i],self.h[:-1],self.power_available_used_ROC_array[:,0])
            power_required = np.interp(h[i],self.h[:-1],self.power_available_used_ROC_array[:,1])
            Rate_of_climb = np.interp(h[i],self.h[:-1],self.power_available_used_ROC_array[:,2])

                
            power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[i],h_cloud,cloud_cover,day_of_year,t[i],0)
            power_gen *= self.S_guess

            if cruise == False:
                power_req = power_required
                power_req_subsys = self.Pavg_climb_subsys

                Energy = np.append(Energy,Energy[i] + (-power_req-power_req_subsys+power_gen)*time_step)
            if cruise == True:
                power_req = self.Pavg_cruise
                Energy = np.append(Energy,Energy[i] + (-power_req+power_gen)*time_step)

            h = np.append(h,h[i] + Rate_of_climb*time_step)

            if not less_than_90percent and Energy[i+1] >= self.E_battery_guess:
                Energy[i+1] = self.E_battery_guess
            
            if Energy[i+1] <= 0:
                timeofday += sunset
                cruise = True
                correctday = True
            
            #if less_than_90percent == False and Energy[i+1] <= self.E_battery_guess*0.95:
                #less_than_90percent = True

            #if Energy[i+1] >= self.E_battery_guess*0.95 and less_than_90percent == True:
                #Energy[i+1] = self.E_battery_guess*0.95
            
            t = np.append(t,t[i]+time_step)

            if h[i+1] >= self.alt:
                h[i+1] = self.alt
                if self.time_cruise_start == 0:
                    self.time_cruise_start = t[i+1]
                    cruise = True

            if timeofday > sunset and cruise == False and correctday == True:
                correctday = False

            if correctday == False and timeofday < sunset:
                correctday = True

            i += 1
            

        if Energy[-1] > 0.9 * self.E_battery_guess or Energy[-2] > 0.9 * self.E_battery_guess:
            Profile_passed = 1
        else:
            Profile_passed = 0

        if Profile_passed == 1:
            days_passed = (t[-1]+start_time) // (24*60*60)
            days_from_solstice = day_of_year + 10 + days_passed
            time_of_day = t[-1] + start_time - days_passed * 24*60*60
            print(self.Pavg_cruise,self.E_battery_guess,self.S_guess,self.lat,days_from_solstice,time_of_day)
            Endurance_class = Endurance(power_consumption=self.Pavg_cruise,init_bat_capacity=self.E_battery_guess,init_bat_charge=90,S=self.S_guess,latitude=self.lat,height=18288,solar_panel=solar,days_from_solstice_start=days_from_solstice,startingtimeofday=time_of_day)
            endurance_pass = Endurance_class.compute_endurance(endurance_limit=86400*2,time_step=1800)
            if not endurance_pass:
                Profile_passed = 0

        Energy_capacity = np.ones_like(Energy)*self.E_battery_guess*0.9
        Energy_capacity_min = np.ones_like(Energy)*self.E_battery_guess*0.1

        if plot:
            plt.plot(t/3600,Energy/self.E_battery_guess*100)
            plt.plot(t/3600,Energy_capacity/self.E_battery_guess*100,color='r',linestyle='dashed')
            plt.plot(t/3600,Energy_capacity_min/self.E_battery_guess*100,color='b',linestyle='dashed')
            plt.fill_between(t/3600,np.ones_like(t)*100,Energy_capacity/self.E_battery_guess*100,color='mistyrose')
            plt.fill_between(t/3600,np.zeros_like(t),Energy_capacity_min/self.E_battery_guess*100,color='mistyrose')
            plt.xscale
            plt.xlabel(f'Time [hours]')
            plt.ylabel(f'Battery level [%]')
            plt.hlines(100,t[0]/3600,t[-1]/3600,colors=['r'],linestyles=['dashed'])
            plt.hlines(0,t[0]/3600,t[-1]/3600,colors=['black'],linestyles=['dashed'])
            plt.ylim(-10,110)
            plt.show()

            plt.plot(t/3600,h)
            plt.xlabel(f'Time [hours]')
            plt.ylabel(f'Altitude [m]')
            plt.vlines(self.time_cruise_start/3600,0,h[-1],colors=["r"],linestyles=['dashed'])
            plt.show()
        

        return Profile_passed, sunrise_time, sunset_time


MTOW = 319.4805700966566
TAS_initial = 25
gamma = 0.0
h_cruise = 60000*0.3048
lat = 30
day_margin = 0
use_batt = True
energy_delta = 0.0
DoD = 0.8
night_time = 0.0

S = 62.93886244512036
Sh_S = 0.0
Sv_S = 0.0

# Choose planform type (uncomment the required one):

# Traditional wing planform:
# Flying wing planform:
fus_geo = fuselage(D=0.0, L1=0.0, L2=0.0, L3=0.0)
nac_geo = nacelles(nr_of_engines=0, pos=[])
planform = airframe(S=S, A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0, twist=-4.675, winglet_h=2.1, fus=fus_geo, nac=nac_geo, display=False, init_polar=True)

aircraft_class = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=30, day_margin=day_margin, DoD=DoD, airframe=planform, use_batt=use_batt, energy_delta=energy_delta)


# DESIGN INPUTS
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
CLIMB_BATTERY_PER_MOTOR = TO_BATTERY_PER_MOTOR * 0.80

dt1 = 450/4
dt2 = 900
dt3 = 1
day_of_year = np.arange(0,365+dt3,dt3)
time_of_day = np.arange(0,86400+dt1,dt1)

deployability = np.zeros((3,len(day_of_year),len(time_of_day)))
sunrise = np.zeros((3,len(day_of_year)))
sunset = np.zeros((3,len(day_of_year)))
latitudes = [30,45,60]

for k in range(len(latitudes)):
    mission_profile = MissionProfile(latitude=latitudes[k], cruise_power_total=cruise_power_total, Propulsion=propulsion,D=D,p_battery_per_motor=CLIMB_BATTERY_PER_MOTOR, solarpower=SolarPower(latitude_deg=latitudes[k]),Aircraft=aircraft_class)

    _, _, _ = mission_profile.climb_profile_init(plot=False,extra_power=0,h_cloud=18500,cloud_cover = 4,day_of_year = day_of_year[0], start_time = time_of_day[0], time_step = dt2)

    for i in range(len(day_of_year)):
        for j in range(len(time_of_day)):
            print(f"Day {day_of_year[i]}, time {time_of_day[j]/3600} hours")
            deployability[k][i][j], sunrise[k][i], sunset[k][i] = mission_profile.climb_profile(plot=False,extra_power=1.2,h_cloud=18500,cloud_cover = 4,day_of_year = day_of_year[i], start_time = time_of_day[j], time_step = dt2)
            print(deployability[k][i][j])


sunrise /= 3600
sunset /= 3600

y = time_of_day/3600
x = day_of_year



from itertools import groupby

def get_windows(row, times):
    """Return list of (open_time, close_time) for all contiguous 1-runs in a row."""
    windows = []
    for val, group in groupby(enumerate(row), key=lambda x: x[1]):
        if val == 1:
            indices = [i for i, _ in group]
            windows.append((times[indices[0]], times[indices[-1]]))
    return windows  # e.g. [(t_open1, t_close1), (t_open2, t_close2)]

fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True,sharey=True,figsize=(12, 5),constrained_layout=True)




# Inspect how many windows each day has
for k in range(3):
    for i, day in enumerate(day_of_year):
        windows = get_windows(deployability[k][i], time_of_day)
        for t_open, t_close in windows:
            axs[k].fill_between(
                [day - dt3/2, day + dt3/2],          # span half a day-step either side
                t_open  / 3600 - dt1/(2*3600),
                t_close / 3600 + dt1/(2*3600),
                color="green", alpha=0.4, edgecolor="none"
            )
            #axs[k].set_title(f"{latitudes[k]}' Latitude")
    
    axs[k].plot(x,sunrise[k],linestyle="dashed",color="black")
    axs[k].plot(x,sunset[k],linestyle="dashdot",color="black")

    ax_r = axs[k].twinx()
    ax_r.set_yticks([])
    ax_r.set_ylabel(f"{latitudes[k]}",rotation=0,labelpad=15)


axs[2].plot([], [], linestyle="dashed", color="black", label="Sunrise Time")
axs[2].plot([], [], linestyle="dashdot", color="black", label="Sunset Time")
axs[2].fill_between([], [], [], color="green", alpha=0.4, label="Deployable Window")

axs[2].legend(loc="lower right")

fig.text(
    1.02, 0.5,              # x, y in figure coordinates
    "Latitude [deg]",  # label text
    rotation=270,           # vertical label (or 90 depending on direction)
    va="center",
    ha="right"
)


plt.ylim(0,24)
plt.xticks(day_of_year[::15])
fig.supylabel("Hour of day")
fig.supxlabel("Day")

#plt.savefig(f"outputs/deployability_{DESIGNS[DESIGN]['name']}.svg",bbox_inches="tight")

plt.show()