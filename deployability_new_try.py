import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import ambiance as am
from Objects.Characteristics.PowerSystem_sizing import solar_incidence
from Objects.Constants import Constants
from Objects.Characteristics.PropulsionSystem import PropulsionSystem
from Objects.Characteristics.ReferenceGeometries import airfoil_e334, airfoil_e387
from Objects.AircraftGeneral.Aircraft import Aircraft
from Objects.Characteristics.Airframe import wing, empennage, fuselage, nacelles


# -----------------------------
# Daylight / Sunrise / Sunset
# -----------------------------
class LightData:
    def __init__(self, latitude_deg, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)
        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)

        self.declination_deg = self.solar_declination(self.days)
        self.daylight_hours = self.compute_daylight_hours(self.days)
        self.sunrise_list, self.sunset_list = self.compute_sunrise_sunset_lists(self.days)

    @staticmethod
    def solar_declination(day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        #return 23.43585 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))
        return 23.43585 * np.sin(np.pi*2/365 * (day_of_year + 284.0))

    def compute_daylight_hours(self, day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        delta_rad = np.deg2rad(self.solar_declination(day_of_year))
        phi = self.latitude_rad
        cos_h0 = -np.tan(phi) * np.tan(delta_rad)
        cos_h0 = np.clip(cos_h0, -1.0, 1.0)
        h0 = np.arccos(cos_h0)
        return (2.0 / 15.0) * np.rad2deg(h0)

    def compute_sunrise_sunset_lists(self, day_of_year):
        daylen = self.compute_daylight_hours(day_of_year)
        sunrise = 12.0 - 0.5 * daylen
        sunset = 12.0 + 0.5 * daylen
        return sunrise, sunset


# -----------------------------
# Solar power (daily average)
# -----------------------------
class SolarPower:
    def __init__(self, latitude_deg=30.0, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)

        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)
        #self.solar_area_m2 = float(solar_area_m2)

        self.efficiency = 0.3*0.97**2*0.95
        self.I0 = 1378.0
        self.powLimS = 270.0

        #self.incidence_angle_rad = self.calc_daily_mean_incidence_angle(self.days)
        #self.power_per_m2_W = self.calc_power_per_m2(self.days)
        #self.power_W = self.power_per_m2_W * self.solar_area_m2

    @staticmethod
    def calc_solar_declination_rad(day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        decl_deg = 23.45 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))
        return np.deg2rad(decl_deg)

    def calc_power_per_m2(self,h,h_cloud=8000,cloud_intensity_effect = 0.35, day_of_year=0,time_passed=0,starting_timeofday=0):
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
            power *= (1-cloud_intensity_effect)

        return power, sunrise_time, sunset_time, time_of_day
    
    # -----------------------------
# Mission profile
# -----------------------------
class MissionProfile:
    def __init__(self, solarpower = SolarPower(), Aircraft=Aircraft()):
        self.m_battery_guess = Aircraft.pow_store.mass
        self.gamma_guess = 2
        self.S_guess = Aircraft.solar.area
        self.solarpower = solarpower

        self.E_battery_guess = self.m_battery_guess * 500 * 3600 * 0.96  # J

        # given
        self.g = 9.81
        self.alt = 60000 * 0.3048
        self.CD0 = Aircraft.CD0
        self.AR = Aircraft.wing.AR
        self.e = Aircraft.e
        self.Pavg_climb_subsys = Aircraft.Pow_req - Aircraft.Pow_motor - 100
        self.Pavg_cruise_subsys = Aircraft.Pow_req - Aircraft.Pow_motor
        self.eta_prop = Aircraft.prop.overall_eff
        self.LD = Aircraft.CL_CD
        self.V_cruise = 25
        self.CL_max = 0.8*Aircraft.wing.CL_max
        self.nr_of_engines = Aircraft.nac.nr_of_engines

        # derived
        self.m_total_guess = Aircraft.MTOW
        #self.density_climb = self.Calc_Density_Climb()
        self.Cl_opt_climb = self.Calc_Cl_opt_climb()
        self.CD_total_climb = self.Calc_CD_total(self.Cl_opt_climb)

        # cruise
        #self.Pprop_cruise = self.Calc_Pprop_cruise()
        #self.Pavg_cruise = self.Pprop_cruise + self.Pavg_cruise_subsys

        self.Pprop_cruise = Aircraft.Pow_motor
        self.Pavg_cruise = Aircraft.Pow_req

    def Calc_Cl_opt_climb(self):
        return np.sqrt(self.CD0 * np.pi * self.AR * self.e)

    def Calc_CD_total(self, CL):
        k = 1.0 / (np.pi * self.AR * self.e)
        return self.CD0 + k * CL**2
    
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
        V = self.V_cruise
        CL = 2*self.m_total_guess*9.81/(am.Atmosphere(self.alt).density * (V)**2 * self.S_guess)
        if CL > self.CL_max:
            CL = self.CL_max
            V = np.sqrt(self.m_total_guess*9.81/self.S_guess * 2/am.Atmosphere(self.alt).density * 1/CL)
            CD = self.Calc_CD_total(CL)

        print(f"Cruise CL:{CL}, Cruise CD:{CD}, Cruise V:{V}")
        return ((self.m_total_guess * self.g) / (CL/CD)) * V / self.eta_prop

    def Calc_Pa(self,h,extra_power):
        V = self.Calc_V_Pr_climb(h)[0]
        propsystem = PropulsionSystem(T=(self.Pprop_cruise+extra_power)*self.eta_prop/V / self.nr_of_engines, velocity=V, alt=h, rpm=1000.0, torque=4.0, motor_temp=-40.0,propeller_diameter=1.5)

        return propsystem.power_required*propsystem.overall_eff * self.nr_of_engines,  propsystem.power_required * self.nr_of_engines
    
    def ROC(self,h,extra_power):
        Pr = self.Calc_V_Pr_climb(h)[1]
        Pa = self.Calc_Pa(h,extra_power)

        return (Pa-Pr)/(self.m_total_guess*9.81)
    
    def climb_profile(self,plot=False,extra_power = 1000,h_cloud=8000,cloud_intensity_effect = 0.35, day_of_year=0, start_time=0,time_step=3600):

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
        
        power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[0],h_cloud,cloud_intensity_effect,day_of_year,t[0],0)

        sunrise_time = sunrise
        sunset_time = sunset

        while timeofday < sunset or cruise == False or correctday == False:
                
            power_gen, sunrise, sunset, timeofday = self.solarpower.calc_power_per_m2(h[i],h_cloud,cloud_intensity_effect,day_of_year,t[i],0)
            power_gen *= self.S_guess

            #print(self.Calc_Pa(h[i],extra_power)[1], self.Pavg_climb_subsys)
            #print(self.Pprop_cruise, self.Pavg_cruise_subsys)

            if cruise == False:
                power_req = self.Calc_Pa(h[i],extra_power)[1]
                power_req_subsys = self.Pavg_climb_subsys

                Energy = np.append(Energy,Energy[i] + (-power_req-power_req_subsys+power_gen)*time_step)
            if cruise == True:
                power_req = self.Pavg_cruise
                Energy = np.append(Energy,Energy[i] + (-power_req+power_gen)*time_step)

            h = np.append(h,h[i] + self.ROC(h[i],extra_power)[0]*time_step)

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

MTOW = 278.8
TAS_initial = 32.5
gamma = 0.0
h_cruise = 60000*0.3048
lat = 30.0
day_margin = 0
use_batt = True
energy_delta = 0.0
DoD = 0.8
night_time = 0.0

S = 55.795
Sh_S = 0.0
Sv_S = 0.0

# Choose planform type (uncomment the required one):

# Traditional wing planform:
wing_geo = wing(S=S,A=20.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0)
fus_geo = fuselage()
emp_geo = empennage(S_h = S*Sh_S, S_v = S*Sv_S,lh=8.0,h_AR=5,v_AR=2)
nac_geo = nacelles(nr_of_engines=4)

# Flying wing planform:
# wing_geo = wing(S=S,A=25.0, qc_sweep=15.0*np.pi/180, taper=1.0, dihedral=0.0*np.pi/180.0, airfoil=airfoil_e334())
# fus_geo = fuselage(D=0.5, L1=0.2, L2=0.6, L3=0.2)
# emp_geo = empennage(S_h = 0.0, S_v = 0.0)
# nac_geo = nacelles(nr_of_engines=4)

aircraft_class = Aircraft(MTOW_guess=MTOW, TAS=TAS_initial, gamma=gamma, lat=30, day_margin=day_margin, DoD=DoD,Sh_S = Sh_S, Sv_S = Sv_S, wing=wing_geo, fus=fus_geo, emp=emp_geo, nac=nac_geo, use_batt=use_batt, energy_delta=energy_delta)

mission_profile = MissionProfile(solarpower=SolarPower(latitude_deg=lat),Aircraft=aircraft_class)

#print(mission_profile.Calc_V_Pr_climb(0))
#print(mission_profile.Calc_Pa(0))
#print(mission_profile.ROC(18500))
#mission_profile.climb_profile()
#print(mission_profile.climb_profile())

dt1 = 3600
dt2 = 3600
day_of_year = np.arange(0,360+30,30)
time_of_day = np.arange(0,86400+dt1,dt1)

deployability = np.zeros((len(day_of_year),len(time_of_day)))
sunrise = np.zeros(len(day_of_year))
sunset = np.zeros(len(day_of_year))

for i in range(len(day_of_year)):
    for j in range(len(time_of_day)):
        print(f"Day {day_of_year[i]}, time {time_of_day[j]/3600} hours")
        deployability[i][j], sunrise[i], sunset[i] = mission_profile.climb_profile(plot=False,extra_power=300,h_cloud=18500,cloud_intensity_effect = 0.35,day_of_year = day_of_year[i], start_time = time_of_day[j], time_step = dt2)
        print(deployability[i][j])

import numpy as np
import matplotlib.pyplot as plt

sunrise /= 3600
sunset /= 3600

y = time_of_day/3600
x = day_of_year

#X,Y = np.meshgrid(x,y)

plt.figure(figsize=(10, 5))

import numpy as np

import numpy as np

lower = []
upper = []

for row in deployability:

    print(row)
    # indices where NOT viable
    forbidden = np.where(row == 0)[0]
    print(forbidden)
    # no forbidden hours
    if len(forbidden) == 0:
        lower.append(np.nan)
        upper.append(np.nan)
        continue
    
    # split into contiguous forbidden regions
    splits = np.where(np.diff(forbidden) > 1)[0]
    regions = np.split(forbidden, splits + 1)

    # choose the largest forbidden block
    largest = max(regions, key=len)

    lower.append(largest[0])
    upper.append(largest[-1])


plt.plot(x,upper,color="red")
plt.plot(x,lower,color='red')
plt.fill_between(x,upper,lower,color="red",alpha=0.5,label="Not Deployable")


plt.plot(x,sunrise,linestyle="dashed",color="black",label="Sunrise Time")
plt.plot(x,sunset,linestyle="dashdot",color="black",label="Sunset Time")
#plt.contourf(X,Y, deployability.T, levels=[-0.5, 0.5, 1.5],cmap=ListedColormap(['red', 'green']))

plt.ylim(0,24)
plt.ylabel("Hour of day")
plt.xlabel("Day")
plt.legend()
plt.title("Launch availability heatmap")

plt.show()
