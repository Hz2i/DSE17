import numpy as np
from matplotlib import pyplot as plt
from Objects.Characteristics.PowerSystem_sizing import solar_incidence
from Objects.Characteristics.Components_Materials import solar_panel, battery, fuel_cell
from Objects.Constants import Constants

class Endurance:
    def __init__(self,power_consumption, init_bat_capacity,init_bat_charge=100,S=35.0,latitude = 40, height = 18000, startingtimeofday = 0,solar_panel = solar_panel(), battery= battery(),days_from_solstice_start = 0):                                 # Initialise with required values (add inputs after self, as per necessity)
        self.init_bat_capacity = init_bat_capacity
        self.cycle_limit_nr = battery.cycle_limit_nr
        self.cycle_limit_degradation = battery.cycle_limit_degradation
        self.days_from_solstice_start = days_from_solstice_start
        self.power_consumption = power_consumption
        self.solar = solar_panel
        self.lat = latitude
        self.h = height
        self.starting_timeofday = startingtimeofday
        self.init_bat_charge = init_bat_charge

        self.S = S


    def P(self,A,h,lat,init_days_from_solstice,time_passed,starting_timeofday=0):
        days_passed = (time_passed+starting_timeofday) // (24*60*60)
        days_from_solstice = init_days_from_solstice + days_passed
        time_of_day = time_passed + starting_timeofday - days_passed * 24*60*60

        daylight_analysis = solar_incidence(lat,days_from_solstice)
        daylight_analysis.daylight_cycle()
        daylight_time = daylight_analysis.daylight_time

        sunrise_time = (24*60*60-daylight_time)/2
        sunset_time = sunrise_time + daylight_time

        if sunrise_time < time_of_day and time_of_day < sunset_time:
            current_incidence = np.arccos(np.cos(lat/180*np.pi)*np.cos(15*np.pi/180*(daylight_time/2-(time_of_day-sunrise_time))/3600)*np.cos(daylight_analysis.eq_inclination)+np.sin(lat/180*np.pi)*np.sin(daylight_analysis.eq_inclination))
        else:
            current_incidence = np.pi/2

        power = A*np.minimum(self.solar.powLimS,self.solar.efficiency*Constants().solar_constant*np.cos(current_incidence))
        # How does efficiency change with height???
        return power
    
    def reduced_capacity_frac(self,nr_of_cycles):
            # SA03 battery from Amprius
            # 100% DOD to 80% SOH: 300 cycles
            # 70% DOD to 90% SOH: 700 cycles
            # Assumption: discharge between -20 and 55 deg C; charge between 10 and 55 deg C
            return (1-nr_of_cycles/self.cycle_limit_nr*self.cycle_limit_degradation)

    def compute_endurance(self,time_step=20,endurance_limit=86400*30):
        Energy = self.init_bat_capacity * self.init_bat_charge/100
        iteration = 0
        cycle = 0
        time_passed = 0

        while Energy > self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.1 and time_passed < endurance_limit:
            iteration += 1
            time_passed = iteration * time_step

            Energy += -self.power_consumption*time_step + self.P(self.S,self.h,self.lat,self.days_from_solstice_start,time_passed,self.starting_timeofday)*time_step

            cycle = (time_passed-self.starting_timeofday)// 86400

            if Energy > self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.9:
                Energy = self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.9
            else:
                continue
    

        if time_passed < endurance_limit:
            print(f'Battery dropped below 10% capacity after {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.')
            return False
        elif time_passed >= endurance_limit:
            print(f'The battery remained sufficiently charged for {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.')
            return True



    
    def plot_endurance(self,total_time, time_step):
        t = np.arange(0,total_time+time_step,time_step)
        Energy = np.zeros_like(t)
        Energy[0] = self.init_bat_capacity * self.init_bat_charge/100
        Energy_capacity = np.zeros_like(t).astype(float)
        Energy_capacity_min = np.zeros_like(t).astype(float)
        Energy_capacity[0] = self.reduced_capacity_frac(0)*0.9 * 100
        Energy_capacity_min[0] = self.reduced_capacity_frac(0)*0.1 * 100
        battery_empty = 1


        for i in range(0,len(t)-1):
            Energy[i+1] = Energy[i] - battery_empty * (self.power_consumption*time_step - self.P(self.S,self.h,self.lat,self.days_from_solstice_start,t[i],self.starting_timeofday)*time_step)
            
            cycle = (t[i+1]-self.starting_timeofday) // 86400


            Energy_capacity[i+1] = self.reduced_capacity_frac(cycle)*0.9 * 100
            Energy_capacity_min[i+1] = self.reduced_capacity_frac(cycle)*0.1 * 100

            if Energy[i+1] >= self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.9:
                Energy[i+1] = self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.9
            #elif Energy[i+1] <= self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.1:
                #Energy[i+1] = self.init_bat_capacity * self.reduced_capacity_frac(cycle)*0.1
                #battery_empty = 0
                #continue
            elif Energy[i+1] <= 0:
                battery_empty = 0
            else:
                continue
        


        t += self.starting_timeofday
        t = t/86400.0

        Energy = Energy/self.init_bat_capacity * 100

        plt.plot(t,Energy)
        plt.plot(t,Energy_capacity,color='r',linestyle='dashed')
        plt.plot(t,Energy_capacity_min,color='black',linestyle='dashed')
        plt.fill_between(t,np.ones_like(t)*100,Energy_capacity,color='mistyrose')
        plt.fill_between(t,np.zeros_like(t),Energy_capacity_min,color='mistyrose')
        plt.xscale
        plt.xlabel(f'Time [days]')
        plt.ylabel(f'Battery level [%]')
        plt.hlines(100,t[0],t[-1],colors=['r'],linestyles=['dashed'])
        plt.hlines(0,t[0],t[-1],colors=['black'],linestyles=['dashed'])
        plt.ylim(-10,110)
        plt.show()

    def plot_daylight(self,total_time,time_step):
        t = np.arange(0,total_time+time_step,time_step)
        y = np.zeros_like(t).astype(float)

        for i in range(len(t)):
            y[i] = self.P(self.S,self.h,self.lat,self.days_from_solstice_start,t[i],self.starting_timeofday)

        plt.plot(t,y)
        plt.show()
