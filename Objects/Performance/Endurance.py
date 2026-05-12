import numpy as np
from matplotlib import pyplot as plt
from Objects.Characteristics.PowerSystem_sizing import solar_incidence
from Objects.Characteristics.Components_Materials import solar_panel, battery, fuel_cell
from Objects.Constants import Constants

class Endurance:
    def __init__(self,power_consumption, init_bat_capacity,S=35.0,latitude = 40, height = 18000, startingtimeofday = 0,solar_panel = solar_panel(), battery= battery(),days_from_solstice_start = 0):                                 # Initialise with required values (add inputs after self, as per necessity)
        self.init_bat_capacity = init_bat_capacity
        self.cycle_limit_nr = battery.cycle_limit_nr
        self.cycle_limit_degradation = battery.cycle_limit_degradation
        self.days_from_solstice_start = days_from_solstice_start
        self.power_consumption = power_consumption
        self.solar = solar_panel
        self.lat = latitude
        self.h = height
        self.starting_timeofday = startingtimeofday

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
            current_incidence = np.pi/2 - daylight_analysis.max_incidence*np.sin((time_of_day-sunrise_time) * np.pi / (daylight_time))
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
        Energy = self.init_bat_capacity
        iteration = 0
        time_passed = 0

        while Energy > 0 and time_passed < endurance_limit:
            iteration += 1
            time_passed = iteration * time_step

            Energy += -self.power_consumption*time_step + self.P(self.S,self.h,self.lat,self.days_from_solstice_start,time_passed,self.starting_timeofday)*time_step

            cycle = (time_passed-self.starting_timeofday) // 86400

            if Energy > self.init_bat_capacity * self.reduced_capacity_frac(cycle):
                Energy = self.init_bat_capacity * self.reduced_capacity_frac(cycle)
            else:
                continue
    

        if time_passed < endurance_limit:
            print(f'Ran out of battery after {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.')
        if time_passed >= endurance_limit:
            print(f'The battery remained sufficiently charged for {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.')

        return time_passed


    
    def plot_endurance(self,total_time, time_step):
        t = np.arange(0,total_time+time_step,time_step)
        Energy = np.zeros_like(t)
        Energy[0] = self.init_bat_capacity

        for i in range(0,len(t)-1):
            Energy[i+1] = Energy[i] - self.power_consumption*time_step + self.P(self.S,self.h,self.lat,self.days_from_solstice_start,t[i],self.starting_timeofday)*time_step
            
            cycle = (t[i+1]-self.starting_timeofday) // 86400

            if Energy[i+1] > self.init_bat_capacity * self.reduced_capacity_frac(cycle):
                Energy[i+1] = self.init_bat_capacity * self.reduced_capacity_frac(cycle)
            elif Energy[i+1] < 0:
                Energy[i+1] = 0
                break
            else:
                continue

        t += self.starting_timeofday
        t = t/86400.0

        Energy = Energy/self.init_bat_capacity * 100

        plt.plot(t,Energy)
        plt.xscale
        plt.xlabel(f'Time [days]')
        plt.ylabel(f'Battery level [%]')
        plt.hlines(100,t[0],t[-1],colors=['r'],linestyles=['dashed'])
        plt.hlines(0,t[0],t[-1],colors=['black'],linestyles=['dashed'])
        plt.ylim(-10,110)
        plt.show()