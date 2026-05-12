from Objects.Performance.Endurance import *

Endurance_class = Endurance(1050,65000000,S=25,latitude=40,height=18000,days_from_solstice_start=-60,startingtimeofday=60*60*17)

Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*20,20)
