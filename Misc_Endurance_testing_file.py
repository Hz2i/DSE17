from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *

#Endurance_class = Endurance(power_consumption=863,init_bat_capacity=42*400*3600,S=29,latitude=40,height=18000,days_from_solstice_start=0,startingtimeofday=0)

#Endurance_class.compute_endurance(time_step=20)

#Endurance_class.plot_endurance(86400*30,50)


scissorplot = ScissorPlot()

scissorplot.minimum_Sh_S(x_cg_min=0.2,x_cg_max=0.4)