from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *

Endurance_class = Endurance(power_consumption=1610.28,init_bat_capacity=80.126*400*3600,S=35.25,latitude=30,height=18500,solar_panel=solar_panel(), battery=battery(),days_from_solstice_start=-14,startingtimeofday=0)

Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*60,50)


#scissorplot = ScissorPlot()
#scissorplot.compute_required_coefs()

#scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)