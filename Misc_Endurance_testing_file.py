from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *

Endurance_class = Endurance(power_consumption=2593,init_bat_capacity=129.02*400*3600,S=23.81,latitude=30,height=18500,days_from_solstice_start=0,startingtimeofday=0)

Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*30,50)


#scissorplot = ScissorPlot()
#scissorplot.compute_required_coefs()

#scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)