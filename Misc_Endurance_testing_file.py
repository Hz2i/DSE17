from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *
from objects_detailed.Characteristics.Components_Materials import solar_panel

solar = solar_panel(0.3*0.97**2*0.95)

Endurance_class = Endurance(power_consumption=3745.538601818465,init_bat_charge=90,init_bat_capacity=250741619.83843064,S=62.93955732922287,latitude=60,height=18288,solar_panel=solar,days_from_solstice_start=282,startingtimeofday=64800)


Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*28,450)


#scissorplot = ScissorPlot()
#scissorplot.compute_required_coefs()

#scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)