from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *
from objects_detailed.Characteristics.Components_Materials import solar_panel

solar = solar_panel(0.3*0.97**2*0.95)

Endurance_class = Endurance(power_consumption=3320.010134953271,init_bat_capacity=129.68196535645288*500*3600*0.96,S=54.19851390297442+1.3,latitude=30,height=18288,solar_panel=solar,days_from_solstice_start=180,startingtimeofday=0)

Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*28,50)


#scissorplot = ScissorPlot()
#scissorplot.compute_required_coefs()

#scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)