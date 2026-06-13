from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *
from objects_detailed.Characteristics.Components_Materials import solar_panel
import pytest

solar = solar_panel(0.3*0.97**2*0.95)

Endurance_class = Endurance(power_consumption=3745.538601818465,init_bat_charge=90,init_bat_capacity=250741619.83843064,S=62.93955732922287,latitude=30,height=18288,solar_panel=solar,days_from_solstice_start=365,startingtimeofday=3600*14)


Endurance_class.compute_endurance(endurance_limit=86400*1,time_step=20)

Endurance_class.plot_endurance(86400*1,450)

