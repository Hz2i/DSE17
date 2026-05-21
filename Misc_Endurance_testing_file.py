from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *

Endurance_class = Endurance(power_consumption=1610.28,init_bat_capacity=80.126*400*3600,S=35.25,latitude=30,height=18500,solar_panel=solar_panel(), battery=battery(),days_from_solstice_start=-14,startingtimeofday=0)

Endurance_class.compute_endurance(time_step=20)

Endurance_class.plot_endurance(86400*60,50)


#scissorplot = ScissorPlot()
#scissorplot.compute_required_coefs()

#scissorplot.plot_scissor_plot(x_cg_min=0.2,x_cg_max=0.4)

DESIGNS = {
    "conventional_batteries": {
        "name": "Conventional Wing with Batteries",
        "endurance": {
            "power_consumption_W": 1939.8220664549292,
            "solar_area_m2": 42.45884708685844,
            "S": 43.60930233284505,
            "capacity_J": 138993304.786096,
            "storage_mass_kg": 96.52312832367777,
        },
    },
    "conventional_FC": {
        "name": "Conventional Wing with Fuel Cell",
        "endurance": {
            "power_consumption_W": 1593.433131761756,
            "solar_area_m2": 34.87708221004417,
            "S": 35.83984024899324,
            "capacity_J": 56833098.40319672,
            "storage_mass_kg": 39.46742944666439,
        },
    },
    "flyingwing_batteries": {
        "name": "Flying Wing with Batteries",
        "endurance": {
            "power_consumption_W": 1610.2879321926978,
            "solar_area_m2": 35.246000270392216,
            "S": 36.18243576930812,
            "capacity_J": 115381325.54686673,
            "storage_mass_kg": 80.12592051865745,
        },
    },
    "flyingwing_FC": {
        "name": "Flying Wing with Fuel Cell",
        "endurance": {
            "power_consumption_W": 1396.9504484709755,
            "solar_area_m2": 30.576467040579246,
            "S": 31.44173199078317,
            "capacity_J": 49825135.877877116,
            "storage_mass_kg": 34.60078880408133,
        },
    },
}