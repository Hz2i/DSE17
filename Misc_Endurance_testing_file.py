from Objects.Performance.Endurance import *
from Objects.Performance.ScissorPlot import *


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
            'type': battery()
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
            'type': fuel_cell()
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
            'type': battery()
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
            'type': fuel_cell()
        },
    },
}

design = "flyingwing_FC"

power_consumption = DESIGNS[design]["endurance"]["power_consumption_W"]
surface = DESIGNS[design]["endurance"]["S"]
bat = DESIGNS[design]["endurance"]["type"]
init_bat_capacity = bat.massEnergy * DESIGNS[design]["endurance"]["storage_mass_kg"]

latitudes = [-60,-45,-30,-15,0,15,30,45,60]

'''
FILE_ID = "outputs/endurance_" + design + ".txt"
out_file = open(FILE_ID, "w")

for latitude in latitudes:
    Endurance_class = Endurance(power_consumption=power_consumption,init_bat_capacity=init_bat_capacity,S=surface,latitude=latitude,height=18500,solar_panel=solar_panel(), battery=bat,days_from_solstice_start=-14,startingtimeofday=0)

    passed, time_passed = Endurance_class.compute_endurance(time_step=20)
    
    if not passed:
        print(f'At a latitude of {latitude} degrees, energy storage dropped below 10% capacity after {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.', file=out_file)
    else:
        print(f'At a latitude of {latitude} degrees, energy storage remained with sufficient capacity for {time_passed//86400} days and {(time_passed - (time_passed//86400) * 86400)/(3600):.2f} hours.', file=out_file)

'''

latitude = 30

Endurance_class = Endurance(power_consumption=power_consumption,init_bat_capacity=init_bat_capacity,S=surface,latitude=latitude,height=18500,solar_panel=solar_panel(), battery=bat,days_from_solstice_start=-14,startingtimeofday=0)

Endurance_class.plot_endurance(86400*60,50)