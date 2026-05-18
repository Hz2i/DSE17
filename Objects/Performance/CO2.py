import numpy as np

def calculate_co2_emissions_MaterialWeight(material_lst, weight_lst):
    # this definition needs the material list with its corresponding weight list to 
    # calculate the CO2 the co2 emissions for each material and then sum them up to
    # get the total co2 emissions for the product. The input is in Kg and the ouptut 
    # is in CO2 kgeq.
    if len(material_lst) != len(weight_lst):
        raise ValueError("Material list and weight list must be of the same length.")
    total_co2_emissions = 0.0
    for i in range(len(material_lst)):
        if material_lst[i] == "Carbon fiber":
            co2_emission = weight_lst[i] * 67.79
        if material_lst[i] == "Zylon":
            co2_emission = weight_lst[i] * 27.0
        if material_lst[i] == "Kevlar":
            co2_emission = weight_lst[i] * 15.22
        if material_lst[i] == "Mylar":
            co2_emission = weight_lst[i] * 3.0
        if material_lst[i] == "Polymer":   
            co2_emission = weight_lst[i] * 5.3
        else:
            raise ValueError("Invalid material type. Please choose from 'Carbon fiber', 'Zylon', 'Kevlar', 'Mylar', or 'Polymer'.")
        total_co2_emissions += co2_emission
    return total_co2_emissions


def calculate_co2_emissions_Solarpanels(energy_consumption, solar_panel_type):
    # this definition needs the energy consumption of the solar panels in kWh to calculate
    #  the CO2 emissions for the solar panels. The output is in CO2 kgeq.
    if solar_panel_type == "GaAs":
        co2_emission = energy_consumption * 0.07
    if solar_panel_type == "SC-Si":
        co2_emission = energy_consumption * 0.06
    if solar_panel_type == "perovskite":
        co2_emission = energy_consumption * 0.06
    if solar_panel_type == "Heterojunction_silicon":
        co2_emission = energy_consumption * 0.02
    else:
        raise ValueError("Invalid solar panel type. Please choose from 'GaAs', 'SC-Si', 'Heterojunction_silicon', or 'perovskite'.")
    return co2_emission


def calculate_co2_emissions_power(battery_fuelcell, consumption):
    # put battery for battery_fuelcell if the power source is battery
    #  and fuelcell if the power source is fuel cell. The consumption
    #  is in kWh for battery, but in KG H2 for fuel cell. The output is in CO2 kgeq.
    if battery_fuelcell == "battery":
        co2_emission = 120.0 * consumption
    if battery_fuelcell == "fuelcell":
        co2_emission = 22.0 * consumption
    return co2_emission

def calculate_co2_emissions_computer(weight_computer_system):
    # this definition needs the weight of the computer system in Kg
    #  to calculate the CO2 emissions for the computer system. The output is in CO2 kgeq.
    co2_emission = weight_computer_system * 30.0
    return co2_emission 

def calculate_total_co2_emissions(material_lst, weight_lst, energy_consumption, solar_panel_type, battery_fuelcell, consumption, weight_computer_system):
    # this definition needs all the inputs for the previous definitions to calculate the total
    #  CO2 emissions for the product. The output is in CO2 kgeq.
    co2_emissions_material = calculate_co2_emissions_MaterialWeight(material_lst, weight_lst)
    co2_emissions_solarpanels = calculate_co2_emissions_Solarpanels(energy_consumption, solar_panel_type)
    co2_emissions_power = calculate_co2_emissions_power(battery_fuelcell, consumption)
    co2_emissions_computer = calculate_co2_emissions_computer(weight_computer_system)
    total_co2_emissions = co2_emissions_material + co2_emissions_solarpanels + co2_emissions_power + co2_emissions_computer
    return total_co2_emissions   