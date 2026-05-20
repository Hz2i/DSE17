def calculate_subassembly_mass(subassembly):
    """
    Calculate the mass of a subassembly based on its components.

    Parameters:
    subassembly (dict): A dictionary containing the components of the subassembly and their masses.

    Returns:
    float: The total mass of the subassembly.
    """
    total_mass = 0.0
    for component, mass in subassembly.items():
        total_mass += mass
    return total_mass

def calculate_subassembly_volume(subassembly):
    """
    Calculate the volume of a subassembly based on its outer components.

    Parameters:
    subassembly (dict): A dictionary containing the outer components of the subassembly and their volumes.

    Returns:
    float: The total volume of the subassembly.
    """
    total_volume = 0.0
    for component, volume in subassembly.items():
        total_volume += volume
    return total_volume

def calculate_degradation_rate(S_wet, S, mass_fraction_airframe, mass_fraction_energy_storage, use_batteries=True):
    """
    Calculate the degradation rate of the ahaps based on the airframe degradation and energy storage degradation.
    battery: use_batteries = True, fuel cell: use_batteries = False
    Returns:
    float: The degradation rate of the product per day
    """
    degradation_rate_airframe = 0.000216074 * mass_fraction_airframe * S_wet/S
    degradation_rate = S_wet/S * degradation_rate_airframe * mass_fraction_airframe
    
    degradation_rate_battery = 0.000150504
    degradation_rate_fuel_cell = 0.00014052
    if use_batteries:
        degradation_rate += degradation_rate_battery * mass_fraction_energy_storage
    else:
        degradation_rate += degradation_rate_fuel_cell * mass_fraction_energy_storage
    return degradation_rate

if __name__ == "__main__":
    # calculate degradation rate
    S_ratio = 2.14                 # Constant 2.14
    S = 0                          # Wing surface area in m^2
    MTOW = 0                       # MTOW in kg
    m_remain = 0                   # remaining mass in kg
    m_storage = 0                  # mass of energy storage in kg
    S_wet = S_ratio * S
    mass_fraction_airframe = m_remain/MTOW
    mass_fraction_energy_storage = m_storage/MTOW
    degradation_rate = calculate_degradation_rate(S_wet, S, mass_fraction_airframe, mass_fraction_energy_storage, use_batteries=True)
    print("Degradation rate per day:", degradation_rate * 100, "%")