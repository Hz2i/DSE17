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