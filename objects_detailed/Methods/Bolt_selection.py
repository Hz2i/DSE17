import numpy as np

def max_tensile_force(diameter, tensile_strength, safety_factor=1):
    """
    Calculate the maximum tensile force a bolt can withstand.

    Parameters:
    diameter (float): The diameter of the bolt in meters.
    tensile_strength (float): The tensile strength of the material in Pascals.

    Returns:
    float: The maximum tensile force in Newtons.
    """
    # Calculate the cross-sectional area of the bolt
    radius = diameter / 2
    area = np.pi * radius ** 2
    
    # Calculate the maximum tensile force
    max_force = area * tensile_strength / safety_factor
    
    return max_force

Mlist = [2,3,4,5,6,8,10,12]
tensile_strength = 900E6
safety_factor = 15

for diameter in Mlist:
    max_force = max_tensile_force(diameter / 1000, tensile_strength, safety_factor)  # Convert diameter to meters
    print(f"Diameter: M{diameter}, Max Tensile Force: {max_force:.2f} N")

lift_est = 200*9.81  # Estimated lift force in Newtons