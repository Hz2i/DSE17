from objects_detailed.Characteristics.Airframe import airframe



def heat_convection(airframe=airframe()):
    altitude = 60000 # Altitude in foot
    c_p = 1005.0  # Specific heat capacity of air in J/(kg*K)
    g = 9.81  # Gravitational acceleration in m/s^2
    k = 0.0257  # Thermal conductivity of air in W/(m*K)
    mu = 1.43226e-5  # Dynamic viscosity of air in kg/(m*s)
    rho = 0.1153  # Density of air in kg/m^3
    t = 0.01 # thickness of the material in meters
    delta_T = 70  # Difference Temperature in Kelvin
    C = 1.32 # formula coefficient
    n = 0.25 # formula exponent
    beta = 1/293 # Volumetric expansion coefficient of air in 1/K


    Gr = (rho**2 * g * beta * c_p *delta_T * t**3) / mu**2 # Grashof number
    Pr = (c_p * mu) / k # Prandtl number
    Ra = Gr * Pr # Rayleigh number
    Nu = C * Ra**n # Nusselt number
    h = (Nu * k) / t # Convective heat transfer coefficient in W/(m^2*K)

    S = airframe.S # Surface area in m^2

    S_reference =  2 * airframe.S # Reference surface area for heat transfer in m^2

    Q = h * S_reference * delta_T # Heat transfer rate in Watts
    return Q

print("Heat transfer rate (Q) in Watts:", heat_convection())
