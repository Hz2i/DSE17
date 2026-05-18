import numpy as np
import ambiance as am
import matplotlib.pyplot as plt
import Objects.Characteristics.PowerSystem_sizing as ps
import Objects.Characteristics.PropulsionSystem as prop


h_range = np.linspace(0, 60000*0.3048, 100000)
density = am.Atmosphere(h_range).density
#integrate density to get mass of air in the climb envelope
int_dens = np.trapz(density, h_range)/(60000*0.3048)

#Define some properties of the climb envelope
mass_nopowersys = 46.8
mass_solar = 9.2
mass_battery_guess = 50
mass_propulsion_guess = 10
mass_total_guess = mass_nopowersys + mass_solar + mass_battery_guess
g = 9.81
alt = 60000*0.3048
climbangle_guess = 5 #degrees
S = 36.5 #m^2
CD0 = 0.010
AR = 24
e = 0.9
Pavg_climb_subsys = 300 #W
Pavg_cruise_subsys = 400+25 #W
Ppower_cruise_initial = 900 #W
Pavg_cruise = Ppower_cruise_initial + Pavg_cruise_subsys
t_night = 827*60 #s
specific_power = 400 #Wh/kg Battery specific power
ROC = 1000 #Initial guess for rate of climb in m/s
ROC_target = 0.42 #Target rate of climb in m/s
while ROC - ROC_target > 0.01:
    climbangle_guess -= 0.1
    CL_opt = np.sqrt(CD0 * np.pi * AR * e)
    print("Optimal CL:", CL_opt)
    gamma = np.radians(climbangle_guess)
    k = 1.0/(np.pi * AR * e)
    CD_total = CD0 + k * CL_opt**2
    Vclimb = np.sqrt(2*mass_total_guess*g*np.cos(gamma)/(int_dens*S*CL_opt))
    D = 0.5*CD_total*S*int_dens*Vclimb**2
    Epot = mass_total_guess * g * alt
    Edrag = D * (alt / np.sin(gamma))
    Esubsys_climb = Pavg_climb_subsys * (alt / (Vclimb * np.sin(gamma)))
    Eclimb = Epot + Edrag + Esubsys_climb
    Ecruise = Pavg_cruise * t_night
    t_climb = alt / (Vclimb * np.sin(gamma))
    battery_mass_usable = (Eclimb + Ecruise) / (specific_power * 3600) #kg
    battery_mass_total = battery_mass_usable / 0.8  # Total pack mass with margin
    mass_total = mass_nopowersys + mass_solar + battery_mass_total
    print(f"Initial Total Mass Guess: {mass_total_guess:.2f} kg, Battery Pack Mass Guess: {mass_battery_guess:.2f} kg")
    print(f"Calculated Mass Total: {mass_total:.2f} kg, Battery Mass Usable: {battery_mass_total:.2f} kg, Battery Mass Total: {battery_mass_total:.2f} kg")
        
    iteration = 0
    while abs(mass_total - mass_total_guess) > 0.1:
        mass_total_guess = mass_total
        CL_opt = np.sqrt(CD0 * np.pi * AR * e)
        print("Optimal CL:", CL_opt)
        k = 1.0/(np.pi * AR * e)
        CD_total = CD0 + k * CL_opt**2
        Vclimb = np.sqrt(2*mass_total*g*np.cos(gamma)/(int_dens*S*CL_opt))
        D = 0.5*CD_total*S*int_dens*Vclimb**2
        Epot = mass_total * g * alt
        Edrag = D * (alt / np.sin(gamma))
        Esubsys_climb = Pavg_climb_subsys * (alt / (Vclimb * np.sin(gamma)))
        Eclimb = Epot + Edrag + Esubsys_climb
        Ecruise = Pavg_cruise * t_night
        t_climb = alt / (Vclimb * np.sin(gamma))
        battery_mass_usable = (Eclimb + Ecruise) / (specific_power * 3600) #kg
        battery_mass_total = battery_mass_usable / 0.8  # Total pack mass with margin
        mass_total = mass_nopowersys + mass_solar + battery_mass_total
        iteration += 1
        print(f"Iteration {iteration}: Total Mass = {mass_total:.2f} kg, Battery Pack Mass = {battery_mass_total:.2f} kg")

    ROC = Vclimb * np.sin(gamma)
    print("Total Mass Initial Guess:", mass_nopowersys + mass_solar + mass_battery_guess, "kg")
    print("Converged Battery Mass (usable - 100%DOD):", battery_mass_usable, "kg")
    print("Converged Battery Mass (total pack - 80% DOD):", battery_mass_total, "kg")
    print("Total Mass after Convergence:", mass_total, "kg")
    print("Climb Energy:", Eclimb/1e6, "MJ")
    print("Cruise Energy:", Ecruise/1e6, "MJ")
    print("Climb Power:", Eclimb / t_climb, "W")
    print("Battery Power Required:", (Eclimb / t_climb) - Pavg_climb_subsys, "W")
    print("Cruise Power:", Ecruise / t_night, "W")
    print("Total Energy:", (Eclimb + Ecruise)/1e6, "MJ")
    print("Climb Velocity:", Vclimb, "m/s")
    print("ROC:", ROC, "m/s")
    print("Time to Climb:", t_climb / 3600, "hours")

    # L_climb = 0.5*CL_opt*S*int_dens*Vclimb**2
    # W_climb = mass_total_guess * g * np.cos(gamma)
    # print("Lift Force:", L_climb, "N")
    # print("Weight Force:", W_climb, "N")
    # print("Lift-to-Weight Ratio:", L_climb/W_climb)

final_climb_angle = climbangle_guess
T_req = D  + (mass_total * g *ROC / Vclimb)
print(f"Final Power Required for Climb: {Eclimb / t_climb:.2f} W")
print(f"Final Climb Angle: {final_climb_angle:.2f} degrees")
print(f"Final Rate of Climb: {ROC:.2f} m/s")
print(f"Final Total Mass: {mass_total:.2f} kg")
print(f"Final Battery Mass (usable - 100%DOD): {battery_mass_usable:.2f} kg")
print(f"Final Battery Mass (total pack - 80% DOD): {battery_mass_total:.2f} kg")
print(f"Final Time to Climb: {t_climb/3600:.2f} hours")
print(f"Thrust Required for Climb: {T_req:.2f} N")

#Does it fit
Pmax = T_req * Vclimb
Np = 16
D_prop = 0.55*(Pmax/(1000*Np))**(1/4)

print(f"Estimated Propeller Diameter: {D_prop:.2f} m")

power_req_cruise = Pavg_cruise