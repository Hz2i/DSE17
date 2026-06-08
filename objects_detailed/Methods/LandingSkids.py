import numpy as np


def m_skid(nr_of_skids=6,L1=0.3,L2=0.3,L3=0.3,LH=0.4,F2_frac=0.5,W=300.0,n=2.0,cf=1.0,Ro=0.02,t=0.005,print_results=False):

    # aluminium 6061-T6
    rho = 2700 # kg/m^3
    E = 68.9e9 # Pa
    sigma_y = 276e6 # Pa
    sigma_ult = 310e6 # Pa
    tau_ult = 207e6 # Pa


    N = W*n
    F1 = N*cf / (L1+L2+L3) * L1
    F2 = N*cf / (L1+L2+L3) * L2
    F2_1 = F2*F2_frac
    F2_2 = F2 - F2_1
    F3 = N*cf / (L1+L2+L3) * L3
    # shear forces
    V1 = F1 + F2_1
    V2 = F2_2 + F3

    Ri = Ro - t # m
    I = np.pi/4 * (Ro**4 - Ri**4)

    n = 0.25 # cantilever
    k = (1/n)**0.5

    Pcr = np.pi**2 * E * I / (k * (LH))**2


    tau_const = 4/3 / (np.pi * (Ro**2 - Ri ** 2)) * (Ro**2 + Ro*Ri + Ri**2) / (Ro**2 + Ri**2)
    tau1_max = tau_const * V1
    tau2_max = tau_const * V2

    sigma_max = N / ((Ro**2-Ri**2)*np.pi*2)


    m = rho*(Ro**2-Ri**2)*np.pi*(LH*2+L1+L2+L3)*1.1

    if print_results:
        print("\n")

        if (Pcr < N/2):
            print(f'Pcr = {Pcr:.2f} N, N/2 = {N/2} N: Column will buckle')
        else:
            print(f'Pcr = {Pcr:.2f} N, N/2 = {N/2} N: Column will not buckle')

        if ((tau1_max or tau2_max) > tau_ult):
            print(f'Tau1_max = {tau1_max*1e-6:.2f} MPa, Tau2_max = {tau2_max*1e-6:.2f} MPa: Bigger than tau_ult = {tau_ult*1e-6} MPa')
        else:
            print(f'Tau1_max = {tau1_max*1e-6:.2f} MPa, Tau2_max = {tau2_max*1e-6:.2f} MPa: Smaller than tau_ult = {tau_ult*1e-6} MPa')

        if (sigma_max > sigma_y):
            print(f'sigma_max = {sigma_max*1e-6:.2f} MPa: Bigger than sigma_y = {sigma_y*1e-6:.2f} MPa')
        else:
            print(f'sigma_max = {sigma_max*1e-6:.2f} MPa: Smaller than sigma_y = {sigma_y*1e-6:.2f} MPa')


        print(f"Mass of {nr_of_skids} skids: {m*nr_of_skids:.2f} kg \n")

    return m*nr_of_skids



mass = m_skid(print_results=True)