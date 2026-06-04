import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

class VnDiagram:
    def __init__(self,CL_max=1.1,CL_min=-0.8,MTOM=200,S=40,CD0=0.015,CL_alpha=2*np.pi,cbar=1.3):
        # Aircraft Properties
        self.CL_max = CL_max
        self.CL_max_flaps = self.CL_max
        self.CL_min = CL_min
        self.MTOM = MTOM # kg
        self.S = S # m2
        self.CD0 = CD0
        self.CL_alpha = CL_alpha # 1/rad
        self.c = cbar # chord length in m

        # Environment properties
        self.sealevel = Atmosphere(0)
        self.rho_0 = self.sealevel.density[0]

        self.VS1 = np.sqrt(2*self.MTOM*9.81 / (self.rho_0*self.S*self.CL_max)) # Equivalent stall speed
        self.VSF = np.sqrt(2*self.MTOM*9.81 / (self.rho_0*self.S*self.CL_max_flaps)) # Equivalent stall speed with flaps
        self.VS11 = np.sqrt(-2*self.MTOM*9.81 / (self.rho_0*self.S*self.CL_min)) # Equivalent negative stall speed

        # CS 22.337 Limit manoeuvring load factors
        self.n1 = 2.0 #+5.3
        self.n2 = 2.0 #+4.0
        self.n3 = 0.0 #-1.5
        self.n4 = -1.0 #-2.65


        # CS 22.335 Design air speeds
        self.VA = self.VS1 * np.sqrt(self.n1) # CS22
        self.VA_min = self.VS11 * np.sqrt(np.abs(self.n4)) # same as VA but for minimum CL
        self.VB = 9.5 # minimum cruise speed
        self.VC = 12.5 # max cruise speed
        self.VD = 15  # dive speed                             

    # CS 25.341(a)(5)(i)
    def gustspeed(self,altitude,dive=False):
        altitude_list = [0,4572,18288]
        gust_velocity = [17.07,13.41,6.36]
        ref_gust = np.interp(altitude,altitude_list,gust_velocity)

        H = 107 # gust gradient distance in m, CS25.341(a)(4)

        if dive:
            gust = 0.5*ref_gust*(H/107)**(1/6)
        else:
            gust = ref_gust*(H/107)**(1/6)
        return gust

    # CS 22.341 Gust load factors
    def gustload(self,altitude, V,dive=False):
        atm = Atmosphere(altitude)
        mu = (2 * self.MTOM/self.S) / (atm.density[0] * self.c * self.CL_alpha)
        k = 0.88*mu / (5.3+mu)

        U = self.gustspeed(altitude,dive)

        n_min = 1 - (k/2 * self.rho_0 * U * V * self.CL_alpha)/(self.MTOM*9.81 / self.S)
        n_max = 1 + (k/2 * self.rho_0 * U * V * self.CL_alpha)/(self.MTOM*9.81 / self.S)

        return n_min, n_max

    # CS 22.333 Flight envelope
    def n_stall(self,V,CL_stall):
        return 1/2 * self.rho_0 * V**2 * self.S * CL_stall * 1/self.MTOM / 9.81


    def n(self,V,h,CL_stall):
        n_array = []
        nm_array = []
        ng_array = []


        if CL_stall < 0:
            n1 = self.n4
            n2 = self.n3
            VA = self.VA_min
            VS1 = self.VS11
            i = 0
            sign = -1
        else:
            n1 = self.n1
            n2 = self.n2
            VA = self.VA
            VS1 = self.VS1
            i = 1
            sign = 1
        VD = self.VD

        for v in V:
            # manoeuvring
            if v <= VA:
                n_m = self.n_stall(v,CL_stall)          # max manoeuvring load factor without stalling
            elif v <= self.VC and v > VA:
                n_m = n1                                # between stall and max cruise velocity
            else:
                n_m = np.interp(v,[self.VC,VD],[n1,n2])    # interpolate between cruise and dive velocity 
            
            # gusts
            if v >= VS1 and v < self.VB:                    # gustload at altitude and airspeed (CS25 gust speeds), limited by stall load times 1.25 in accordance with CS22
                if sign == 1:
                    n_g = min(self.gustload(h,v)[i],sign*1.25*(v/VS1)**2)
                elif sign == -1:
                    n_g = max(self.gustload(h,v)[i],sign*1.25*(v/VS1)**2)
            elif v >= self.VB:                              # gustload at altitude and airspeed, interpolating between dive speed (where gust load is halved), and design gust speed.
                if sign == 1:
                    n_g = min(np.interp(v,[self.VB,VD],[self.gustload(h,self.VB)[i],self.gustload(h,VD,dive=True)[i]]),sign*1.25*(v/VS1)**2)
                else:
                    n_g = max(np.interp(v,[self.VB,VD],[self.gustload(h,self.VB)[i],self.gustload(h,VD,dive=True)[i]]),self.gustload(h,v)[i],sign*1.25*(v/VS1)**2)

            else:
                n_g = 0

            if sign == 1:
                n_array.append(max(n_g,n_m))
            else:
                n_array.append(min(n_g,n_m))

            nm_array.append(n_m)
            ng_array.append(n_g)


        
        n_array = np.array(n_array)
        nm_array = np.array(nm_array)
        ng_array = np.array(ng_array)

        return n_array,nm_array,ng_array
    


    



Vn = VnDiagram()


dV = 0.1

V = np.arange(0,Vn.VD+dV,dV)




n_max_0,nm_max_0,ng_max_0 = Vn.n(V,0,Vn.CL_max)
n_min_0,nm_min_0,ng_min_0 = Vn.n(V,0,Vn.CL_min)

n_max_100,nm_max_100,ng_max_100 = Vn.n(V,10000*0.3048,Vn.CL_max)
n_min_100,nm_min_100,ng_min_100 = Vn.n(V,10000*0.3048,Vn.CL_min)

n_max_200,nm_max_200,ng_max_200 = Vn.n(V,20000*0.3048,Vn.CL_max)
n_min_200,nm_min_200,ng_min_200 = Vn.n(V,20000*0.3048,Vn.CL_min)

n_max_300,nm_max_300,ng_max_300 = Vn.n(V,30000*0.3048,Vn.CL_max)
n_min_300,nm_min_300,ng_min_300 = Vn.n(V,30000*0.3048,Vn.CL_min)

n_max_400,nm_max_400,ng_max_400 = Vn.n(V,40000*0.3048,Vn.CL_max)
n_min_400,nm_min_400,ng_min_400 = Vn.n(V,40000*0.3048,Vn.CL_min)

n_max_500,nm_max_500,ng_max_500 = Vn.n(V,50000*0.3048,Vn.CL_max)
n_min_500,nm_min_500,ng_min_500 = Vn.n(V,50000*0.3048,Vn.CL_min)

n_max_600,nm_max_600,ng_max_600 = Vn.n(V,60000*0.3048,Vn.CL_max)
n_min_600,nm_min_600,ng_min_600 = Vn.n(V,60000*0.3048,Vn.CL_min)

plt.plot(V,n_max_0,color='r',label='0 ft')
plt.plot(V,n_min_0,color='r')
plt.vlines(V[-1],n_min_0[-1],n_max_0[-1],colors=['r'])

plt.plot(V,n_max_100,color='orange',label='10000 ft')
plt.plot(V,n_min_100,color='orange')
plt.vlines(V[-1],n_min_100[-1],n_max_100[-1],colors=['orange'])

plt.plot(V,n_max_200,color='y',label='20000 ft')
plt.plot(V,n_min_200,color='y')
plt.vlines(V[-1],n_min_200[-1],n_max_200[-1],colors=['y'])

plt.plot(V,n_max_300,color='g',label='30000 ft')
plt.plot(V,n_min_300,color='g')
plt.vlines(V[-1],n_min_300[-1],n_max_300[-1],colors=['g'])

plt.plot(V,n_max_400,color='blue',label='40000 ft')
plt.plot(V,n_min_400,color='blue')
plt.vlines(V[-1],n_min_400[-1],n_max_400[-1],colors=['blue'])

plt.plot(V,n_max_500,color='purple',label='50000 ft')
plt.plot(V,n_min_500,color='purple')
plt.vlines(V[-1],n_min_500[-1],n_max_500[-1],colors=['purple'])

plt.plot(V,n_max_600,color='black',label='60000 ft')
plt.plot(V,n_min_600,color='black')
plt.vlines(V[-1],n_min_600[-1],n_max_600[-1],colors=['black'])



# plt.plot(V,nm_max_500,color='pink',linestyle="dashed",label="Manoeuvre Loads")
# plt.plot(V,nm_min_500,color='pink',linestyle="dashed")
# plt.plot(V,ng_max_500,color='grey',linestyle="dashed",label="Gust Loads")
# plt.plot(V,ng_min_500,color='grey',linestyle="dashed")
plt.legend()

plt.show()