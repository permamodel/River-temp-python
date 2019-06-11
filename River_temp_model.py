# coding=utf-8
# River water temperature model
# Irina Overeem, Lei Zheng, Kang Wang
# Institute of Arctic and Alpine Research
# CSDMS, Institute of Arctic and Alpine Research, University of Colorado
# Jun. 2019

#--------------------------Inputs---------------------------
# His:    incident solar (shortwave) radiation (w m-2)
# Hdl:    Downward thermal (Longwave) radiation (w m-2)
# Ta:     Air temperature (°C)
# Uw:     Wind speed (m/s)
# rh:     Relative humidity (0-1)
# P:      Surface pressure (pa)
# h:      River stage (m)
# iceoff: First day of river induation (DOY)
# iceon:  Last day of river induation (DOY)
# days:   Number of days

#--------------------------Output---------------------------
# Tw:     Water temperature (°C)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import math

reload(sys)
sys.setdefaultencoding('utf-8')
#os.chdir(os.getcwd() + "\\River_temp\\")  # direct to current dir

# River water temperature model
def River_temp_model(His, Hdl, Ta, Uw, rh, P, h, iceoff, iceon):
    # --------------------Declare constants----------------------
    c      = 4210  # Heat capacity of fresh water (J kg-1 K-1)
    rho    = 1000  # Density of fresh water (kg m-3)
    R      = 0.2  # White water albedo-McMahon et al., 2017
    sigma  = 5.67e-8  # Stefan-Boltzmann constant (W m-2 K-4)
    #               ----Riverbed module----
    dz     = 1  # Each depth step is 1 meter
    Nz     = 30  # Choose the number of depth steps
    Nt     = 365  # Choose the number of time steps
    dt     = (365*24*60*60)/Nt  # Length of each time step in seconds
    K      = 2.3  # Thermal conductivity of sediment ( W m-1 K-1 )
    density= 1250  # Density of sediment(kg m-3)
    C      = 660  # Heat capacity of sediment
    dsoil  = 30  # Depth without soil temperature variations (m)
    Tc     = -5  # Soil temperature keeps -5°C at dsoil (30 m)
    T      = Tc*np.ones([Nz, Nt-1])  # Create ground temperature matrix with Nz+1 rows, and Nt+1 columns
    r      = 0.6  # fraction of shortwave absorbed in water surface layer
    d      = 0.2  # fraction of shortwave reflected by riverbed-Web and Zhang, 1997
    f      = 0.05  # Attenuation coefficient (m-1)-Web and Zhang, 1997

    # -------------------River inundation timing------------------
    days           = len(Ta)
    DOY            = np.array(range(1, days+1))
    riverice       = np.zeros(days-1)
    riverice_index = np.where((DOY>=iceoff-1)&(DOY<=iceon-1))
    riverice[np.array(riverice_index)] = 1

    # -------------------DAILY WATER TEMPERATURE-----------------
    Tw     = np.zeros(days)
    Hsr    = np.zeros(days)
    Hlr    = np.zeros(days)
    Hlh    = np.zeros(days)
    Hsc    = np.zeros(days)
    Hht    = np.zeros(days)
    Hba    = np.zeros(days)
    Hbc    = np.zeros(days)
    DeltaH = np.zeros(days)  # Heat balance
    wt = 0

    for i in range(iceoff-1, iceon-1):
        # -------------------Solar radiation heat gain-------------------
        Hsr[i-1] = (1-R)*His[i-1]

        # -------Longwave radiation heat (Gao and Merrick, 1996)---------
        Hlr[i-1] = Hdl[i-1]-(Tw[i-1]+273.15)**4*0.97*sigma

        # ----------Evaporation heat gain (Hebert et al., 2011)----------
        # Saturated vapor pressure at the water temperature (mm Hg)
        Es = 25.374*math.exp(17.62-5271/(Tw[i-1]+273.15))
        # Atmospheric water vapour pressure (mm Hg)
        Ea = rh[i-1]*25.374*math.exp(17.62-5271/(Ta[i-1]+273.15))
        Hlh[i-1] = (3*Uw[i-1]+6)*(Ea-Es)

        # ----------Convective heat gain (Hebert et al., 2011)-----------
        Hsc[i-1] = (3.66+1.83*Uw[i-1])*(P[i-1]*0.0075/1000)*(Ta[i-1]-Tw[i-1])

        # ---------riverbed heat conduction (Web and Zhang, 1997)--------
        T[0,i-1] = Tw[i-1]  # First layer temperature (same as water temperature)
        # ------------Finite difference approximations--------------------
        depth_2D = (T[0:-3,i-2]-2*T[1:-2,i-2]+T[2:-1,i-2])/dz**2
        time_1D  = (K/density/C)*(depth_2D)
        T[1:-2,i-1] = time_1D*dt+T[1:-2,i-2]
        # ----------------------------------------------------------------
        T[dsoil/dz-1,i]= Tc  # Enforce bottom BC
        # --------Riverbed heat transfer-------
        Hht[i-1] = -K*(T[0,i-1]-T[1,i-1])/ dz
        # -----Riverbed absorbed radiation-----
        Hba[i-1] = -(1-R)*(1-r)*(1-d)*Hsr[i-1]*math.exp(-f* h[i-1])
        # --------Riverbed heat balance---------
        Hbc[i-1] = Hht[i-1]+Hba[i-1]

        # -----------------------Total heat gain-------------------------
        DeltaH[i-1] = Hsr[i-1]+Hlr[i-1]+Hlh[i-1]+Hsc[i-1]+Hbc[i-1]

        # -----------------------Water temperature-----------------------
        dwt = (1/((h[i-1]*rho*c)))*DeltaH[i-1]*riverice[i-1]
        wt  = wt+(dwt*24*3600)
        wt  = max(0,wt)
        if riverice[i-1]==0:
            wt = 0
        Tw[i] = wt

    DeltaH[(DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Tw[(    DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Hlh[(   DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Hsc[(   DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Hsr[(   DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Hlr[(   DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    Hbc[(   DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan
    h[(     DOY >= iceon - 1) | (DOY <= iceoff - 1)] = np.nan

    # ----------------------Display----------------------
    plt.figure(1)
    p1 = plt.plot(DOY, h)
    plt.ylabel('River stage (m)')
    plt.title('River stage')
    plt.grid(True)
    plt.xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365],
               ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))
    plt.axis([121, 335, 0, 3])
    plt.show()

    plt.figure(2)
    plt.plot(DOY, Hlh, linewidth='1')
    plt.plot(DOY, Hsc, linewidth='1')
    plt.plot(DOY, Hlr, linewidth='1')
    plt.plot(DOY, Hsr, linewidth='1')
    plt.plot(DOY, Hbc, linewidth='1')
    plt.plot(DOY, DeltaH, linewidth='2')
    plt.ylabel('Heat flux (W m$^{-2}$)')
    plt.title('Heat flux')
    plt.grid(True)
    plt.xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365],
               ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))
    plt.axis([121, 335, -300, 300])
    plt.legend(labels=['Latent heat flux (H$_{lh}$)', 'Convective flux (H$_{ch}$)', \
                       'Net longwave radiation (H$_{lr}$)', 'Net solar radiation (H$_{sr}$)', \
                       'Streambed heat flux (H$_{bh}$)', 'Heat balance (DeltaH$_w$)'], loc=3, ncol=2, fontsize='9')
    plt.show()

    plt.figure(3)
    plt.plot(DOY, Tw, linewidth='2')
    plt.plot(DOY, Ta, linewidth='2')
    plt.ylabel('Tw (Celsius)')
    plt.title('Temperature')
    plt.grid(True)
    plt.xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365],
               ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'))
    plt.axis([121, 335, -15, 20])
    plt.legend(labels=['River water temperature','Air temperature'], loc=1, ncol=1, fontsize='9')
    plt.show()

    return Tw

# Read atmospheric parameters
Mdata = pd.read_csv('Mfile.csv', header=0)
Mdata = np.array(Mdata)
DOY   = Mdata[:,0]  # Day of the year
His   = Mdata[:,1]  # Downward solar(shortwave)radiation(w m-2)
Hdl   = Mdata[:,2]  # Downward thermal(Longwave)radiation(w m-2)
Ta    = Mdata[:,3]  # Air temperature(°C)
Uw    = Mdata[:,4]  # Wind speed(m/s)
rh    = Mdata[:,5]  # Relative humidity(0-1)
P     = Mdata[:,6]  # Surface pressure(pa)

# Read river stage
Sdata = pd.read_csv('Sfile.csv', header=0)
Sdata = np.array(Sdata)
h     = Sdata[:,1]  # River stage (m)

# DOY of river ice off and ice on
iceoff=137
iceon =317

# Call the function and simulate river water temperature
Tw = River_temp_model(His,Hdl,Ta,Uw,rh,P,h,iceoff,iceon)
