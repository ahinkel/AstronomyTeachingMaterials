#Austin Hinkel
#2023-04-22
#
#A program to simulate light curves and Radial Velocity data
#for an exoplanetary system at 10 minute intervals.
#
#I need to clean up the numbers here to make the code more readable
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd

#Constants:
G = 6.67e-11
twoPi = 2*np.pi
days2years = 1/365.24
AU2meters = 1.496e11
sigmaSB = 5.67e-8 #Stefan-Boltzmann Constant in mks units
R_sol = 6.957e8 #meters
M_sol = 1.989e30 #kg
L_sol = 3.827e26 #Watts
R_Earth = 6.371e6 #meters
M_Earth = 5.792e24 #kg


#Input Planetary System Simulation Parameters:
tmax = 75 #days 
P_planet = 21.64 #days
Mplanet = 7.26 #Earth masses
Lstar = 0.0637 #solar units 
Tstar = 4200 #Kelvin
phase = 1.57 #Radians 0 - 2pi
transitDepth = 0.0031
transitDuration = 7 #hours


#Derived Parameters:
NumPoints = tmax*24*2 #every 30 minutes (make more general later)
Mstar = Lstar ** (1.0/4)

a_planet = (Mstar * P_planet**2 * days2years**2)**(1.0/3) #AU

L_star_Watts = Lstar * L_sol
R_star_mks =  np.sqrt(L_star_Watts/(sigmaSB * 4 * np.pi * Tstar**4))
Rstar = R_star_mks/R_sol

R_planet_EarthRadii = np.sqrt(transitDepth * Rstar**2 * R_sol**2) / R_Earth
R_planet_mks = np.sqrt(transitDepth * Rstar**2 * R_sol**2)
a_planet_mks = a_planet * AU2meters #meters                                
Mstar_mks = Mstar * M_sol #kg                                         
Mplanet_mks = Mplanet * M_Earth #kg 

planetDensity = Mplanet_mks / (4 * np.pi * R_planet_mks**3 / 3.0)
AmpRV = np.sqrt(G)*Mplanet_mks/np.sqrt(a_planet_mks * Mstar_mks) #m/s


#Parameters translated to discretized time steps
pts_P_planet = int(P_planet*24*2)
pts_phase = int(pts_P_planet * phase / twoPi)
pts_transitDuration = int(transitDuration*2)



#Begin Simulation:
t = np.linspace(0, tmax, NumPoints)

#Add 'observation' noise:
noiseRV = []
noiseLC = []
for i in range(0, NumPoints):
    noiseRV.append(random.gauss(0, 0.140111))
    noiseLC.append(random.gauss(0, 0.000041))


RV_clean = AmpRV*np.sin(t*twoPi/P_planet + phase)
RV_messy = RV_clean + noiseRV

LC_start = np.ones(NumPoints) + noiseLC

LC_final = []
for i in range(0, NumPoints):
    #param = (i - pts_phase) % pts_P_planet
    if (i-pts_phase) % pts_P_planet == 0:
        #print i
        LC_final.append(LC_start[i] - transitDepth)
    elif abs((i-pts_phase) % pts_P_planet) < int(pts_transitDuration/2):
        LC_final.append(LC_start[i] - transitDepth)
    elif abs((i-pts_phase) % pts_P_planet) < pts_transitDuration:
        #param = (i-pts_phase) % 
        LC_final.append(LC_start[i] - transitDepth*(1.0 + random.random()/100))
    else:
        LC_final.append(LC_start[i])





plt.plot(t, RV_messy)
plt.show()

plt.plot(t, LC_final)
plt.show()



df = pd.DataFrame()
df['time'] = t
df['Relative_Luminosity'] = LC_final
df['Radial_Velocity'] = RV_messy

df.to_csv('exoplanetLab_v1.csv', index=False)

print " "
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "Star System Parameters and Solutions:   "
print "________________________________________________"
print "Star Mass (Solar masses):            ", Mstar
print "Star Luminosity (Solar units):       ", Lstar
print "Star Temperature (Kelvin):           ", Tstar
print "Star Radius (Solar Radii):           ", Rstar
print "Maximal Star Velocity (m/s):         ", AmpRV
print "Planet Mass (Earth masses):          ", Mplanet
print "Planet Orbital Period (days):        ", P_planet
print "Planet Orbital Semi-major Axis (AU): ", a_planet
print "Planet Radius (Earth radii):         ", R_planet_EarthRadii
print "Planet Density (kg/m^3):             ", planetDensity
print "Transit Depth:                       ", transitDepth
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " "
