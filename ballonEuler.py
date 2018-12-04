#-*-coding:utf-8-*-
"""
MBSE MitX Week 2 project
Balloon modelling
@author:Tim Gaggstatter
This models the movement of a balloon filled with helium from the bottom
to the ceiling of a building, it uses the forward Euler approximation method.
If we set the flag 'refinedModel' it'll recalculate values for air pressure
in each step and will consider a meandering motion path for the balloon.
"""

import numpy as np
from math import pi as PI
from matplotlib import pyplot as plt
import sys

refinedModel = True

#conversion
m2inch = 0.0254
ft2m = 0.3048

#  USER INPUT
#1 28inches 21.2s #2 34inches 26.2s
Di = 10.828 #Diameter [inches]
m_b_g = 2.9 #Mass [grams]
T_e = 21.2 #Experimental time [seconds]
Lf = 87.487 #Height [feet]
per_air = 0.45 #Air percentage [%]

# SYSTEM CONSTANTS
p_he = 0.1785 #Density of pure helium [kg/m^3]
p_air = 1.293 #Density of air [kg/m^3]
Cd = 0.47 #Estimated drag coeficient
g = 9.81 #Gravitational acceleration [m/s^2]


# COMPUTED VARIABLES
Lm = Lf * ft2m #Height [m]
#3.5% path length increase, meandering motion
Lm = Lm * 1.035 if refinedModel else Lm
Dm = Di * m2inch #Diameter [m]
r = (Dm/2.) #Radius [m]
m_b_kg = m_b_g / 1000. #Mass [kg]
V =( 4* PI * r**3 ) / 3. #Volume [m^3]
A = PI * r**2 #Frontal area [m^2]
#Terminal velocity [m/s]
vt = ((2*g*(p_air * V -(m_b_kg + p_he * V))) / (p_air * Cd * A) )**0.5
p_eff = (per_air * p_air) + ((1-per_air)*p_he) #Density mixed gases [kg/m^3]
m_tot = m_b_kg + ( p_eff * V ) #Total mass [kg]
alpha_C = m_tot #Alpha constant
beta_C = 0.5*p_air*Cd*A #Beta constant
gamma_C = -1 * ((p_air-p_eff)*V*g-m_b_kg*g) #Gamma constant

print ("Lm={0:.3f}, Dm={1:.3f}, m_b_kg={2}, V={3:.4f}, A={4:.4f}").format(
        Lm, Dm, m_b_kg, V, A)
print ("vt={0:.3f}, m_tot={1:.4f}, alpha_C={2:.4f}, beta_C={3:.4f}, gamma_C={4:.4f}").format(
                    vt, m_tot, alpha_C, beta_C, gamma_C)

t, T_mod = 0., 0.
ind = 1
dt = 0.001
n = int(T_e*1.2 / dt) #Max 20% longer time than experimental allowed
ddz, dz, z, error = np.zeros([n]), np.zeros([n]), np.zeros([n]), 0.

for i in range(0,n):
    ind+=1
    t += dt
    
    if refinedModel:
        #Include the fact that air pressure is non-constant, 
        #Here we'll have additional variable air pressure, beta and gamma
        p_air = 4.175e-11*(288.14-0.00649*z[i-1])**4.256 #Take care b/c 4*10e-11 is wrong
        beta_C = 0.5*p_air*Cd*A
        gamma_C = -1*((p_air-p_eff)*V*g - m_b_kg*g)        
        
        
    ddz[i] = (-1./alpha_C) * beta_C * dz[i-1]**2 -(1./alpha_C) * gamma_C
    dz[i] = dz[i-1] + ddz[i-1] * dt
    z[i] = z[i-1] + dz[i-1] * dt + 0.5 * ddz[i-1] * dt**2
    
    # Uncomment the following lines to see each iterations values
    #print ("t={0:.3f}, ddz={1:.3f}, dz={2:.3f}, z={3:.3f}, ind={4}").format(
    #                t, ddz[i], dz[i], z[i], ind)
    #print ("alpha={0:.3f}, beta={1:.3f}, gamma={2:.3f}").format(
    #                alpha_C, beta_C, gamma_C)
    
    if z[i] > Lm:   #reached the ceiling
        T_mod = t
        break

if T_mod == 0: # never reached the ceiling
    print ("Iterations exhausted before reaching the ceiling, check parameters")
    sys.exit()

print ("\nRESULTS: t={0:.3f}, ddz={1:.3f}, dz={2:.3f}, z={3:.3f}, ind={4}").format(
                    t, ddz[ind-2], dz[ind-2], z[ind-2], ind)

print ("ERROR: {0:.2f}% Model time={1} Experimental time={2}").format(
                100.-T_mod*100./T_e, T_mod, T_e)


#Plot configuration
dt_arr = range(0, (int)(t*1000), (int)(dt*1000))
ind = len(dt_arr)
plt.plot(np.array(dt_arr)/1000., z[:ind], label="Position [m]")
plt.plot(np.array(dt_arr)/1000., dz[:ind], label="Velocity [m/s]")
plt.plot(np.array(dt_arr)/1000., ddz[:ind],label="Acceleration [m/s^2]")
plt.axhline(y=Lm, color='y', linestyle='-', label="Ceiling [m]")

plt.legend(loc='upper left')
plt.xlabel("Time[s]")
#plt.ylabel("position[m], vel.[m/s] & acc.[m/s^2]")
plt.title("Approximate Solution with Forward Euler's Method 2")
plt.show()