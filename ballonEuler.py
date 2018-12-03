#-*-coding:utf-8-*-
"""
MBSE MitX Week 2 project
Balloon modelling
@author:Tim Gaggstatter
"""

import numpy as np
from math import pi as PI
from matplotlib import pyplot as plt
def example():
        
    x0=0.
    y0=1.
    xf=10.
    n=101
    deltax=(xf-x0)/(n-1)
    x=np.linspace(x0,xf,n)
    y=np.zeros([n])
    y[0]=y0
    
    for i in range(1,n):
        y[i]=deltax*(-y[i-1]+np.sin(x[i-1]))+y[i-1]
        
    for i in range(n):
        print(x[i],y[i])
        
    plt.plot(x,y,'o')
    plt.xlabel("Value of x")
    plt.ylabel("Value of y")
    plt.title("Approximate Solution with Forward Euler's Method")
    plt.show()


#example()


refinedModel = False

#conversion
m2inch = 0.0254
ft2m = 0.3048

#  USER INPUT
Di = 10.828 #inches
m_b_g = 2.9
T_e = 21.2
Lf = 87.487 #feet
per_air = 0.45

# SYSTEM CONSTANTS
p_he = 0.1785
p_air = 1.293
Cd = 0.47
g = 9.81


# COMPUTED VARIABLES
Lm = Lf * ft2m
#3.5% path length increase, meandering motion
Lm = Lm * 1.035 if refinedModel else Lm
Dm = Di * m2inch
r = (Dm/2.)
m_b_kg = m_b_g / 1000.
V =( 4* PI * r**3 ) / 3.
A = PI * r**2 #Frontal area
vt = ((2*g*(p_air * V -(m_b_kg + p_he * V))) / (p_air * Cd * A) )**0.5
p_eff = (per_air * p_air) + ((1-per_air)*p_he)
m_tot = m_b_kg + ( p_eff * V )
alpha_C = m_tot #alpha constant
beta_C = 0.5*p_air*Cd*A
gamma_C = -1 * ((p_air-p_eff)*V*g-m_b_kg*g)

print ("Lm={0}, Dm={1}, m_b_kg={2}, V={3}, A={4}").format(Lm, Dm, m_b_kg, V, A)
print ("vt={0}, m_tot={1}, alpha_C={2}, beta_C={3}, gamma_C={4}").format(
                    vt, m_tot, alpha_C, beta_C, gamma_C)

t, T_mod = 0., 0.
ind = 1
dt = 0.001
n = int(T_e*1.10 / dt)
ddz, dz, z, error = np.zeros([n]), np.zeros([n]), np.zeros([n]), 0.

for i in range(0,n):
    ind+=1
    t += dt
    
    if refinedModel:
        #Include the fact that air pressure is non-constant, 
        #Here we'll have additional variable air pressure, beta and gamma
        p_air = 4.175e-11*(288.14-0.00649*z[i-1])**4.256 #Take care bc 4*10e-11 is wrong
        beta_C = 0.5*p_air*Cd*A
        gamma_C = -1*((p_air-p_eff)*V*g - m_b_kg*g)        
        
        
    ddz[i] = (-1./alpha_C) * beta_C * dz[i-1]**2 -(1./alpha_C) * gamma_C
    dz[i] = dz[i-1] + ddz[i-1] * dt
    z[i] = z[i-1] + dz[i-1] * dt + 0.5 * ddz[i-1] * dt**2
    
    
    #if t > 2.008:
    #    break
    #print t, i
    if z[i] > Lm:   #reached the ceiling
        #print ("\nRESULTS: t={0}, ddz={1}, dz={2}, z={3}, ind={4}").format(
        #            t, ddz[i], dz[i], z[i], ind, error)
        T_mod = t
        break

    
print ("\nRESULTS: t={0:.3f}, ddz={1:.3f}, dz={2:.3f}, z={3:.3f}, ind={4}").format(
                    t, ddz[ind-2], dz[ind-2], z[ind-2], ind, error)

print ("ERROR: {0:.2f}% Model time={1} Experimental time={2}").format(
                100.-T_mod*100./T_e, T_mod, T_e)


dt_arr = range(0, (int)(t*1000), (int)(dt*1000))
ind = len(dt_arr)
plt.plot(np.array(dt_arr)/1000., z[:ind])
plt.plot(np.array(dt_arr)/1000., dz[:ind])
plt.plot(np.array(dt_arr)/1000., ddz[:ind])
plt.axhline(y=Lm, color='y', linestyle='-')

plt.xlabel("Time[s.]")
plt.ylabel("position[m.], vel. & acc.")
plt.title("Approximate Solution with Forward Euler's Method 2")
plt.show()






