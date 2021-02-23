"""
This thing should find the reflectance and absorbtance as a function of 
angle of incidence. Then it fits these functions to a 4th order polynomial
because that's what EnergyPlus does for some reason.
"""

import numpy as np
from wpv import Layer,Stack
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This whole thing uses microns for length

degree = np.pi/180
inc_angles = np.linspace(0,89,num=20)*degree

'''        
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams)
lamrange = [min(lams),max(lams)]
'''

lam = 0.55

Glass = Layer(4000,'nkLowFeGlass','i')
TiO2 = Layer(0.05,'nkTiO2','c')
FTO = Layer(0.3,'nkFTO','c')
MAPI = Layer(0.5,'nkMAPI','c')
ITO = Layer(0.4,'nkITO','c')
SnO2 = Layer(0.5,'nkSnO2','c')
NiO = Layer(0.05,'nkNiO','c')
Ag = Layer(0.01,'nkAg','c')
TiO2lowE = Layer(0.02,'nkTiO2','c')
Bleach = Layer(0.5,'nkTiO2','c')
EVA = Layer(1500,'nkEVA','i')


#MAPI.plotnk(lams)


layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#layers = [MAPI]

stack = Stack(layers)

Rs=[]
As=[]
Ts=[]
for iang in inc_angles:
    [R,A,T] = stack.get_RAT(lam,iang)
    Rs.append(R)
    As.append(A)
    Ts.append(T) 

'''
for iang in inc_angles:
    Rs = []
    Ts = []
    for lam in lams:
        [R,A,T] = stack.get_RAT(lam,iang)
        Rs.append(R)
        Ts.append(T)
'''   
Rs = np.array(Rs)
As = np.array(As)
Ts = np.array(Ts)

taus = Ts/Ts[0]

def taubar(theta,t0,t1,t2,t3,t4):
    return t0+t1*np.cos(theta)+t2*np.cos(theta)**2+t3*np.cos(theta)**3+t4*np.cos(theta)**4

        
popt, pcov = curve_fit(taubar, inc_angles, taus)


plt.figure()
plt.plot(inc_angles/degree,Rs,label="$R$")
plt.plot(inc_angles/degree,As,label="$A$")
plt.plot(inc_angles/degree,Ts,label="$T$")
plt.plot(inc_angles/degree, Ts[0]*taubar(inc_angles, *popt),label="$T_{fitted}$")
plt.plot(inc_angles/degree,Rs+As+Ts,label="$R+A+T$")
plt.xlabel(r"$\theta$")
plt.legend()
plt.show()

print(taus)
