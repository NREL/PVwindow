"""
This script should make a big ol' data file with the angular and specular 
resolved T, R_f, and R_b for a PV window in a format edible by EnergyPlus
"""

import numpy as np
from wpv import Layer,Stack
import matplotlib.pyplot as plt


# This whole thing uses microns for length

degree = np.pi/180
inc_angles = np.linspace(0,89,num=20)*degree

      
num_lams = 50
lams = np.linspace(0.3,2.5,num=num_lams)
lamrange = [min(lams),max(lams)]


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
stack_b = stack.reverse()



Rfs = []
Rbs = []
Ts = []
outlams = []
outangs = []


for iang in inc_angles:

    for lam in lams:
        outlams.append(lam)
        outangs.append(iang)
        [Rf,A,T] = stack.get_RAT(lam,iang)
        Rfs.append(Rf)
        Ts.append(T)
        [Rb,A,T] = stack_b.get_RAT(lam,iang)
        Rbs.append(Rb)
        
 
outlams = np.array(outlams)
outangs = np.array(outangs)
Rfs = np.array(Rfs)
Rbs = np.array(Rbs)
Ts = np.array(Ts)

X = np.transpose([outangs/degree,outlams,Ts])
np.savetxt('./Output/T_spectral_angular.txt',X,delimiter=',',header="angle [rad], wavelength [micron], T [1]")

Y = np.transpose([outangs/degree,outlams,Rfs])
np.savetxt('./Output/Rf_spectral_angular.txt',Y,delimiter=',',header="angle [rad], wavelength [micron], Rf [1]")

Z = np.transpose([outangs/degree,outlams,Rbs])
np.savetxt('./Output/Rb_spectral_angular.txt',Z,delimiter=',',header="angle [rad], wavelength [micron], Rb [1]")


'''
plt.figure()
plt.plot(inc_angles/degree,Rs,label="$R$")
plt.plot(inc_angles/degree,Ts,label="$T$")
plt.plot(inc_angles/degree, Ts[0]*taubar(inc_angles, *popt),label="$T_{fitted}$")
plt.plot(inc_angles/degree,Rs+As+Ts,label="$R+A+T$")
plt.xlabel(r"$\theta$")
plt.legend()
plt.show()
'''