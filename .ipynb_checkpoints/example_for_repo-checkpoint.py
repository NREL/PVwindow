"""
The purpose of this script is to give a basic demonstration of PVwindow
Author: vmwheeler
"""

import numpy as np
import wpv
import matplotlib.pyplot as plt



Glass = wpv.Layer(4000,'nkLowFeGlass','i')
FTO = wpv.Layer(0.3,'nkFTO','c')
MAPI = wpv.Layer(0.06,'nkMAPI','c')
Ag = wpv.Layer(0.01,'nkAg','c')
TiO2lowE = wpv.Layer(0.02,'nkTiO2','c')
EVA = wpv.Layer(1500,'nkEVA','i')


layers = [Glass,FTO,MAPI,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

stack = wpv.Stack(layers)


num_lams = 200
lams = np.linspace(0.3,2.5,num=num_lams)


Rfs = []
As = []
Ts = []

iang = 0.

for lam in lams:

    [Rf,A,T] = stack.get_RAT(lam,iang)
    Rfs.append(Rf)
    As.append(A)
    Ts.append(T)


Rfs = np.array(Rfs)
As = np.array(As)
Ts = np.array(Ts)

plt.figure()
plt.plot(lams,Rfs,label=r"$R$")
plt.plot(lams,As,label=r"$A$")
plt.plot(lams,Ts,label=r"$T$")
plt.plot(lams,Rfs+As+Ts,label=r"$R+A+T$")
plt.xlabel(r"$\lambda$, micron")
plt.ylabel(r"R, A, or T")
plt.legend(loc='upper right')
plt.show()

eta = 1 #electron-hole pair extraction efficienc
Ti = 298 #inside temperature
To = 311 #outside temperature
Ui = 8.3 #overall heat transfer coefficient of layers inside active layer
Uo = 17 #same for outside
Rs = 0 #series resistence
Rsh = 1e5 #shunt resistence
AbsorberLayer = 3 #which layer is PV absorber layer


stuff = wpv.get_performance_characteristics(stack,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,iang)

print(stuff)


