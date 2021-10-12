# PVwindow: Photovoltaic Window Simulation Software

## Description

PVwindow integrates existing open source simulation codes to create a novel tool for modeling photovoltaic windows. Windows are modeled as a stack with solar irradiation impinging the outer surface with an arbitrary angle of incidence. The amount of absorption is determined in each layer. Photovoltaic layers convert some of the absorbed solar power into electrical power according to an adapted Shockley–Queisser model with a set internal quantum efficiency.[^1] The code outputs:
 - the power conversion efficiency of the photovoltaic window
 - the solar heat gain coefficient
 - The visible light transmission
 - the color of transmitted light and apparent window color window

## Implementation

PVwindow is written entirely in Python. The intensity of light that is absorbed in each layer is determined by a solution to Maxwell's equations for a stack of layers. The solution is obtained using the transfer matrix method (TMM) with an implementation that accounts for nanometer-scale layers—where interference effects are important—as well as macroscopic layers where interference effects are negligible.[^2] 

## Installation

PVwindow is a pure python code. It can be cloned into your local directory and run as long as a local installation of python is available to run it. It does depend on a few libraries, all readily available.

### Dependencies
 - numpy (https://pypi.org/project/numpy/)
 - tmm (https://pypi.org/project/tmm/)
 - colorpy (pypi repo out-dated... full python 3 version included in this repo)
 - vegas (https://pypi.org/project/vegas/)
 - pandas (https://pypi.org/project/pandas/)
 - scipy (https://pypi.org/project/scipy/)
 - matplotlib (https://pypi.org/project/matplotlib/)

## Usage

The code is organized such that a simulation script can be constructed use two types of objects: layer and stack. Layer objects store all information necessary for it to be fully described for TMM solution, its: location in the stack, thickness, and complex refractive index. A stack object defines the order of layers in the set of layers defining the window. Properties of the photovoltaic window can be processed once the stack is constructed, including:
- visible light transmission (VLT)
- solar heat gain coefficient (SHGC)
- power conversion efficiencey (PCE)
- cell temperature

### Example

The example below shows an analysis of an example window containing 8 layers of 5 different materials. First some useful tools are imported, then the library containing most of the functionality of PVwindow is imported, wpv.

```
import numpy as np
import matplotlib.pyplot as plt
import wpv
```
Layers can now be created from library of materials included in the repository. It is trivially extendable. The layers are formed into a stack.
```
Glass = wpv.Layer(4000,'nkLowFeGlass','i')
FTO = wpv.Layer(0.3,'nkFTO','c')
MAPI = wpv.Layer(0.06,'nkMAPI','c')
Ag = wpv.Layer(0.01,'nkAg','c')
TiO2lowE = wpv.Layer(0.02,'nkTiO2','c')
EVA = wpv.Layer(1500,'nkEVA','i')

layers = [Glass,FTO,MAPI,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

stack = wpv.Stack(layers)
```
The user can now, for example, pick an array of wavelengths and calculate the spectral reflectance, absorbtance, and transmittance of the window then plot the result. 
```
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
```
Further, the user can calculate numerous interesting quantities with a single function call after defining a number of parameters defining the PV absorber type and the environment surrounding the window.
```
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
```
The result is

[^1]: https://pubs.acs.org/doi/full/10.1021/acsenergylett.9b01316
[^2]: Code available here: https://github.com/sbyrnes321/tmm. Theory described here: https://arxiv.org/abs/1603.02720
