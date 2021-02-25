# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:29:14 2021

@author: aduell
"""

import numpy as np
import scipy
import scipy.optimize as sci
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
'''
x = np.linspace(-5,3,num = 50)
f = lambda x:  (x+3)**4


def W(params):
    x,y,z = params
    return (x)**2+y**2- 4*z


#W = lambda x, y : ((x)**2+y**2)
print('W = ',W([3,5,2]))
print(('f = ',sci.minimize(f, -4)))
initial_guess = [8,1,6]
print('While W is ',sci.minimize(W, initial_guess, bounds = ((-2,3),(-1,4),(-10,10))))

#plt.plot(x, f(x))

def Ult(eta): 
    return(eta*W())


def eggholder(x):
   return (-(x[1] + 47) * np.sin(np.sqrt(abs(x[0]/2 + (x[1]  + 47))))-x[0] * np.sin(np.sqrt(abs(x[0] - (x[1]  + 47)))))

bounds = [(-512, 512), (-512, 512)]

x = np.arange(-512, 513)
y = np.arange(-512, 513)
xgrid, ygrid = np.meshgrid(x, y)
xy = np.stack([xgrid, ygrid])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45, -45)
ax.plot_surface(xgrid, ygrid, eggholder(xy), cmap='terrain')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('eggholder(x, y)')
plt.show()



import numpy as np
from scipy.optimize import minimize

def objective(x):
    return x[0]*x[3]*(x[0]+x[1]+x[2])+x[2]

def constraint1(x):
    return x[0]*x[1]*x[2]*x[3]-25.0

def constraint2(x):
    sum_eq = 40.0
    for i in range(4):
        sum_eq = sum_eq - x[i]**2
    return sum_eq

# initial guesses
n = 4
x0 = np.zeros(n)
x0[0] = 1.0
x0[1] = 5.0
x0[2] = 5.0
x0[3] = 1.0

# show initial objective
print('Initial Objective: ' + str(objective(x0)))

# optimize
b = (1.0,5.0)
bnds = (b, b, b, b)
con1 = {'type': 'ineq', 'fun': constraint1}
con2 = {'type': 'eq', 'fun': constraint2}
cons = ([con1,con2])
solution = minimize(objective,x0,method='SLSQP',\
                    bounds=bnds,constraints=cons)
x = solution.x

# show final objective
print('Final Objective: ' + str(objective(x)))

# print solution
print('Solution')
print('x1 = ' + str(x[0]))
print('x2 = ' + str(x[1]))
print('x3 = ' + str(x[2]))
print('x4 = ' + str(x[3]))

'''


import numpy as np
import tmm
import pandas as pd
#import tmm_vw as tmm
import matplotlib.pyplot as plt
from wpv import Layer, Stack
import scipy.interpolate, scipy.integrate, pandas, sys
import scipy
#from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e, A, ohm
import sympy
import sympy.solvers.solvers
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
import pvlib

from pvlib import pvsystem


#import numericalunits
#help(numericalunits)

# This whole thing uses microns for length

degree = np.pi/180
#inc_angle = 10.*degree
inc_angle = 0.*degree
        
num_lams = 500

lams = np.linspace(0.3,2.5,num=num_lams) #um

Tio2Thick = 0.050
FTOThick = 0.250

def Glass(Thickness = 6000):
    return Layer(Thickness,'nkLowFeGlass','i')
def TiO2(Thickness = .05):
    return Layer(Thickness,'nkTiO2','c')
def FTO(Thickness = .250):
    return Layer(Thickness,'nkFTO','c')
def MAPI(Thickness):
    return Layer(0.130,'nkMAPI','c')
def AZO(Thickness):
    return Layer(0.200,'nkAZO','c')
def ITO(Thickness): 
    return Layer(0.200,'nkITO','c')
def ITOlowE(Thickness): 
    return Layer(0.075,'nkITO','c')
def SnO2(Thickness): 
    return Layer(0.05,'nkSnO2','c')
def SnO2lowE(Thickness):
    return Layer(0.030,'nkSnO2','c')
def SnO2lowEfat(Thickness): 
    return Layer(0.050,'nkSnO2','c')
def SiO2(Thickness): 
    return Layer(0.024,'nkSiO2','c')
def NiO(Thickness): 
    return Layer(0.050,'nkNiO','c')
def Ag(Thickness): 
    return Layer(0.015,'nkAg','c')
def TiO2lowE(Thickness): 
    return Layer(0.030,'nkTiO2','c')
def TiO2lowEfat(Thickness): 
    return Layer(0.060,'nkTiO2','c')
def Bleach(Thickness): 
    return Layer(0.370,'nkBleach','c')
def ClAlPc(Thickness): 
    return Layer(0.300,'nkClAlPc','c')
def C60(Thickness): 
    return Layer(0.200,'nkC60','c')
def IR(Thickness):
    return Layer(0.060,'nkPTB7_ThIEICO_4F','c')
def MAPBr(Thickness): 
    return Layer(0.500,'nkMAPbBr3','c')
def EVA(Thickness): 
    return Layer(3000,'nkEVA','i')

Glass(400)

Thickness = [6000,0.05,0.25]
def GiveLayers(Thickness,layer1,layer2,layer3):
    return [layer1(Thickness[0]),layer2(Thickness[1]),layer3(Thickness[2])]

    
#MAPI.plotnk(lams)
#Glass.plotnk(lams)

#Triple silver low-E
#layers = [Glass,SnO2lowE,Ag,SnO2lowEfat,Ag,SnO2lowEfat,Ag,SnO2lowE]

#Double silver low-E (45,15,90,15,45)
#layers = [Glass,SnO2lowE,Ag,SnO2lowEfat,Ag,SnO2lowE]

#Double silver low-E (30,15,60,15,30)
#layers = [Glass,TiO2lowE,Ag,TiO2lowEfat,Ag,TiO2lowE]

#Single silver (30,15,30)
#layers = [Glass,TiO2lowE,Ag,TiO2lowE]

#Solar cell + Low-E on surface 4

# 50% VLT with wavelength-selective absorber, IR = 60 nm
#layers = [Glass,FTO,TiO2,IR,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
# 50% VLT with wavelength-selective absorber, C60 = 100 nm, ClAlPc = 200 nm
#layers = [Glass,FTO,TiO2,C60,ClAlPc,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

# 50% VLT with non-wavelength-selective absorber, MAPbBr3 = 500 nm

#layers = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

layers = GiveLayers(Thickness, Glass,FTO,TiO2)

# Different thicknesses of MAPI: 50% VLT = 40 nm, 25% VLT = 130 nm, 5% VLT = 370 nm, 0.5% VLT = 775 nm
#layers = [Glass,FTO,TiO2,MAPI,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
# Here's the corresponding bleached layers for 5 and 0.5%
#layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

# Different thicknesses of bleach: 5% VLT = 370 nm, 0.5% VLT = 775 nm
#layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

#layers = [Glass,FTO,TiO2,MAPI,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE,Ag,TiO2lowE,Ag,TiO2lowE]
#layers = [Glass,FTO,TiO2,MAPI,NiO,AZO,EVA,Glass,SnO2lowE,Ag,SnO2lowEfat,Ag,SnO2lowEfat,Ag,SnO2lowE]
#layers = [Glass,FTO,TiO2,Bleach,NiO,AZO,EVA,Glass,SnO2lowE,Ag,SnO2lowEfat,Ag,SnO2lowE]
#Single silver (30,15,30)
#layers = [Glass,FTO,TiO2,FTO,MAPI,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#Double silver low-E (30,15,60,15,30)
#layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowEfat,Ag,TiO2lowE]

#Solar cell + Low-E on surface 2

#layers = [Glass,FTO,TiO2,MAPI,NiO,ITO,EVA,TiO2lowE,Ag,TiO2lowE,Glass]
#layers = [Glass,TiO2lowE,Ag,TiO2lowE,Ag,TiO2lowE,Ag,TiO2lowE,EVA,Glass,ITO,NiO,MAPI,TiO2,FTO,Glass]
#layers = [Glass,TiO2lowE,Ag,TiO2lowE,Ag,TiO2lowE,EVA,Glass,ITO,NiO,Bleach,TiO2,FTO,Glass]          

#Tandem transparent solar cells
#layers = [Glass,ClAlPc,EVA,Glass]
#layers = [Glass,FTO,SnO2,IR,NiO,ITO,EVA,Glass]
#layers = [Glass,FTO,SnO2,IR,NiO,ITO,EVA,ITO,SnO2,MAPI,NiO,FTO,Glass]
#layers = [Glass,FTO,SnO2,IR,NiO,ITO,EVA,ITO,SnO2,Bleach,NiO,FTO,Glass]
#layers = [Glass,FTO,SnO2,IR,NiO,ITO,EVA,Glass]

'''
Ttests = []
for lam in lams:
    Ttests.append(np.exp(-4*np.pi*MAPI.k(lam)/lam*MAPI.d))

plt.figure()
plt.plot(lams,Ttests)
plt.show()
'''

thicks = [tmm.inf]
iorcs = ['i']
for layer in layers:
    thicks.append(layer.d)
    iorcs.append(layer.i_or_c)
thicks.append(tmm.inf)
iorcs.append('i')

thicks_bw = thicks[::-1]
iorcs_bw = iorcs[::-1]

Ts = []
Rfs = []
Rbs = []
AbsByAbsorbers = []
#EQEs2 = []
#IREQEs = []

#layerchoice = 4
layerchoice = 4
layerchoice2 = 5

for lam in lams:

    nks = [1]
    for layer in layers:
        nks.append(layer.nk(lam))
    nks.append(1)

    nks_bw = nks[::-1]
    
    front_spol = tmm.inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
    front_ppol = tmm.inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)
    back_spol = tmm.inc_tmm('s',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
    back_ppol = tmm.inc_tmm('p',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
    
    AbsByAbsorber_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice]
    AbsByAbsorber_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
    AbsByAbsorbers.append( (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2. )
    
   # EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
   # EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
   # EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
    Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
    Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
    Ts.append( (front_spol['T']+front_ppol['T']) / 2. )


Ts = np.array(Ts)
Rfs = np.array(Rfs)
Rbs = np.array(Rbs)
As = 1-Ts-Rfs
sanities = Ts+Rfs+As

AbsByAbsorbers = np.array(AbsByAbsorbers)
#EQEs2 = np.array(EQEs2)
#IREQEs=EQEs+EQEs2

# Here I calculate VLT and spit it out to the screen
def GiveVLTSpectrum(layers):
    return Stack(layers)
def GiveVLT(layers):
    VLTstack=Stack(layers)
    return VLTstack.get_visible_light_transmission(lams,inc_angle)
#VLTstack=Stack(layers)
#VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
print("VLT =",GiveVLT(layers))
#print("VLT =",VLT)
#
'''
def VLTconstraint(layers):
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT - 0.5
def FTOconstraint(layers):
    return layers[1]-.15
VLTc = {'type': 'ineq', 'func': VLTconstraint}
FTOc = {'type': 'ineq', 'func': FTOconstraint}

initial_guess = layers# [6000, .15,.25]
sci.minimize(VLT, initial_guess)
'''


X = np.transpose([lams,AbsByAbsorbers])
np.savetxt('./Output/AbsByAbsorber.txt',X,delimiter=',',header="wavelength [micron], AbsByAbsorber [1]")

Y = np.transpose([lams,Ts,Rfs,Rbs])
np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")

# This is for when there are 2 layers contributing to the AbsByAbsorber:
#Z = np.transpose([lams,IREQEs])
#np.savetxt('./Output/IREQE.txt',Z,delimiter=',',header="wavelength [micron], EQE [1]")

plt.figure()
plt.plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
plt.plot(lams,Ts,color='green',marker=None,label="$T$")
plt.plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(lams,As,color='black',marker=None,label="A")
plt.plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
# This is for when there are 2 layers contributing to the EQE:
#plt.plot(lams,IREQEs,color='gray',linestyle='--',marker=None,label="EQE")
plt.plot(lams,sanities,color='gold',marker=None,label="R+A+T")
# This is the photopic eye response
plt.plot(lams,GiveVLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
# This is the solar spectrum
 #plt.plot(lams,VLTstack.Is(lams)/max(VLTstack.Is(lams)),color='gray',marker=None,label="AM1.5")
plt.xlabel('wavelength, $\mu$m')
plt.legend(loc = 'upper right')
plt.show()



# ******************** Here I add PCE calculation *********************#

# Solar cell temperature is 300 kelvin:

c0 = 299792458 #m/s
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K              
# Tack on units
lams *= 1000 #* nm

worksheet = pandas.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
#worksheet = pandas.read_excel('/Users/lwheeler/Code/pv-window-bem/Data/astmg173.xls')
downloaded_array = np.array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)


# Tack on the appropriate units:
AM15[:,0] #*= 1000 #from um to nm
AM15[:,1] #*= W / m**2 / nm
    


Ephoton = hPlanck * c0 / lams *1e9 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck


plt.plot(Ephoton, Ts, color='magenta',marker=None,label="$T$")
plt.plot(Ephoton, Rfs,color='green',marker=None,label="$R_f$")
plt.plot(Ephoton, Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(Ephoton, AbsByAbsorbers,color='black',marker=None,label="Abs")
#plt.plot(Ephoton, DarkAeVnointerp[:,1], color='blue',marker=None,label="$Abs$")
plt.legend(loc = 'upper right')
plt.xlabel('Energy, J')
plt.show()


# Interpolate to get a continuous function which I will be able to do integrals on:

AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])#, fill_value="extrapolate")
#This requires nm scale 300-2500

# Here’s the plot, it looks correct:


#λs = np.linspace(λ_min, λ_max, num=500)
y_values = np.array([AM15interp(x) for x in lams])
#y_values = np.array([AM15interp(x) for x in λs])
plt.figure()
plt.plot(lams , y_values)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral intensity (W/m$^2$/nm)")
plt.title("Light from the sun");
plt.show()




def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton *1e9  #nm
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
    #units photons/m^2       *    1/J       *    J*s   * m/s   /   J^2 idk
    #Units = should be photons/(s*m^2)
    #Ephoton must convert to 300-2500 nm from J

#SPTest = SPhotonsPerTEA(EphotonTest)
#print('SPTest =', SPTest)



PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
print('Solar constant =',solar_constant) #/ (W/m**2))

# I need to get a continuous function of the fraction of the photons absorbed
# in the absorber layer. Here I tack on units and interpolate:

'''
# Round AbsByAbsorber to make really small numbers equal to zero
AbsByAbsorbers = AbsByAbsorbers.round(8)
AbsInterp = scipy.interpolate.interp1d(lams, AbsByAbsorbers)#, fill_value="extrapolate")
EInterp = scipy.interpolate.interp1d(Ephoton, AbsByAbsorbers)
'''

def GiveAbsInterp(AbsByAbsorbers):
    Abs = AbsByAbsorbers.round(8)
    return scipy.interpolate.interp1d(lams, Abs)#, fill_value="extrapolate")

def GiveEInterp(AbsByAbsorbers):
    AbsByAbsorbers = AbsByAbsorbers.round(8)
    return scipy.interpolate.interp1d(Ephoton, AbsByAbsorbers)



#λs = np.linspace(λ_min, λ_max, num=500)
#Abs_values = np.array([AbsInterp(x) for x in λs])
Abs_values = np.array([GiveAbsInterp(AbsByAbsorbers)(x) for x in lams])
plt.figure()
#plt.plot(λs, Abs_values )
plt.plot(lams , Abs_values )
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorptance")
plt.title("Absorbed in absorber layer");
plt.show()

#AbsByAbsorbers

#def AbsVsE(Ephoton):
#    λ = hPlanck * c0 / Ephoton 
#    return AbsInterp(λ)
#return AbsInterp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2)
    
# def EQE(eta,AbsInterp):
#    lam = hPlanck * c0 / Ephoton 
#    return eta * AbsInterp(lam)




#Tcell = 300 #* K
Rs = .02 #* ohm #lumped series resistance parameter
Rsh = 10 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
q = 1.602176634e-19 #elementary charge C
Absorbed = GiveEInterp(AbsByAbsorbers)
#Absorbed = EInterp
#Absorbed = AbsInterp




# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

def RR0(eta,Absorbed,Tcell =300):
    integrand = lambda E : eta * Absorbed(E) * (E)**2 / (np.exp(E / (kB * Tcell)) - 1)
    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral# / 1.60218e-19 #J/eV
#units = 1/(s*m**2)


def Generated(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(E) * SPhotonsPerTEA(E)
#    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
#units 1/(s*m**2)


#RR = RR0(eta,Absorbed,Tcell)
#Gen = Generated(eta,Absorbed)
#print('Gen =', Gen * q,'. Example value is ~2-5')
#print('Gen =', Gen * q * c0**2 * hPlanck**2 * 1 / (kB**2 * Tcell**2),'. Example value is ~2-5')
#print('GenDA in photons*C/s =', Gen * q * c0**2 * hPlanck / (kB * Tcell),'. Example value is ~2-5')
#print('RR0 =', RR * q, '. Example value is ~10^-8')
#print('RR0 =', RR * q* c0**2 * hPlanck**2 * 1 / (kB**2 * Tcell**2), '. Example value is ~10^-8')
#print('RR0 and Gen need to be converted to Amps = C/s')




def Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell = 300, n = 1, Ns = 1):
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    return data['p_mp']

AbsTotal = As.round(8)
#AbsT = scipy.interpolate.interp1d(lams, AbsByAbsorbers)#, fill_value="extrapolate")
AbsTotal = scipy.interpolate.interp1d(Ephoton, AbsTotal)

Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K) 

def Qabs(eta, AbsTotal):
    def LowerB():
        return E_min
    def UpperB():
        return E_max
    def integrand(self,E):
        return eta * AbsTotal(E) * SPhotonsPerTEA(E)
    Qabs = scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]
    return Qabs

#Caluclate equilibrium Tcell

def TcellCalc(Ti,To, eta, Absorbed, AbsTotal, Tcell = 300):
    Tcell = (Qabs(eta,AbsTotal) - Give_Pmp(eta,Absorbed,Rs,Rsh) + Ui*Ti + Uo*To)/(Ui + Uo)
    return Tcell

Tcell = TcellCalc(300,300,0.6,Absorbed,AbsTotal)
print('Tcell = ',Tcell)

data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed, Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)

Isc = data['i_sc']
Voc = data['v_oc']
Imp = data['i_mp']
Vmp = data['v_mp']
Pmp = data['p_mp']
Vvalues = np.array(data['v'])
Ivalues = np.array(data['i'])
print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)

plt.figure()
plt.plot(Vvalues,Ivalues, label = 'IV')
plt.xlabel('Voltage, (V)')
plt.ylabel('Current (A)')
P_values = np.array([Ivalues * Vvalues])
plt.plot(Vvalues , P_values.T, label = 'Power')
plt.ylim(-1, 150)
plt.legend(loc = 'upper right')
plt.show()


PCE = Pmp / solar_constant
print('PCE =', PCE)

def Give_PCE():
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)  
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    Vvalues = np.array(data['v'])
    Ivalues = np.array(data['i'])
    


    
 
    
 
    
 
    
 
    
MoopsEnd
 
  ##This stuff might be important later____________

def VLTconstraint(layers):
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT - 0.5


layers 
#thicks


def PCE():
    return Pmp/solar_constant  

Layers0 = [6000,.250,.050,.5,.1,.1,.1,6000,.1,.1,.1]
#VLTc = {'type': 'ineq', 'func': VLTconstraint}


#def DoBigCalcOptimize(layers):
    
Opti = scipy.optimize.minimize(BigoptimizePCE, Layers0)#, constraints = VLTc)#, method = 'SLSQP')
print(Opti)
 
    
 
#print(PCE())


 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
    