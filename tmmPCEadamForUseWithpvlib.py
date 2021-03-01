# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:30:47 2021

@author: aduell
"""


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



def Glass(Thickness = 6000):
    return Layer(Thickness,'nkLowFeGlass','i')
#def Glass(Thickness = 6000):
 #   return Layer(Thickness,'nkLowFeGlass','i')
def TiO2(Thickness = 0.050):
    return Layer(Thickness,'nkTiO2','c')
def FTO(Thickness = 0.250):
    return Layer(Thickness,'nkFTO','c')

#Glass = Layer(6000,'nkLowFeGlass','i')
#TiO2 = Layer(0.050,'nkTiO2','c')
#FTO = Layer(0.250,'nkFTO','c')
def MAPI(Thickness = 0.130): 
    return Layer(Thickness,'nkMAPI','c')
def AZO(Thickness = 0.200):
    return Layer(Thickness,'nkAZO','c')
def ITO(Thickness = 0.200):
    return Layer(Thickness,'nkITO','c')
def ITOlowE(Thickness = 0.075):
    return Layer(Thickness,'nkITO','c')
def SnO2(Thickness = 0.05):
    return Layer(Thickness,'nkSnO2','c')
def SnO2lowE(Thickness = 0.030):
    return Layer(Thickness,'nkSnO2','c')
def SnO2lowEfat(Thickness = 0.050):
    return Layer(Thickness,'nkSnO2','c')
def SiO2(Thickness = 0.024):
    return Layer(Thickness,'nkSiO2','c')
def NiO(Thickness = 0.050):
    return Layer(Thickness,'nkNiO','c')
def Ag(Thickness = 0.015):
    return Layer(Thickness,'nkAg','c')
def TiO2lowE(Thickness = 0.030):
    return Layer(Thickness,'nkTiO2','c')
def TiO2lowEfat(Thickness = 0.060):
    return Layer(Thickness,'nkTiO2','c')
def Bleach(Thickness = 0.370):
    return Layer(Thickness,'nkBleach','c')
def ClAlPc(Thickness = 0.300):
    return Layer(Thickness,'nkClAlPc','c')
def C60(Thickness = 0.200):
    return Layer(Thickness,'nkC60','c')
def IR(Thickness = 0.060):
    return Layer(Thickness,'nkPTB7_ThIEICO_4F','c')
def MAPBr(Thickness = 0.500):
    return Layer(Thickness,'nkMAPbBr3','c')
def EVA(Thickness = 3000):
    return Layer(Thickness,'nkEVA','i')



GlassBound = (5999,6001)
TiO2Bound = (0.025,.1)
FTOBound = (0.1,0.5)
MAPIBound = (.06,.260)
AZOBound = (.1,.4)
ITOBound = (.1,.4)
ITOlowEBound = (0.03,.15)
SnO2Bound = (.025,.1)
SnO2lowEBound = (.015,.06)
SnO2lowEfatBoun = (0.025,.1)
SiO2Bound = (.012,.05)
NiOBound = (.025,.1)
AgBound = (.007, .03)
TiO2lowEBound = (.015, .06)
TiO2lowEfatBound = (.03,.12)
BleachBound = (.180, .500)
ClAlPcBound = (.150, .600)
C60Bound = (.100,.400)
IRBound = (.030, .12)
MAPBrBound = (.250,1)
EVABound = (2999,3001)

#Thickness = [6000,0.05,0.25]
LayersMaterials = [Glass,FTO,TiO2]

def GiveLayers(Thickness,LayersMaterials):
    x = len(LayersMaterials)
    if x == len(Thickness):
        Layers = []
        for i in range(x):
            Layers.append(LayersMaterials[i](Thickness[i]))
        return Layers
    else:  
        return print('Error, Number of layers and Thickness values do not match')
#GiveLayers(Thickness, LayersMaterials)

def GiveBounds(LayersMaterials):
    x = len(LayersMaterials)
    Bounds = []
    for i in range(x):
        Bounds.append(LayersMaterials[i].__name__ + 'Bound')
    return Bounds

Bounded = GiveBounds(LayersMaterials)
print(Bounded)
print('GlassBound')

'''
def GiveLayers(Thickness,layer1,layer2 = None ,layer3 = None):
    return [layer1(Thickness[0]),layer2(Thickness[1]),layer3(Thickness[2])]
'''

#layers = GiveLayers(Thickness, Glass,FTO,TiO2)

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
layers = [Glass(),FTO(),TiO2(),MAPBr(),NiO(),ITO(),EVA(),Glass(),TiO2lowE(),Ag(),TiO2lowE()]

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

#Calculates Spectra Based on the layers of the cell
def Spectra(layers):
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
    Spectra = {'AbsByAbsorbers':AbsByAbsorbers, 'Ts':Ts,'Rfs':Rfs,'Rbs':Rbs,'As':As,'Total':sanities}
    return Spectra

spectra = Spectra(layers)



AbsByAbsorbers = spectra['AbsByAbsorbers']
Ts = spectra['Ts']
Rfs = spectra['Rfs']
Rbs = spectra['Rbs']
As = spectra['As']
sanities = spectra['Total']



# Here I calculate VLT and spit it out to the screen
def VLTSpectrum(layers):
    return Stack(layers)
def VLT(layers):
    VLTstack=Stack(layers)
    return VLTstack.get_visible_light_transmission(lams,inc_angle)
#VLTstack=Stack(layers)
#VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
print("VLT =",VLT(layers))
#

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
plt.plot(lams,VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
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
#lams *= 1000 #* nm
#++++++++++This line is breaking the calculation of absorbed later on++++++++++++#####




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
    

#λ_min = 300 #* nm
#λ_max = 2500 #* nm
#λ_min = 280 * nm
#λ_max = 4000 * nm
Ephoton = hPlanck * c0 / lams *1e6 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck
#E_min = hPlanck * c0 / λ_max * 1e9 #J   energy units from hPlanck
#E_max = hPlanck * c0 / λ_min * 1e9 #J   energy units from hPlanck



#EphotonTest = np.linspace(E_max, E_min, num=500)
#EphotonTestReverse = np.linspace(E_min, E_max, num=500)
plt.plot(Ephoton, Ts, color='magenta',marker=None,label="$T$")
plt.plot(Ephoton, Rfs,color='green',marker=None,label="$R_f$")
plt.plot(Ephoton, Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(Ephoton, AbsByAbsorbers,color='black',marker=None,label="Abs")
#plt.plot(Ephoton, DarkAeVnointerp[:,1], color='blue',marker=None,label="$Abs$")
plt.legend(loc = 'upper right')
plt.xlabel('Energy, J')
plt.show()


# Interpolate to get a continuous function which I will be able to do integrals on:

AM15interp = scipy.interpolate.interp1d(AM15[:,0]/1000, AM15[:,1])#, fill_value="extrapolate")
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
    λ = hPlanck * c0 / Ephoton *1e6  #nm
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
    #units photons/m^2       *    1/J       *    J*s   * m/s   /   J^2 idk
    #Units = should be photons/(s*m^2)
    #Ephoton must convert to 300-2500 nm from J


PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
print('Solar constant =',solar_constant) #/ (W/m**2))

# I need to get a continuous function of the fraction of the photons absorbed
# in the absorber layer. Here I tack on units and interpolate:
    
# X[:,0] *= 1000 * nm

# AbsInterp = scipy.interpolate.interp1d(X[:,0], X[:,1])

# Tack on units
#lams *= 1000 #* nm

# Round AbsByAbsorber to make really small numbers equal to zero
#AbsByAbsorbers = AbsByAbsorbers.round(8)
#AbsInterp = scipy.interpolate.interp1d(lams, AbsByAbsorbers)#, fill_value="extrapolate")
#EInterp = scipy.interpolate.interp1d(Ephoton, AbsByAbsorbers)

def GivelamsInterp(Parameter):
    Curve = Parameter.round(8)
    return scipy.interpolate.interp1d(lams, Curve)#, fill_value="extrapolate")

#def GiveEInterp(AbsByAbsorbers):
 #   AbsByAbsorbers = AbsByAbsorbers.round(8)
  #  return scipy.interpolate.interp1d(Ephoton, AbsByAbsorbers)
def GiveEInterp(Parameter):
    Curve = Parameter.round(8)
    return scipy.interpolate.interp1d(Ephoton, Curve)

#λs = np.linspace(λ_min, λ_max, num=500)
#Abs_values = np.array([AbsInterp(x) for x in λs])
Abs_values = np.array([GivelamsInterp(AbsByAbsorbers)(x) for x in lams])
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
Rs = .02 #* ohm #series resistance
Rsh = 10 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
q = 1.602176634e-19 #elementary charge C

Absorbed = GiveEInterp(Spectra(layers)['AbsByAbsorbers'])
#Absorbed = GiveAbsInterp(Spectra(layers)['AbsByAbsorbers'])

#Absorbed = GiveEInterp(spectra['AbsByAbsorbers'])
#Absorbed = AbsInterp


# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

def RR0(eta,Absorbed,Tcell):
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




def Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell, n = 1, Ns = 1):
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    return data['p_mp']


Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K) 


#AbsTotal = GiveInterpCurve(As)
#AbsTotal = As.round(8)
#AbsTotal = scipy.interpolate.interp1d(Ephoton, AbsTotal)

#Calculate equilibrium Tcell
def TcellCalc(TotalAbs, Ti,To, eta, Absorbed):
    AbsTotal = GiveEInterp(TotalAbs)
    def Qabs(eta, AbsTotal):
        def LowerB():
            return E_min
        def UpperB():
            return E_max
        def integrand(self,E):
            return eta * AbsTotal(E) * SPhotonsPerTEA(E)
        return scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]        
    Temp = lambda Tcell: (Qabs(eta,AbsTotal) - Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell) + Ui*Ti + Uo*To)/(Ui + Uo)-Tcell
    return scipy.optimize.fsolve(Temp, 300)[0]

#TrueTempMaybe = ImplicitTcellCalc(Ti,To,eta,Absorbed,AbsTotal)
#print('True temp =',TrueTempMaybe)

#def TcellCalc(Ti,To, eta, Absorbed, AbsTotal, Tcell):
#    Temp = (Qabs(eta,AbsTotal) - Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell) + Ui*Ti + Uo*To)/(Ui + Uo)
#    return Temp

Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
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
plt.ylabel('Current (A) or Power (W/m^2)')
P_values = np.array([Ivalues * Vvalues])
plt.plot(Vvalues , P_values.T, label = 'Power')
plt.ylim(-1, 150)
plt.legend(loc = 'upper right')
plt.show()





def SHGC(eta, Ts, Ti, To, Rtot, Tcell):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    TransTotal = GiveEInterp(Ts)
    def Qtrans(eta, TransTotal):
        def LowerB():
            return E_min
        def UpperB():
            return E_max
        def integrand(self,E):
            return eta * TransTotal(E) * SPhotonsPerTEA(E)
        return scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]
    return (Qtrans(eta, TransTotal) + Ui*(Tcell-Ti) - ((To-Ti)/Rtot))/solar_constant




def max_efficiency(eta,Absorbed,Tcell):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    return Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell) / solar_constant
#these are termination points in the chain of functions
print('SHGC = ',SHGC(eta, Ts, Ti, To, 8, Tcell))
print("PCE =",max_efficiency(eta,Absorbed,Tcell))


#+++++++++Start optimization parts+++++++++++++++++++++++#


def VLTconstraint(Thickness):
    layers = GiveLayers(Thickness, LayersMaterials)
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT - 0.5

VLTc = {'type': 'eq', 'func': VLTconstraint}





def MediumOptimize(Thickness):
    layerss = GiveLayers(Thickness, LayersMaterials)
    Abbsorbed = GiveEInterp(Spectra(layerss)['AbsByAbsorbers'])
    Tcell = TcellCalc(Spectra(layerss)['As'],Ti,To,eta,Abbsorbed)
    return max_efficiency(eta,Abbsorbed,Tcell)

def dotheoptimize(Thickness):
    #layerss = GiveLayers(Thicknesses, Glass,FTO,TiO2)
    func_to_minimize = lambda x : -MediumOptimize(x)
    return scipy.optimize.minimize(func_to_minimize, Thickness,bounds = ((5999,6001),(.02,.1),(.15,.5)))#,(.2,1),(.02,.07),(.1,.4)))#, constraints = (VLTc))

                                   
Thickness = [6000,.05,.25]
#Thickness = [6000,.05,.25,0.5,.050, 0.2]
LayersMaterials = [Glass, FTO,TiO2]
#LayersMaterials = [Glass, FTO,TiO2,MAPBr, NiO, ITO]
#LayersMaterials = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
Boundeds = GiveBounds(LayersMaterials)
print('Sim PCE for Optimization =',MediumOptimize(Thickness))
WERT = dotheoptimize(Thickness)
print(WERT)
#print(WERT['x'])

#With 6 layers, dotheoptimize took 57 minutes to complete.No constraint on VLT


#Tried to plot the effect of FTO vs TiO2.
#MediumOptimize doesn't work if 2 layers or less



'''
#LayersMaterials = [FTO,TiO2]
x = np.linspace(.02, 1,num =100)
y = np.linspace(.2, 1,num =100)
xgrid, ygrid = np.meshgrid(x, y)
xy = np.stack([xgrid, ygrid])
layerss = [FTO(x),TiO2(y)]
Abbsorbed = GiveEInterp(Spectra(layerss)['AbsByAbsorbers'])
TcellA = TcellCalc(Spectra(layerss)['As'],Ti,To,eta,Abbsorbed)
Function = max_efficiency(eta,Abbsorbed,TcellA)

#def function(x):
# return np.sin(x[0])*x[1]**2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(38, -30)
ax.plot_surface(xgrid, ygrid, Function, cmap='terrain')
ax.set_xlabel('FTO Thick')
ax.set_ylabel('TiO2 Thick')
ax.set_zlabel('PCE(x, y)')
plt.show()
'''


#Fix this to combine mediumoptimize and dotheoptimize with inputs thicknesses and layersmaterials
def CompleteOptimize(Thicknesses, LayersMaterials):
    LayersMaterials = LayersMaterials
    def OptimizeThis(Thicknesses):
        layerss = GiveLayers(Thicknesses, LayersMaterials) 
        Abbsorbed = GiveEInterp(Spectra(layerss)['AbsByAbsorbers'])
        return max_efficiency(eta,Abbsorbed)
    def OptimizationStep(Thicknesses):
        func_to_minimize = lambda x : -OptimizeThis(x)
        return scipy.optimize.minimize(func_to_minimize, Thicknesses)
    return OptimizationStep(Thicknesses)['x']


'''
def GiveIVdata(n,Ns,Tcell,Absorbed):
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    Vvalues = np.array(data['v'])
    Ivalues = np.array(data['i'])
    plt.figure()
    plt.plot(Vvalues,Ivalues, label = 'IV')
    plt.xlabel('Voltage, (V)')
    plt.ylabel('Current (A)')
    P_values = np.array([Ivalues * Vvalues])
    plt.plot(Vvalues , P_values.T, label = 'Power')
    plt.ylim(-1, 150)
    plt.legend(loc = 'upper right')
    plt.show()
    print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)
    return #[Isc, Voc, Imp, Vmp, Pmp]
GiveIVdata(1,1,Tcell,EInterp)
'''


'''
def Give_PCE(eta,Absorbed, Ti, To, n = 1, Ns = 1):
    #Tcell = 300
    #data1 = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)  
    #Pmp = data1['p_mp']
    #print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)  
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
    #print("PCE =", Pmp/solar_constant)
    return Pmp/solar_constant
print("PCE =", Give_PCE(0.6,EInterp,300,300))
'''
#def current_density_ideal(voltage, eta,Absorbed):
#    return e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * voltage / (kB * Tcell)))
#def current_density_ideal(voltage, eta,Absorbed):
#    return e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * voltage / (kB * Tcell)))

#def current_density(voltage, eta,Absorbed):
#    Current_dens = lambda current:(e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * (voltage + current * A * Rs) / (kB * Tcell))) - current * A)
#    Solution = scipy.optimize.fsolve(Current_dens, current_density_ideal(voltage, eta, Absorbed))
#    return Solution

#def current_density(voltage, eta,Absorbed):
#    Current_dens = lambda current:(e * ((Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * (voltage + current * A * Rs) / (kB * Tcell)) - (voltage + current * A * Rs)/Rsh)) - current * A)
#    Solution = scipy.optimize.fsolve(Current_dens, current_density_ideal(voltage, eta, Absorbed))
#    return Solution


#print('current without resistance = ', (current_density_ideal(6, 0.8, AbsInterp)))
#print('current with resistance = ', (current_density(6, 0.8, AbsInterp)))


#print(AM15[:,0])
#print(lams)
#data2 = pvlib.pvsystem.singlediode(12,RR0(eta, Absorbed),Rs,Rsh,2*kB*Tcell/q,ivcurve_pnts = 5)
#print(data2)
