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
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
import PVColor as pvc
import tmmPCECalc as tpc


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
def TiO2(Thickness = 0.050):
    return Layer(Thickness,'nkTiO2','c')
def FTO(Thickness = 0.250):
    return Layer(Thickness,'nkFTO','c')
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



#GlassBound = scipy.optimize.Bounds(5999,6001)
GlassBound = (5999,6001)
TiO2Bound = (0.025,.1)
FTOBound = (0.1,0.5)
MAPIBound = (.06,.260)
AZOBound = (.1,.4)
ITOBound = (.1,.4)
ITOlowEBound = (0.03,.15)
SnO2Bound = (.025,.1)
SnO2lowEBound = (.015,.06)
SnO2lowEfatBound = (0.025,.1)
SiO2Bound = (.012,.05)
NiOBound = (.025,.1)
AgBound = (.0149, .0151)
TiO2lowEBound = (.015, .070)
TiO2lowEfatBound = (.03,.12)
BleachBound = (.180, .500)
ClAlPcBound = (.150, .600)
C60Bound = (.100,.400)
IRBound = (.030, .12)
MAPBrBound = (.250,1)
EVABound = (2999,3001)
'''
Thickness = (6000,.250,0.050,0.500,0.050,0.200,3000,6000,0.030,0.015,0.030)
LayersMaterials = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
len(Thickness) == len(LayersMaterials)
'''

#Thickness = [6000,0.05,0.25]
#LayersMaterials = [Glass,FTO,TiO2]

def GiveLayers(Thickness,LayersMaterials):
    x = len(LayersMaterials)
    if x == len(Thickness):
        Layers = []
        for i in range(x):
            Layers.append(LayersMaterials[i](Thickness[i]))
        return Layers
    else:  
        raise ValueError ('layers and Thickness lengths do not match')


def GiveBounds(LayersMaterials):
    x = len(LayersMaterials)
    Bounds = []
    for i in range(x):
        Bounds.append(LayersMaterials[i].__name__ + 'Bound')
    #Bounds = [i.replace("'", "") for i in Bounds]
    return Bounds




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
#layers = [Glass(),FTO(),TiO2(),MAPBr(),NiO(),ITO(),EVA(),Glass(),TiO2lowE(),Ag(),TiO2lowE()]


layers = [Glass(),FTO(),TiO2(),MAPI(),NiO(),ITO(),EVA(),Glass(),TiO2lowE(),Ag(),TiO2lowE()]


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
def Spectra(layers, AbsorberLayer):
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


    layerchoice = AbsorberLayer 
    #layerchoice2 = 5

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

spectra = Spectra(layers,4)


AbsByAbsorbers = spectra['AbsByAbsorbers']
Ts = spectra['Ts']
Rfs = spectra['Rfs']
Rbs = spectra['Rbs']
As = spectra['As']
sanities = spectra['Total']


QAs = tpc.GiveQ(tpc.GiveEInterp(As))
QTs = tpc.GiveQ(tpc.GiveEInterp(Ts))
QRfs = tpc.GiveQ(tpc.GiveEInterp(Rfs))
QRbs = tpc.GiveQ(tpc.GiveEInterp(Rbs))
print(QAs,QTs,QRfs,QRbs, QAs+QTs+QRfs)

Moopsway





lamsnm = np.array(lams)
lamsnm*=1000
spectrumT = np.vstack((lamsnm, Ts)).T
spectrumRf = np.vstack((lamsnm, Rfs)).T
plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()
plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()
pvc.plot_xy_on_fin(spectrumT, spectrumRf)



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
    λ = hPlanck * c0 / Ephoton *1e6  #um
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
    #units photons/m^2       *    1/J       *    J*s   * m/s   /   J^2 idk
    #Units = should be photons/(s*m^2)
    #Ephoton must convert to 0.3-2.5 nm from J


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

Absorbed = GiveEInterp(AbsByAbsorbers)#GiveEInterp(Spectra(layers,4)['AbsByAbsorbers'])
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





'''
def Qabs(eta, AbsTotal):
        def LowerB():
            return E_min
        def UpperB():
            return E_max
        def integrand(self,E):
            return eta * AbsTotal(E) * SPhotonsPerTEA(E)
        return scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]        

def Qabs2(eta, AbsTotal):
        def integrand(E):
            return eta * AbsTotal(E) * SPhotonsPerTEA(E)
        return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]        

wert = Qabs(1,Absorbed)
wert2 = Qabs2(1,Absorbed)
print(wert, wert2)
Moopsbrgd
'''

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

#Constraint on VLT
def VLTconstraint(Thickness):
    layers = GiveLayers(Thickness, LayersMaterials)
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT - 0.5
VLTc = {'type': 'eq', 'func': VLTconstraint}

type(VLTconstraint)



def MediumOptimize(Thickness):
    AbsorberLayer = 4
    layerss = GiveLayers(Thickness, LayersMaterials)
    SpectraCurves = Spectra(layerss,AbsorberLayer)
    Abbsorbed = GiveEInterp(SpectraCurves['AbsByAbsorbers'])
    Tcell = TcellCalc(SpectraCurves['As'],Ti,To,eta,Abbsorbed)
    return max_efficiency(eta,Abbsorbed,Tcell)

def dotheoptimize(Thickness):
    #layerss = GiveLayers(Thicknesses, Glass,FTO,TiO2)
    func_to_minimize = lambda x : -MediumOptimize(x)
    #bnd = scipy.optimize.Bounds(.02, .1, keep_feasible=False)#If testing a single layer use this line
    #return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP',bounds = bnd )
    return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP', bounds = (GlassBound,FTOBound,TiO2Bound,MAPBrBound,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound))#, constraints = (VLTc))

                                   
#Thickness = [6000,.05,.25, 0.5]
#Thickness = [6000,.050,.25,0.8,.050, 0.2]
#Thickness = [6000,.03,.3,0.5,.069, 0.39]
#Thickness = [0.02]
#[Glass,FTO,TiO2,C60,ClAlPc,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
Thickness = (6000,.250,0.050,0.500,0.050,0.200,3000,6000,0.030,0.015,0.030)
#LayersMaterials = [FTO]
#LayersMaterials = [Glass, FTO,TiO2, MAPBr]
#LayersMaterials = [Glass, FTO,TiO2,MAPBr, NiO, ITO]
LayersMaterials = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#Bounds = GiveBounds(LayersMaterials)
print('Sim PCE for Optimization =',MediumOptimize(Thickness))
WERT = dotheoptimize(Thickness)
print(WERT)
print('VLT = ',VLTconstraint(WERT['x'])+.5)
#print(WERT['x'])

#With 6 layers, dotheoptimize took 57 minutes to complete.No constraint on VLT
#With 6 layers and tighter bounds, only took 18 minutes. No constraint on VLT.  x: array([6.00000000e+03, 5.19785253e-02, 1.88190522e-01, 6.92350335e-01,5.89070162e-02, 1.15748085e-01])
#With 6 layers after changing starting values took 23:40. No constraint on VLT x: array([6.00000000e+03, 4.99533492e-02, 2.97889624e-01, 5.04817002e-01,6.91214981e-02, 3.89139617e-01])
#With 6 layers and reducing number of calculations too 12:50. No Constraint x: array([6.00000000e+03, 5.19785253e-02, 1.88190522e-01, 6.92350335e-01,5.89070162e-02, 1.15748085e-01])
#Adding decimal places changes the result. 6 layers yadda yadda 11:35 x: array([6.00000000e+03, 4.98528198e-02, 2.61261043e-01, 5.88426467e-01,6.27028296e-02, 1.00000000e-01]
#With 6 layers after changing starting values +better calcs took 17:24. No constraint on VLT x: array([6.00000000e+03, 4.99533492e-02, 2.97889624e-01, 5.04817002e-01,6.91214981e-02, 3.89139617e-01])
#With 4 layers and reducing number of caluclations took 6 minutes. No Constraint.
#With all 11 layers took ~30 minutes.PCE = 0.04282822642972328,nfev: 157, x: array([6.00000000e+03, 1.63842049e-01, 2.50000000e-02, 6.41565714e-01, 7.06368506e-02, 1.76317816e-01, 3.00000000e+03, 6.00000000e+03,4.93834317e-02, 1.49000000e-02, 4.69851778e-02])
#Tried again w 11 layers. 21 minutes.PCE = .04282822642972328, nfev: 157, x: array([6.00000000e+03, 1.63842049e-01, 2.50000000e-02, 6.41565714e-01, 7.06368506e-02, 1.76317816e-01, 3.00000000e+03, 6.00000000e+03,4.93834317e-02, 1.49000000e-02, 4.69851778e-02])


#Tried to plot the effect of FTO vs TiO2.
#MediumOptimize doesn't work if 2 layers or less
#Need to add variable to change which layer is the absorber
#_^ This fixes issue with 2 or fewer layers
#In SHGC need to figure out what Rtot actually is

moopsbrgd

#LayersMaterials = [FTO,TiO2]
x = np.linspace(.02, 1,num =100)
y = np.linspace(.2, 1,num =100)
Thick = [x,y]
LayersMaterials = [FTO,TiO2]
xgrid, ygrid = np.meshgrid(x, y)
xy = np.stack([xgrid, ygrid])

layerss = GiveLayers(Thick, LayersMaterials)
Abbsorbed = GiveEInterp(Spectra(layerss,1)['AbsByAbsorbers'])
Tcell = TcellCalc(Spectra(layerss,1)['As'],Ti,To,eta,Abbsorbed)
function2 = max_efficiency(eta,Abbsorbed,Tcell)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(38, -30)
ax.plot_surface(xgrid,ygrid,function2, cmap='terrain')
ax.set_xlabel('FTO Thick')
ax.set_ylabel('TiO2 Thick')
ax.set_zlabel('PCE(x, y)')
plt.show()



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
