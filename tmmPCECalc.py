# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:29:21 2021

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
import tmmPVColor as pvc


# This whole thing uses microns for length

degree = np.pi/180
inc_angle = 0.*degree   
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams) #um

#Rs = .02 #* ohm #series resistance
#Rsh = 10 #* ohm #shunt resistance
#eta = 0.6
#n = 1
#Ns = 1
q = 1.602176634e-19 #elementary charge C
#Ti = 300
#To = 300
#Ui = 8.3 #W/(m**2 *K) 
#Uo = 17 #W/(m**2 *K) 
c0 = 299792458 #m/s
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  



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


# Here I calculate VLT and spit it out to the screen
def VLTSpectrum(layers):
    return Stack(layers)
def VLT(layers):
    VLTstack=Stack(layers)
    return VLTstack.get_visible_light_transmission(lams,inc_angle)

#Requires input of single absorber layer with a tuple of (lb,ub)
def GiveMinMaxVLT(AbsorberType, Bounds):
    minThick = GiveLayers([Bounds[0]], [AbsorberType]) 
    maxThick = GiveLayers([Bounds[1]], [AbsorberType])
    minimum = VLT(maxThick)
    maximum = VLT(minThick)
    return {'Material':AbsorberType.__name__,'minVLT':minimum, 'maxVLT':maximum, 'minThick':Bounds[0],
            'maxThick':Bounds[1]}

def GiveMinMaxVLTFromMaterials(LayersMaterials, AbsorberLayer, Bounds):
    AbsorberType = LayersMaterials[AbsorberLayer-1]
    minThick = GiveLayers([Bounds[0]], [AbsorberType]) 
    maxThick = GiveLayers([Bounds[1]], [AbsorberType])
    minimum = VLT(maxThick)
    maximum = VLT(minThick)
    return {'Material':AbsorberType.__name__,'minVLT':minimum, 'maxVLT':maximum, 'minThick':Bounds[0],
            'maxThick':Bounds[1]}



# ******************** Here I add PCE calculation *********************#
            


worksheet = pandas.read_excel('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Data/ASTMG173.xls')#('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
#worksheet = pandas.read_excel('/Users/lwheeler/Code/pv-window-bem/Data/astmg173.xls')
downloaded_array = np.array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)

    
Ephoton = hPlanck * c0 / lams *1e6 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck


# Interpolate to get a continuous function which I will be able to do integrals on:

AM15interp = scipy.interpolate.interp1d(AM15[:,0]/1000, AM15[:,1])
#This requires nm scale 300-2500

# Here’s the plot, it looks correct:
'''
y_values = np.array([AM15interp(x) for x in lams])
plt.figure()
plt.plot(lams , y_values)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral intensity (W/m$^2$/nm)")
plt.title("Light from the sun");
plt.show()
'''



def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton *1e6  #um
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
def PowerPerTEA(Ephoton):
    return Ephoton * SPhotonsPerTEA(Ephoton)
def Solar_Constant(Ephoton):
    #PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    return scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.

solar_constant = Solar_Constant(Ephoton)


def GivelamsInterp(Parameter):
    Curve = Parameter.round(8)
    return scipy.interpolate.interp1d(lams, Curve)

def GiveEInterp(Parameter):
    Curve = Parameter.round(8)
    return scipy.interpolate.interp1d(Ephoton, Curve)


def GiveQ(Spectra, eta = 1):#Spectra must be an interpolated function
        def integrand(E):
            return eta * Spectra(E) * PowerPerTEA(E)
        return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]        
'''
#trapz calcs
def GiveQ(Spectra, eta = 1):#Spectra must be an array
        integrand = eta*Spectra*PowerPerTEA(Ephoton)
        return -np.trapz(integrand, Ephoton)     
'''

'''
def GivePhotons(Spectra, eta):#Spectra must be an interpolated function
        def integrand(E):
            return eta * Spectra(E) * SPhotonsPerTEA(E)
        return scipy.integrate.quad(integrand, E_min, E_max)[0]        
'''
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
'''
#Using trapezoidal rule for integration instaed of quad
#AbsByAbsorbers is an aray of intensities, not an interpolated function.
def RR0(eta,Absorbed,Tcell):
    AbsByAbsorbers = AbsByAbsorbers.round(8)
    integrand = eta * AbsByAbsorbers * (Ephoton)**2 / (np.exp(Ephoton / (kB * Tcell)) - 1)
    integral = scipy.integrate.trapz(integrand, Ephoton)
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

def Generated(eta,Absorbed):
    Absorbed = Absorbed.round(8)
    integrand = eta * Absorbed * SPhotonsPerTEA(Ephoton)
#    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return np.trapz(integrand, Ephoton)
'''
def Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell, n = 1, Ns = 1):
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    return data['p_mp']

#Calculate equilibrium Tcell
def TcellCalc(TotalAbs, eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh):
    AbsTotal = GiveEInterp(TotalAbs)
    #Absorbed = GiveEInterp(Absorbed)
    #eta is included in GiveQ for simplicity but should not be used for calculating Tcell. eta is set to 1
    Qabs = GiveQ(AbsTotal)
    #def Qabs(eta, AbsTotal):
    #    def LowerB():
    #        return E_min
    #    def UpperB():
    #        return E_max
    #    def integrand(self,E):
    #        return eta * AbsTotal(E) * SPhotonsPerTEA(E)
    #    return scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]        
    Temp = lambda Tcell: (Qabs - Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell) + Ui*Ti + Uo*To)/(Ui + Uo)-Tcell
    return scipy.optimize.fsolve(Temp, 300)[0]



def GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1):
    data = pvlib.pvsystem.singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed, Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)

    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    Vvalues = np.array(data['v'])
    Ivalues = np.array(data['i'])
    #print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)

    plt.figure()
    plt.plot(Vvalues,Ivalues, label = 'IV')
    plt.xlabel('Voltage, (V)')
    plt.ylabel('Current (A) or Power (W/m^2)')
    P_values = np.array([Ivalues * Vvalues])
    plt.plot(Vvalues , P_values.T, label = 'Power')
    plt.ylim(-1, 150)
    plt.legend(loc = 'upper right')
    plt.show()
    return data




def SHGC(Ts, Ti, To, Tcell, solar_constant, Ui):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    Rtot = 1/Ui #This is approximate because Ui is assumed
    #Included in GiveQ for simplicity but should not be used for calculating SHGC
    TransTotal = GiveEInterp(Ts)
    Qtrans = GiveQ(TransTotal,1)
    #def Qtrans(eta, TransTotal):
    #    def LowerB():
    #        return E_min
    #    def UpperB():
    #        return E_max
    #    def integrand(self,E):
    #        return eta * TransTotal(E) * SPhotonsPerTEA(E)
    #    return scipy.integrate.dblquad(integrand, E_min, E_max, LowerB(), UpperB())[0]
    return (Qtrans + Ui*(Tcell-Ti) - ((To-Ti)/Rtot))/solar_constant

def max_efficiency(eta,Absorbed,Tcell, solar_constant, Rs, Rsh):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    return Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell) / solar_constant

def GiveImportantInfo(WERT, LayersMaterials,eta,Ti,To,Ui,Uo,Rs,Rsh,solar_constant):
    
    NeoThickness = WERT#['x'] 
    layers = GiveLayers(NeoThickness,LayersMaterials)

    spectra = Spectra(layers ,4)
    AbsByAbsorbers = spectra['AbsByAbsorbers']
    Ts = spectra['Ts']
    Rfs = spectra['Rfs']
    Rbs = spectra['Rbs']
    As = spectra['As']
    sanities = spectra['Total']
    Absorbed = GiveEInterp(AbsByAbsorbers)
    VLTcalc = VLT(layers)
    Tcell = TcellCalc(As,eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
    #Absorbed = tpc.GiveEInterp(tpc.Spectra(tpc.GiveLayers(Thickness, LayersMaterials),4)['AbsByAbsorbers'])
    data = GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1)
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    SHGCcalc = SHGC(Ts, Ti, To, Tcell, solar_constant, Ui)
    PCE = max_efficiency(eta,Absorbed,Tcell, solar_constant, Rs, Rsh)
    #TimePCE = (end1-start1)
    #TimeOptimize = (end2 - start2)


    #Spectral Curves
    #X = np.transpose([lams,Absorbed])
    #np.savetxt('./Output/AbsByAbsorber.txt',X,delimiter=',',header="wavelength [micron], AbsByAbsorber [1]")
    #Y = np.transpose([lams,Ts,Rfs,Rbs])
    #np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")
    plt.figure()
    plt.plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
    plt.plot(lams,Ts,color='green',marker=None,label="$T$")
    plt.plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
    plt.plot(lams,As,color='black',marker=None,label="A")
    plt.plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
    plt.plot(lams,sanities,color='gold',marker=None,label="R+A+T")
    plt.plot(lams,VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    plt.xlabel('wavelength, $\mu$m')
    plt.legend(loc = 'upper right')
    plt.show()

    plt.figure()
    plt.plot(Ephoton, Ts, color='magenta',marker=None,label="$T$")
    plt.plot(Ephoton, Rfs,color='green',marker=None,label="$R_f$")
    plt.plot(Ephoton, Rbs,color='purple',marker=None,label="$R_b$")
    plt.plot(Ephoton, AbsByAbsorbers,color='black',marker=None,label="Abs")
    #plt.plot(Ephoton,tpc.VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    plt.legend(loc = 'upper right')
    plt.xlabel('Energy, J')
    plt.show()


    '''
    lamsnm = np.array(lams)
    lamsnm*=1000
    spectrumT = np.vstack((lamsnm, Ts)).T
    spectrumRf = np.vstack((lamsnm, Rfs)).T
    '''
    '''
    plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    '''
    pvc.GiveColorSwatch(Ts, Rfs)
    pvc.plot_xy_on_fin(Ts, Rfs)





    print('PCE = ',PCE,'VLT = ', VLTcalc, 'SHGC = ',SHGCcalc, 'Tcell = ',Tcell)#,'time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
    return {'PCE':PCE, 'VLT':VLTcalc, 'SHGC':SHGCcalc, 'Tcell':Tcell,'Isc':Isc, 'Voc': Voc, 'Imp': Imp, 'Vmp': Vmp,'Pmp': Pmp}
