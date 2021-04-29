# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:29:21 2021

@author: aduell
"""


#import numpy as np
from numpy import pi, linspace, array, exp
#import tmm
from tmm import inc_tmm, inc_absorp_in_each_layer, inf
#import pandas as pd
#import tmm_vw as tmm
#import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,figure,xlabel,ylabel,show,ylim,legend
from wpv import Layer, Stack
#import scipy.interpolate, scipy.integrate, pandas, sys
from scipy.interpolate import interp1d
from scipy.integrate import quad, trapz
from scipy.optimize import fsolve#, Bounds
import scipy.optimize
from pandas import read_excel
import sys
#import scipy
#from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e, A, ohm
#import sympy
#import sympy.solvers.solvers
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
from pvlib.pvsystem import singlediode
import tmmPVColor as pvc
import CalculateVLTFromSpectrum as cvs
from CalculateVLTFromSpectrum import AM15G, cieplf
import vegas



# This whole thing uses microns for length

'''We determine the incident angle of the sun shining on the cell. Input is in degrees'''
def giveincangle(angle):
    degree = pi/180
    return angle*degree
inc_angle = giveincangle(0)  
'''We determine the size and scaling of the photon wavelength scale. Units are um'''   
num_lams = 500
lams = linspace(0.3,2.5,num=num_lams) #um

'''We are constants and help control units'''
q = 1.602176634e-19 #coulombs. elementary charge  
c0 = 299792458 #m/s #Speed of light
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  

'''Some units and terms'''
'''Tcell, Ti, To are cell temperature, inside temp and outside temp. Always in kelvin'''
'''Ui and Uo are overall heat-transfer coefficient ofr in side and outside. W/(m**2 *K)'''
'''AbsorberLayer is a number indicating the photoactive layer. If the fourth layer is the PV layer, input is 4'''
''''Rs is series resistance, Rsh is shunt resistance in ohms. See pveducation.org for more info''' 
'''eta is the electron-hole pair extraction efficiency term. eta times all absorbed light in the PV layer gives the EQE'''
'''n = diode ideality factor. Used in singlediode equation
Ns = number of cells in series. Used in singlediode equation'''
'''Rtot is total thermal resistance of the window'''





'''We are all the different materials currently available
Thickness is in microns'''
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


'''We are boundary conditions corresponding to each material type
Can be changed to tune optimization range'''
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


'''I assemble a list of layer objects using Thicknesses and Materials''' 
def GiveLayers(Thickness,Materials):
    x = len(Materials)
    if x == len(Thickness):
        Layers = []
        for i in range(x):
            Layers.append(Materials[i](Thickness[i]))
        return Layers
    else:  
        raise ValueError ('layers and Thickness lengths do not match')

'''I give a list of boundaries from a list of materials. Dict is a dictionary containing the boundary conditions
All items in the dicitonary are labelled as 'Material'+'Bound'  '''
'''
def GiveBounds(Materials, DictBound):
    x = len(Materials)
    Bounds = []
    for i in range(x):
        Bounds.append(DictBound[Materials[i].__name__ + 'Bound'])
    Bounds = array(Bounds)
    return Bounds
'''

'''I produce a Bounds object that defines the boundary conditions for optimization
The version above can be used to produce a list of bounds rather than an object'''

def GiveBounds(Materials, DictBound):
    x = len(Materials)
    lb = []
    ub = []
    for i in range(x):
        lb.append(DictBound[Materials[i].__name__ + 'Bound'][0])
    for i in range(x):
        ub.append(DictBound[Materials[i].__name__ + 'Bound'][1])
    bounds = scipy.optimize.Bounds(lb,ub)
    return bounds

'''I give a list of thicknesses from a list of materials. Dict is a dictionary containing the thickness values
All items in the dicitonary are labelled as 'Material'+'Th'  '''
def GiveThicks(Materials, DictTh):
    x = len(Materials)
    Th = []
    for i in range(x):
        Th.append(DictTh[Materials[i].__name__ + 'Th'])
    return Th

'''Calculates Spectra Based on the layers of the cell
AbsorberLayer is an integer giving the position of the PV layer in the stack. Currently supports 1 PV layer'''
def Spectra(layers, AbsorberLayer):
    thicks = [inf]
    iorcs = ['i']
    for layer in layers:
        thicks.append(layer.d)
        iorcs.append(layer.i_or_c)
    thicks.append(inf)
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
        
        front_spol = inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
        front_ppol = inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)
        back_spol = inc_tmm('s',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
        back_ppol = inc_tmm('p',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
    
        AbsByAbsorber_spol = inc_absorp_in_each_layer(front_spol)[layerchoice]
        AbsByAbsorber_ppol = inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
        AbsByAbsorbers.append( (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2. )
    
        # EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
        # EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
        # EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
        Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
        Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
        Ts.append( (front_spol['T']+front_ppol['T']) / 2. )


    Ts = array(Ts)
    Rfs = array(Rfs)
    Rbs = array(Rbs)
    As = 1-Ts-Rfs
    sanities = Ts+Rfs+As

    AbsByAbsorbers = array(AbsByAbsorbers)
    Spectra = {'AbsByAbsorbers':AbsByAbsorbers, 'Ts':Ts,'Rfs':Rfs,'Rbs':Rbs,'As':As,'Total':sanities}
    return Spectra




''' Here I calculate VLT and spit it out to the screen'''

'''Gives a spectrum of VLT. Used for plotting'''
def VLTSpectrum(layers):
    return Stack(layers)
'''Gives VLT as a single number'''
def VLT(layers):
    VLTstack=Stack(layers)
    return VLTstack.get_visible_light_transmission(lams,inc_angle)

'''This gives VLT as a single number. eliminates
need to recalculate AM15G and cieplf every iteration. Unclear if this will work for 
optimization'''
def getFancyVLT(layers):#,lamrange,inc_angle):
    integ = vegas.Integrator([lams])
    Trans=Stack(layers)
    numerator = integ(lambda lam: AM15G(lam)*cieplf(lam)*Trans.get_RAT(lam,inc_angle)[2], nitn=10, neval=100)[0]
    denominator = integ(lambda lam: AM15G(lam)*cieplf(lam), nitn=10, neval=100)[0]
    VLT = numerator/denominator
    return VLT.mean

'''Gives minimum and maximum VLT based exclusively on the PV layer. 
Only useful for judging VLT constraint for a given PV material
Requires input of single absorber layer with a tuple of (lb,ub)'''
def GiveMinMaxVLT(AbsorberType, Bounds):
    minThick = GiveLayers([Bounds[0]], [AbsorberType]) 
    maxThick = GiveLayers([Bounds[1]], [AbsorberType])
    minimum = VLT(maxThick)
    maximum = VLT(minThick)
    return {'Material':AbsorberType.__name__,'minVLT':minimum, 'maxVLT':maximum, 'minThick':Bounds[0],
            'maxThick':Bounds[1]}

'''Gives minimum and maximum VLT based exclusively on the PV layer. 
Requires list of materials, absorbing layer, and absorber bounds'''
def GiveMinMaxVLTFromMaterials(Materials, AbsorberLayer, Bounds):
    AbsorberType = Materials[AbsorberLayer-1]
    minThick = GiveLayers([Bounds[0]], [AbsorberType]) 
    maxThick = GiveLayers([Bounds[1]], [AbsorberType])
    minimum = VLT(maxThick)
    maximum = VLT(minThick)
    return {'Material':AbsorberType.__name__,'minVLT':minimum, 'maxVLT':maximum, 'minThick':Bounds[0],
            'maxThick':Bounds[1]}



# ******************** Here I add PCE calculation *********************#
            
'''This stuff imports a spreadsheet of the solar spectrum'''
#worksheet = pandas.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
worksheet = read_excel('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Data/ASTMG173.xls')#('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
#worksheet = pandas.read_excel('/Users/lwheeler/Code/pv-window-bem/Data/astmg173.xls')
downloaded_array = array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)

# Interpolate to get a continuous function which I will be able to do integrals on:
'''Interpolated solar spectrum
when using, inputs must be within 300-2500 nm'''
AM15interp = interp1d(AM15[:,0]/1000, AM15[:,1])


# Here’s the plot, it looks correct:
'''Plot of the solar spectrum for verification'''
'''
y_values = np.array([AM15interp(x) for x in lams])
figure()
plot(lams , y_values)
xlabel("Wavelength (nm)")
ylabel("Spectral intensity (W/m$^2$/nm)")
title("Light from the sun");
show()
'''


'''I convert wavelength to energy. E_min and max are used for integration limits '''
Ephoton = hPlanck * c0 / lams *1e6 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck


'''I give the number of photons per......'''
def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton *1e6  #um
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
'''I give the power for each......'''
def PowerPerTEA(Ephoton):
    return Ephoton * SPhotonsPerTEA(Ephoton)
'''I give the solar constant which is the W/m*2 emitted by the sun. Should be ~1000'''
def Solar_Constant(Ephoton):
    #PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    return quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
'''This is the solar constant value. It is called by optimization and used in a variety of functions here
Should always be ~1000'''
solar_constant = Solar_Constant(Ephoton)

'''I return an interpolated function of a spectrum relative to photon wavelength. Used for plotting'''
def GivelamsInterp(Parameter):
    Curve = Parameter.round(8)
    return interp1d(lams, Curve)

'''I return an interpolated function of a spectrum relative to photon energy'''
def GiveEInterp(Parameter):
    Curve = Parameter.round(8)
    return interp1d(Ephoton, Curve)

'''I give Q based on a given spectrum. Units are W/m^2
Input is a spectrum interpolated with respect to energy, E
eta should only be used if looking at a PV layer. Otherwise it is set to 1'''
def GiveQ(Spectra, eta = 1):#Spectra must be an interpolated function
        def integrand(E):
            return eta * Spectra(E) * PowerPerTEA(E)
        return quad(integrand, E_min, E_max, full_output=1)[0]       
    
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
        return quad(integrand, E_min, E_max)[0]        
'''
# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

'''I give the rate of recombination for the solar cell, Units are photons/(s*m**2)'''
def RR0(eta,Absorbed,Tcell):
    integrand = lambda E : eta * Absorbed(E) * (E)**2 / (exp(E / (kB * Tcell)) - 1)
    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return ((2 * pi) / (c0**2 * hPlanck**3)) * integral# / 1.60218e-19 #J/eV
#units = photons/(s*m**2)

'''I give the amount of energy converted to electricity in terms of photons, units are photons(s/m**2)'''
def Generated(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(E) * SPhotonsPerTEA(E)
#    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return quad(integrand, E_min, E_max, full_output=1)[0]
#units photons/(s*m**2)
'''
#Using trapezoidal rule for integration instaed of quad
#AbsByAbsorbers is an aray of intensities, not an interpolated function.
def RR0(eta,Absorbed,Tcell):
    AbsByAbsorbers = AbsByAbsorbers.round(8)
    integrand = eta * AbsByAbsorbers * (Ephoton)**2 / (np.exp(Ephoton / (kB * Tcell)) - 1)
    integral = trapz(integrand, Ephoton)
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

def Generated(eta,Absorbed):
    Absorbed = Absorbed.round(8)
    integrand = eta * Absorbed * SPhotonsPerTEA(Ephoton)
#    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return np.trapz(integrand, Ephoton)
'''

'''I use the single diode equation to return the max power of the cell in watts
Check PVlib documentation for details'''
def Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell, n = 1, Ns = 1):
    data = singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    return data['p_mp']

'''I calculate equilibrium tmperature of the cell assuming the cell is infinitely thin
TotalAbs is the full absorptance of the stack as an array of intensities, uninterpolated. 
Absorbed is PV layer absorptance interpolated
Temperature calculation is implicit so the numerical solver fsolve is used.
This equation is derived from Wheeler and Wheeler Detailed Balance Analysis of Photovoltaic Windows'''
def TcellCalc(TotalAbs, eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh):
    AbsTotal = GiveEInterp(TotalAbs)
    Qabs = GiveQ(AbsTotal)
    Temp = lambda Tcell: (Qabs - Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell) + Ui*Ti + Uo*To)/(Ui + Uo)-Tcell
    return fsolve(Temp, 300)[0]


'''I use the single diode equation to produce an IV curve and power plot
I also return related values such as Voc, Isc, and Pmp in units volts, amps, and watts
See pvlib singlediode equation for more information'''
def GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1):
    data = singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed, Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)

    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    Vvalues = array(data['v'])
    Ivalues = array(data['i'])
    #print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)

    figure()
    plot(Vvalues,Ivalues, label = 'IV')
    xlabel('Voltage, (V)')
    ylabel('Current (A) or Power (W/m^2)')
    ylabel('Power (W/m^2)')
    P_values = array([Ivalues * Vvalues])
    plot(Vvalues , P_values.T, label = 'Power')
    ylim(-1, 150)
    legend(loc = 'upper right')
    show()
    return data



'''I give the solar heat gain coefficient. unitless numebr between 0 and 1
Ts is the transmission spectra. Must be a list of intensities, not an interpolated function
This equation comes form a combination of Wheeler and Wheeler Detailed Balance Analysis of Photovoltaic Windows
and equation 3.18 from Fundamentals of Heat and Mass Transfer 6ed Incropera'''
def SHGC(Ts, Ti, To, Tcell, Ui):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    Rtot = 1/Ui #This is approximate because Ui is assumed
    #Included in GiveQ for simplicity but should not be used for calculating SHGC
    TransTotal = GiveEInterp(Ts)
    Qtrans = GiveQ(TransTotal,1)
    return (Qtrans + Ui*(Tcell-Ti) - ((To-Ti)/Rtot))/solar_constant

'''I give max efficiency also called PCE'''
'''Absorbed must be an interpolated function of the absorption spectrum of the PV layer'''
def max_efficiency(eta,Absorbed,Tcell, Rs, Rsh):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    return Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell) / solar_constant

'''I give important info about a solar cell such as PCE, SHGC, Temperature, etc'''
def GiveImportantInfo(Thickness, Materials,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,Angle=0):
    global inc_angle
    inc_angle = giveincangle(Angle)
    
    layers = GiveLayers(Thickness,Materials)
    
    spectra = Spectra(layers ,AbsorberLayer)
    AbsByAbsorbers = spectra['AbsByAbsorbers']
    Ts = spectra['Ts']
    Rfs = spectra['Rfs']
    Rbs = spectra['Rbs']
    As = spectra['As']
    sanities = spectra['Total']
    Absorbed = GiveEInterp(AbsByAbsorbers)
    VLTcalc =  cvs.getVLT(Ts,lams)#VLT(layers)
    Tcell = TcellCalc(As,eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
    #Absorbed = tpc.GiveEInterp(tpc.Spectra(tpc.GiveLayers(Thickness, Materials),4)['AbsByAbsorbers'])
    data = GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1)
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    SHGCcalc = SHGC(Ts, Ti, To, Tcell, Ui)
    PCE = max_efficiency(eta,Absorbed,Tcell, Rs, Rsh)


    #Spectral Curves
    figure()
    plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
    plot(lams,Ts,color='green',marker=None,label="$T$")
    plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
    plot(lams,As,color='black',marker=None,label="A")
    plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
    plot(lams,sanities,color='gold',marker=None,label="R+A+T")
    plot(lams,VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    xlabel('wavelength, $\mu$m')
    ylabel('Intensity')
    legend(loc = 'upper right')
    show()
    
    EphotoneV = Ephoton*6.241509e+18 
    figure()
    plot(EphotoneV, Ts, color='magenta',marker=None,label="$T$")
    plot(EphotoneV, Rfs,color='green',marker=None,label="$R_f$")
    plot(EphotoneV, Rbs,color='orange',marker=None,label="$R_b$")
    plot(EphotoneV, AbsByAbsorbers,color='black',marker=None,label="Abs")
    #plot(Ephoton,tpc.VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    legend(loc = 'upper right')
    xlabel('Energy, eV')
    ylabel('Intensity')
    show()

    pvc.GiveColorSwatch(Ts, Rfs)
    pvc.plot_xy_on_fin(Ts, Rfs)

    print('PCE = ',PCE,'VLT = ', VLTcalc, 'SHGC = ',SHGCcalc, 'Tcell = ',Tcell)#,'time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
    return {'PCE':PCE, 'VLT':VLTcalc, 'SHGC':SHGCcalc, 'Tcell':Tcell,'Isc':Isc, 'Voc': Voc, 'Imp': Imp, 'Vmp': Vmp,'Pmp': Pmp}
