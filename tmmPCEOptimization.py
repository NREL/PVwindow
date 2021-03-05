# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:24:50 2021

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

import tmmPCECalc as tpc
from tmmPCECalc import Glass,TiO2, FTO, MAPI,AZO,ITO,ITOlowE,SnO2,SnO2lowE,SnO2lowEfat,SiO2,NiO,Ag,TiO2lowE,TiO2lowEfat,Bleach,ClAlPc,C60,IR,MAPBr,EVA
import PVColor as pvc
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once


degree = np.pi/180
inc_angle = 0.*degree   
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams) #um

Rs = .02 #* ohm #series resistance
Rsh = 10 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K)

q = 1.602176634e-19 #elementary charge C
c0 = 299792458 #m/s
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  

#Bounds
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

#Thicknesses
GlassTh = 6000
TiO2Th = .050
FTOTh = .250
MAPITh = .130
AZOTh = 0.200
ITOTh = 0.200
ITOlowETh = 0.075
SnO2Th = 0.05
SnO2lowETh = 0.030
SnO2lowEfatTh = 0.050
SiO2Th = 0.024
NiOTh = 0.050
AgTh = 0.015
TiO2lowETh = 0.030
TiO2lowEfatTh = .060
BleachTh = 0.370
ClAlPcTh = 0.300
C60Th = 0.200
IRTh = 0.060
MAPBrTh = 0.500
EVATh = 3000




worksheet = pandas.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
downloaded_array = np.array(worksheet)
AM15 = downloaded_array[1:, [0,2]]

Ephoton = hPlanck * c0 / lams *1e6 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck

solar_constant = tpc.Solar_Constant(Ephoton)
#+++++++++Start optimization parts+++++++++++++++++++++++#

#Constraint on VLT
def VLTconstraint(Thickness):
    layers = tpc.GiveLayers(Thickness, LayersMaterials)
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT - 0.5
VLTc = {'type': 'eq', 'func': VLTconstraint}



def MediumOptimize(Thickness):
    #AbsorberLayer = 4
    layerss = tpc.GiveLayers(Thickness, LayersMaterials)
    SpectraCurves = tpc.Spectra(layerss,AbsorberLayer)
    Absorbed = tpc.GiveEInterp(SpectraCurves['AbsByAbsorbers'])
    Tcell = tpc.TcellCalc(SpectraCurves['As'],Ti,To,eta,Absorbed,Ui, Uo, Rs, Rsh)
    return tpc.max_efficiency(eta,Absorbed,Tcell, solar_constant, Rs, Rsh)

def dotheoptimize(Thickness):
    #layerss = GiveLayers(Thicknesses, Glass,FTO,TiO2)
    func_to_minimize = lambda x : -MediumOptimize(x)
    #bnd = scipy.optimize.Bounds(.02, .1, keep_feasible=False)#If testing a single layer use this line
    #return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP',bounds = bnd )
    return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP', bounds = Boundary)#, constraints = (VLTc))

def TotalOptimize(eta, Thickness, LayersMaterials, Boundaries, AbsorberLayer, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 100, n = 1, Ns = 1):
    AbsorberLayer = AbsorberLayer
    Rs = Rs
    Rsh =  Rsh
    eta = eta
    n = n
    Ns = Ns
    Ti = Ti
    To = To
    Ui = Ui 
    Uo = Uo
    return dotheoptimize(Thickness)







#Thickness = [6000,.05,.25, 0.5]
#Thickness = [6000,.050,.25,0.8,.050, 0.2]
#Thickness = [6000,.03,.3,0.5,.069, 0.39]
#Thickness = [0.02]
#Thickness = (6000,.250,0.050,0.500,0.070,0.206,3000,6000,0.040,0.022,0.038)
#Thickness = (6000,.250,0.050,0.500,0.050,0.200,3000,6000,0.030,0.015,0.030)
Thickness = [GlassTh,FTOTh,TiO2Th,MAPITh,NiOTh,ITOTh,EVATh,GlassTh,TiO2lowETh,AgTh,TiO2lowETh]
#LayersMaterials = [FTO]
#LayersMaterials = [Glass, FTO,TiO2, MAPBr]
#LayersMaterials = [Glass, FTO,TiO2,MAPBr, NiO, ITO]
LayersMaterials = [Glass,FTO,TiO2,MAPI,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#Bounds = GiveBounds(LayersMaterials)
Boundary = [GlassBound,FTOBound,TiO2Bound,MAPIBound,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound]
AbsorberLayer = 4
Rs = .02 #* ohm #series resistance
Rsh = 10 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K)
Rtot = 8
print('Sim PCE for Optimization =',MediumOptimize(Thickness))
WERT = dotheoptimize(Thickness)
print(WERT)
print('VLT = ',VLTconstraint(WERT['x'])+.5)
#WERT2 = TotalOptimize(eta, Thickness, LayersMaterials, AbsorberLayer, Boundary, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 100, n = 1, Ns = 1)
#print(WERT2)



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


NeoThickness = WERT['x'] 
layers = tpc.GiveLayers(NeoThickness,LayersMaterials)

spectra = tpc.Spectra(layers ,4)
AbsByAbsorbers = spectra['AbsByAbsorbers']
Ts = spectra['Ts']
Rfs = spectra['Rfs']
Rbs = spectra['Rbs']
As = spectra['As']
sanities = spectra['Total']
Absorbed = tpc.GiveEInterp(AbsByAbsorbers)
VLT = tpc.VLT(layers)
Tcell = tpc.TcellCalc(As, Ti,To, eta, Absorbed, Ui, Uo, Rs, Rsh)
#Absorbed = tpc.GiveEInterp(tpc.Spectra(tpc.GiveLayers(Thickness, LayersMaterials),4)['AbsByAbsorbers'])
data = tpc.GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1)
SHGC = tpc.SHGC(eta, Ts, Ti, To, Rtot, Tcell, solar_constant, Ui)
PCE = tpc.max_efficiency(eta,Absorbed,Tcell, solar_constant, Rs, Rsh)
print('PCE = ',PCE,'VLT = ', VLT, 'SHGC = ',SHGC, 'Tcell = ',Tcell)



#Spectral Curves
X = np.transpose([lams,Absorbed])
#np.savetxt('./Output/AbsByAbsorber.txt',X,delimiter=',',header="wavelength [micron], AbsByAbsorber [1]")
Y = np.transpose([lams,Ts,Rfs,Rbs])
#np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")
plt.figure()
plt.plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
plt.plot(lams,Ts,color='green',marker=None,label="$T$")
plt.plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(lams,As,color='black',marker=None,label="A")
plt.plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
plt.plot(lams,sanities,color='gold',marker=None,label="R+A+T")
plt.plot(lams,tpc.VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
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


lamsnm = np.array(lams)
lamsnm*=1000
spectrumT = np.vstack((lamsnm, Ts)).T
spectrumRf = np.vstack((lamsnm, Rfs)).T
'''
plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()
plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()
'''
pvc.GiveColorSwatch(spectrumRf, spectrumT)
pvc.plot_xy_on_fin(spectrumT, spectrumRf)




moopsbrgd

#LayersMaterials = [FTO,TiO2]
x = np.linspace(.02, 1,num =100)
y = np.linspace(.2, 1,num =100)
Thick = [x,y]
LayersMaterials = [FTO,TiO2]
xgrid, ygrid = np.meshgrid(x, y)
xy = np.stack([xgrid, ygrid])

layerss = tpc.GiveLayers(Thick, LayersMaterials)
Abbsorbed = tpc.GiveEInterp(tpc.Spectra(layerss,1)['AbsByAbsorbers'])
Tcell = tpc.TcellCalc(tpc.Spectra(layerss,1)['As'],Ti,To,eta,Abbsorbed, Ui, Uo, Rs, Rsh)
function2 = tpc.max_efficiency(eta,Abbsorbed,Tcell)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(38, -30)
ax.plot_surface(xgrid,ygrid,function2, cmap='terrain')
ax.set_xlabel('FTO Thick')
ax.set_ylabel('TiO2 Thick')
ax.set_zlabel('PCE(x, y)')
plt.show()

