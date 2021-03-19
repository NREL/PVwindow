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
import tmmPVColor as pvc
#import tmmPCETemperature as tpt
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
import time


degree = np.pi/180
inc_angle = 0.*degree   
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams) #um
'''
Rs = .02 #* ohm #series resistance
Rsh = 10 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K)
'''
q = 1.602176634e-19 #elementary charge C
c0 = 299792458 #m/s
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  

#Bounds
GlassBound = (5999.9,6000.1)
TiO2Bound = (0.025,.1)
FTOBound = (0.1,0.5)
MAPIBound = (.06,.900)#.260)
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
#PTAApolymer

#Thicknesses
GlassTh = 6000
TiO2Th = 0.050
FTOTh = 0.250
MAPITh = 0.130  #.800 
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

Ephoton = hPlanck * c0 / lams *1e6 #J tpc.Ephoton#
#E_min = min(Ephoton) #J   energy units from hPlanck
#E_max = max(Ephoton) #J   energy units from hPlanck

solar_constant = tpc.Solar_Constant(Ephoton)#tp.solar_constant
#+++++++++Start optimization parts+++++++++++++++++++++++#

#Constraint on VLT
def VLTconstraint(Thickness):
    layers = tpc.GiveLayers(Thickness, LayersMaterials)
    VLT = tpc.VLT(layers)
    return VLT - 0.5
VLTc = {'type': 'ineq', 'fun': VLTconstraint}



def MediumOptimize(Thickness):
    #AbsorberLayer = 4
    layers = tpc.GiveLayers(Thickness, LayersMaterials)
    SpectraCurves = tpc.Spectra(layers,AbsorberLayer)
    Absorbed = tpc.GiveEInterp(SpectraCurves['AbsByAbsorbers'])
    Tcell = tpc.TcellCalc(SpectraCurves['As'],eta, Ti,To,Absorbed,Ui, Uo, Rs, Rsh)
    return tpc.max_efficiency(eta,Absorbed,Tcell, solar_constant, Rs, Rsh)

def dotheoptimize(Thickness):
    #layerss = GiveLayers(Thicknesses, Glass,FTO,TiO2)
    func_to_minimize = lambda x : -MediumOptimize(x)
    #bnd = scipy.optimize.Bounds(.02, .1, keep_feasible=False)#If testing a single layer use this line
    #return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP',bounds = bnd )
    return scipy.optimize.minimize(func_to_minimize, Thickness,method='SLSQP', bounds = Boundary, constraints = (VLTc))#, options={'ftol': 1e-06, 'eps': 1.4901161193847656e-08})#, 'finite_diff_rel_step': None})

def TotalOptimize(eta, Thickness, LayersMaterials, Boundaries, AbsorberLayer, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 1000, n = 1, Ns = 1):
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


#Global Min
def VLTGconstraint(Thickness):
    layers = tpc.GiveLayers(Thickness, LayersMaterials)
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
    return VLT
VLTGc = scipy.optimize.NonlinearConstraint(VLTGconstraint, 0.5, 1)

def GlobalOptimize(Thickness):
    func_to_minimize = lambda x : -MediumOptimize(x)
    return scipy.optimize.differential_evolution(func_to_minimize, bounds = Boundary)#, constraints = (VLTGc)) #Thickness,method='SLSQP', bounds = Boundary, constraints = (VLTc), options={'ftol': 1e-06, 'eps': 1.4901161193847656e-08})#, 'finite_diff_rel_step': None})
def GlobalOptimize2(Thickness):
    func_to_minimize = lambda x : -MediumOptimize(x)
    return scipy.optimize.dual_annealing(func_to_minimize, bounds = Boundary)



Thickness = [GlassTh,FTOTh,TiO2Th,MAPBrTh,NiOTh,ITOTh,EVATh,GlassTh,TiO2lowETh,AgTh,TiO2lowETh]
LayersMaterials = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#Bounds = GiveBounds(LayersMaterials)
Boundary = [GlassBound,FTOBound,TiO2Bound,MAPBrBound,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound]
AbsorberLayer = 4
AbsorberBoundary = MAPBrBound
Rs = .002 #* ohm #series resistance
Rsh = 1000 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K)
Rtot = 1/Ui
'''
layers = tpc.GiveLayers(Thickness,LayersMaterials)
spectres = tpc.Spectra(layers,AbsorberLayer)
Absorbed = tpc.GiveEInterp(spectres['AbsByAbsorbers'])
'''
'''
RR0 = tpc.RR0(eta,Absorbed,300)
gen = tpc.Generated(eta,Absorbed)
Q = tpc.GiveQ(Absorbed)
print('RR0 = ',RR0)
print('gen = ',gen)
print('Q = ',Q)
'''


minmax = tpc.GiveMinMaxVLTFromMaterials(LayersMaterials,AbsorberLayer, AbsorberBoundary)
print('VLT range and thicknesses of absorber are',minmax)

#___________________________________Pre Optimization______________________
'''
layers=tpc.GiveLayers(Thickness,LayersMaterials)
Spectres = tpc.Spectra(layers, AbsorberLayer)
As = Spectres['As']
Absorbed = tpc.GiveEInterp(Spectres['AbsByAbsorbers'])
start0 = time.time()
TcellPreOp = tpc.TcellCalc(As, eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
end0 = time.time()
print('Tcell PreOp = ',TcellPreOp)
print('TIme to calculate Tcell in sec = ', end0-start0)
'''
PreOpInfo = tpc.GiveImportantInfo(Thickness, LayersMaterials,eta,Ti,To,Ui,Uo,Rs,Rsh,solar_constant)
#print('PreOp Calculations give',PreOpInfo)





#_________________________________Optimization here__________________________
#tpc.Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell)
start1 = time.time()
print('Sim PCE for Optimization =', MediumOptimize(Thickness))
end1 = time.time()
TimePCE = (end1-start1)
print('Time to calculate PCE in sec', TimePCE)
start2 = time.time()
WERT = dotheoptimize(Thickness)
end2 = time.time()
print(WERT)
#TimePCE = (end1-start1)
TimeOptimize = (end2 - start2)
print('time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
#print('VLT = ',VLTconstraint(WERT['x'])+.5)
#WERT2 = TotalOptimize(eta, Thickness, LayersMaterials, AbsorberLayer, Boundary, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 100, n = 1, Ns = 1)
#print(WERT2)
'''
start3 = time.time()
GlobalWERT = GlobalOptimize(Thickness)
end3 = time.time()
GlobalTime = (end3-start3)
print(GlobalWERT)
print('Time to optimize globally in minutes = ',GlobalTime/60)

#This runs forever. Do not use--> Messed up VLT range.
#If VLT Range cannot be reached with absorber type then it will not solve.
start4 = time.time()
GlobalWERT2 = GlobalOptimize2(Thickness)
end4 = time.time()
GlobalTime2 = (end4-start4)
print(GlobalWERT2)
print('Time to optimize globally in minutes = ',GlobalTime2/60)
'''




WERTinfo = tpc.GiveImportantInfo(WERT['x'], LayersMaterials, eta,Ti,To,Ui,Uo,Rs,Rsh,solar_constant)
print('WERTinfo = ', WERTinfo)
#GlobalWERTinfo = tpc.GiveImportantInfo(GlobalWERT['x'], LayersMaterials)
#print('GlobalWERTinfo = ', GlobalWERTinfo)
#GlobalWERTinfo2 = tpc.GiveImportantInfo(GlobalWERT2['x'], LayersMaterials)
#print('GlobalWERTinfo2 = ', GlobalWERTinfo2)


'''
#calculated using regular minimization
BERT = {'x':[6.00000000e+03, 1.00000000e-01, 8.87880248e-02, 7.81564055e-01]}
BERTinfo = tpc.GiveImportantInfo(BERT, LayersMaterials,eta,Ti,To,Ui,Uo,Rs,Rsh,solar_constant)
#calculated using global optimizer differential evolution
YERT = {'x':[5.99992328e+03, 1.00000000e-01, 3.44122617e-02, 9.00000000e-01]}
YERTinfo = tpc.GiveImportantInfo(YERT, LayersMaterials,eta,Ti,To,Ui,Uo,Rs,Rsh,solar_constant)
print('BERTinfo = ', BERTinfo)
print('YERTinfo = ', YERTinfo)
'''


'''
#GLoblawert data
PCE =  0.08320342823723645 VLT =  0.5240058022060908 SHGC =  0.9394718110561349 Tcell =  312.05678678191924
[6.00004059e+03, 1.75730672e-01, 2.50000000e-02, 9.97337406e-01])
    #regular
PCE =  0.0823026767378474 VLT =  0.5368890922289812 SHGC =  1.023465230604362 Tcell =  309.82218041481286
[6.00000000e+03, 1.00000000e-01, 2.50000000e-02, 6.97669968e-01])
'''

#moopsbrgd

'''
#Tried to plot the effect of FTO vs TiO2.

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
'''
