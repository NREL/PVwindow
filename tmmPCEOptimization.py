# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:24:50 2021

@author: aduell
"""


#import numpy as np
#from numpy import pi, linspace, array, exp

#import tmm
#import pandas as pd
#import tmm_vw as tmm
#import matplotlib.pyplot as plt
from wpv import Layer, Stack
#import scipy.interpolate, scipy.integrate, pandas, sys
#import scipy
from scipy.optimize import minimize, differential_evolution, NonlinearConstraint, dual_annealing#, Bounds
#from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e, A, ohm
#import sympy
#import sympy.solvers.solvers
import sys
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
#import pvlib
#from pvlib import pvsystem

import tmmPCECalc as tpc
from tmmPCECalc import Glass,TiO2, FTO, MAPI,AZO,ITO,ITOlowE,SnO2,SnO2lowE,SnO2lowEfat,SiO2,NiO,Ag,TiO2lowE,TiO2lowEfat,Bleach,ClAlPc,C60,IR,MAPBr,EVA
from tmmPCECalc import VLT,GiveLayers,Spectra,GiveEInterp,TcellCalc,max_efficiency,GiveThicks,GiveBounds,GiveImportantInfo,inc_angle,lams,SHGC
#import tmmPVColor as pvc
#import tmmPCETemperature as tpt
#from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
from time import time
from statistics import stdev
import CalculateVLTFromSpectrum as cvs

#degree = np.pi/180

#inc_angle = giveincangle(0)   
#inc_angle = 0.*tpc.degree   
#num_lams = 500
#lams = np.linspace(0.3,2.5,num=num_lams) #um


'''We are boundary conditions corresponding to each material type
Can be changed to tune optimization range. All values are in um'''
GlassBound = (5999.9,6000.1)
TiO2Bound = (0.0250,.0750)#(0.025,.1)
FTOBound = (0.10,0.30)
MAPIBound = (.06,.900)#.260)
AZOBound = (.10,.30)
ITOBound = (.10,.30)
ITOlowEBound = (0.03,.09)
SnO2Bound = (.025,.075)
SnO2lowEBound = (.015,.045)
SnO2lowEfatBound = (0.025,.075)
SiO2Bound = (.012,.036)
NiOBound = (.025,.075)#(.025,.1)
AgBound = (.0130, .0170)
TiO2lowEBound = (.020, .040)
TiO2lowEfatBound = (.03,.09)
BleachBound = (.180, .500)
ClAlPcBound = (.150, .600)
C60Bound = (.100,.300)
IRBound = (.030, .12)
MAPBrBound = (.250,1)
EVABound = (2999.9,3000.1)
#PTAApolymer
'''This dictionary is used for auto-generating a list of bounds using tpc.GiveBounds'''
DictBound={'GlassBound':GlassBound,'TiO2Bound':TiO2Bound,'FTOBound':FTOBound,'MAPIBound':MAPIBound,'AZOBound':AZOBound,'ITOBound':ITOBound,'ITOlowEBound':ITOlowEBound,'SnO2Bound':SnO2Bound,'SnO2lowEBound':SnO2lowEBound,'SnO2lowEfatBound':SnO2lowEfatBound,'SiO2Bound':SiO2Bound,'NiOBound':NiOBound,'AgBound':AgBound,'TiO2lowEBound':TiO2lowEBound,'TiO2lowEfatBound':TiO2lowEfatBound,'BleachBound':BleachBound,'ClAlPcBound':ClAlPcBound,'C60Bound':C60Bound,'IRBound':IRBound,'MAPBrBound':MAPBrBound,'EVABound':EVABound}

'''We are the input thicknesses for each material. These can be changed to any desired value,
however they must be within the bounds above for the program to run correctly. All values are in um'''
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
'''This dictionary is used for auto-generating a list of thicknesses using tpc.GiveThicks'''
DictTh={'GlassTh':GlassTh,'TiO2Th':TiO2Th,'FTOTh':FTOTh,'MAPITh':MAPITh,'AZOTh':AZOTh,'ITOTh':ITOTh,'ITOlowETh':ITOlowETh,'SnO2Th':SnO2Th,'SnO2lowETh':SnO2lowETh,'SnO2lowEfatTh':SnO2lowEfatTh,'SiO2Th':SiO2Th,'NiOTh':NiOTh,'AgTh':AgTh,'TiO2lowETh':TiO2lowETh,'TiO2lowEfatTh':TiO2lowEfatTh,'BleachTh':BleachTh,'ClAlPcTh':ClAlPcTh,'C60Th':C60Th,'IRTh':IRTh,'MAPBrTh':MAPBrTh,'EVATh':EVATh}



#_________________________________Here I add Optimization Functionanilty___________________________________#

'''Constraint on VLT
To change the acceptable VLT range, change the values in VLTc. VLTc is a scipy.optimize nonlinear constraint object '''
def VLTconstraint(Thickness):
    '''V1''' 
    layers = GiveLayers(Thickness, Materials)
    #VLTc = VLT(layers)-0.5
    VLTcalc = VLT(layers)
    '''V2'''''''
    #VLTc = cvs.getVLT(Transmission,lams)-.5
    VLTc = cvs.getVLT(Transmission,lams)'''
    print('VLT=',VLTcalc)

    return VLTcalc

#VLTc = {'type': 'ineq', 'fun': VLTconstraint}
VLTc =  NonlinearConstraint(VLTconstraint, 0.5, 1.0)


'''Constraint on SHGC. Does not seem to function
 Input thickness is a placeholder to keep things consistent
with the optimization function.Currently it operates based on global variables
since the only variables that change are Tcell and Ts.
Both vary with each iteration of MediumOptimize'''
'''To change the constraint range, change the values in SHGCc'''
'''
def SHGCconstraint(Thickness):
    return #SHGC(Ts, Ti, To, Tcell, Ui)
SHGCc =  NonlinearConstraint(SHGCconstraint, 0.0, 0.5)
'''

'''This function gets optimized. Returns PCE as a function of layer thickness. Units are watts
It prints the PCE vlaue everytime it is run. Allows for monitoring of optimization progression'''
def MediumOptimize(Thickness):
    #startcalc = time()
    
    layers = GiveLayers(Thickness, Materials)
    SpectraCurves = Spectra(layers,AbsorberLayer)
    Absorbed = GiveEInterp(SpectraCurves['AbsByAbsorbers'])
    #global Transmission
    #Transmission = SpectraCurves['Ts']
    #global Tcell
    Tcell = TcellCalc(SpectraCurves['As'],eta, Ti,To,Absorbed,Ui, Uo, Rs, Rsh)
    #print('...thinking...')
    MaxEff = max_efficiency(eta,Absorbed,Tcell, Rs, Rsh)
    print(MaxEff)
    #endcalc=time()
    #print('sec/it',endcalc-startcalc)
    return  MaxEff   #  max_efficiency(eta,Absorbed,Tcell, Rs, Rsh)  




'''This function does the optimization step on MediumOptimize'''
def dotheoptimize(Thickness):
    #Using minimize function so make the equation negative to find a maximum.
    func_to_minimize = lambda x : -MediumOptimize(x)
    #bnd = scipy.optimize.Bounds(.02, .1, keep_feasible=False)#If testing a single layer use this line to designate the bounds
    return minimize(func_to_minimize, Thickness,method='SLSQP', bounds = Boundary, constraints = (VLTc),  options={'ftol': 1e-5, 'eps': 1.4901161193847656e-07,'disp': True, 'finite_diff_rel_step': None})
    #return scipy.optimize.minimize(func_to_minimize, Thickness,method='trust-constr', bounds = Boundary, constraints = (VLTc), options={'verbose':3})

'''This function is intended to be a standalone optimizaiton function. DOES NOT WORK AS INTENDED
All needed parameters are given as inputs rather than defined somewhere else in the code.
This needs to be rewritten because the variable definitions don't work as intended'''
def TotalOptimize(eta, Thickness, Materials, Boundaries, AbsorberLayer, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 1000, n = 1, Ns = 1):
    global AbsorberLayer
    AbsorberLayer = AbsorberLayer
    global Rs
    Rs = Rs
    global Rsh
    Rsh =  Rsh
    global eta
    eta = eta
    global n
    n = n
    global Ns
    Ns = Ns
    global Ti
    Ti = Ti
    global To
    To = To
    global Ui
    Ui = Ui 
    global Uo
    Uo = Uo
    return dotheoptimize(Thickness)


#___________________________Global Optimization______________________________

'''A VLT constraint for constraining the first global optimizer
Does the same thing as VLTConstraint'''
def VLTGconstraint(Thickness):
    layers = GiveLayers(Thickness, Materials)
    #VLTstack=Stack(layers)
    VLTgc=VLT(layers)
    return VLTgc
VLTGc = NonlinearConstraint(VLTGconstraint, 0.5, 1)

'''This function uses differential evolution to find a global maximum of the MediumOptimize function. Gives PCE.'''
def GlobalOptimize(Thickness):
    func_to_minimize = lambda x : -MediumOptimize(x)
    return differential_evolution(func_to_minimize, bounds = Boundary, constraints = (VLTGc)) 

'''This function uses differential annealing to to find a global maximum of MediumOptimize. Gives PCE.'''
def GlobalOptimize2(Thickness):
    func_to_minimize = lambda x : -MediumOptimize(x)
    return dual_annealing(func_to_minimize, bounds = Boundary)

#Things to look at-> boundary definitions


#_____________________________________________Operational Area__________________________________________________#

'''This stuff is used to control the optimization'''
Materials = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
Thickness = GiveThicks(Materials,DictTh)           #[GlassTh,FTOTh,TiO2Th,MAPBrTh,NiOTh,ITOTh,EVATh,GlassTh,TiO2lowETh,AgTh,TiO2lowETh]
Boundary = GiveBounds(Materials,DictBound)
#Boundary = [GlassBound,FTOBound,TiO2Bound,MAPBrBound,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound]

#Boundary = (GlassBound,FTOBound,TiO2Bound,MAPBrBound,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound)
AbsorberLayer = 4
AbsorberBoundary = MAPBrBound
Rs = .002 #* ohm #series resistance
Rsh = 1000 #* ohm #shunt resistance
eta = 0.6
n = 1 # diode ideality factor. Used in singlediode equaiton
Ns = 1 #number of cells in series. Used in singlediode equaiton
Ti = 300
To = 300
Ui = 8.3 #W/(m**2 *K) . overall heat-transfer coefficient,
Uo = 17 #W/(m**2 *K) . overall heat-transfer coefficient,
Rtot = 1/Ui

'''
layers = GiveLayers(Thickness, Materials)
startnewvlt =time()
NewVLT = tpc.getFancyVLT(layers)
endnewvlt=time()
print(NewVLT)
print(endnewvlt-startnewvlt)
'''
#Issues
#Opening up the Ag boundary cuases VLT to go super low in older version



#6 layers.time to calculate PCE from scratch in seconds =  8.253074884414673 Time to run optimizer in minutes =  10.35908164183298. nfev: 49
#8 layers.time to calculate PCE from scratch in seconds =  5.444987773895264 Time to run optimizer in minutes =  12.90131440560023. nfev: 81
#10 layers seemed to freeze up
#9 Layers also seems to stop. Suspect is TiO2LowE
#all1 layers. Actually worked this time. time to calculate PCE from scratch in seconds =  7.389192342758179 Time to run optimizer in minutes =  46.965387948354085.nfev: 234.  x: array([6.00001240e+03, 1.00000000e-01, 2.50000000e-02, 1.00000000e+00, 4.55586054e-02, 1.00629629e-01, 2.99982870e+03, 5.99995192e+03,1.50000241e-02, 1.51000000e-02, 1.50000000e-02])
#11 layes 6.42 seconds for PCE. After 2 hours still no result.
#removed solar constant as input. 11 layers 6.49 seconds for pce. 2hr10min later no result still
#8 layers 16 minutes

#6 layers 18 minutes
#6 layers after changing import style 16.49 min
#8 layers 20 min
#8 layers after changing import style 14 min
#8 layers after changing import style 12.4 min
#9 layers switching tio2 for glass 49 min
#9 294 minutes. But also failed
#9 layers. tightened tio2 bounds, 130 minutes, success,nfev: 615, nit: 51
#9 layers. Very tight TiO2 bounds, 50 min, success, nfev 196, nit 16
#10 layers Failure. 285 min, nfev 1345, nit 100
#11 layers last 3 layers MAPBr,NiO,ITO, success, nit 34, nfev 490, 102 min
#11 layers after changing constraint definition. 43 min. success. nit 14, nfev 181
#11 layers after changing constraint definition. 58 min. success. nit 23, nfev 264
#11 layers after changing Boundary definition. 36.6 min. success. nit 14, nfev 168
#" 23 and 62 minutes
#New and imporoved VLT calc and paramters. 11 layers 20 minutes. nit 15 nfev 169. ftol e-4 eps e-7
#+SHGCc, 18 min. fun =.08207, nit = 10, nfev -121, ftol e-4, eps e-6

#___________________________________Pre Optimization______________________


#This big block times the individual functions calculated during the optimization.
#Used to identify which calculations take the most time.
TcellT = 300
startlayer=time()

layerst = GiveLayers(Thickness, Materials)

endlayer = time()
startspectra= time()
spectrat = Spectra(layerst, AbsorberLayer)
endspectra= time()

startvlt =time()
VLTtest = VLT(layerst)
endvlt=time()

'''
startgiveeinterp=time()
Absby=spectra['AbsByAbsorbers']
Absorbed = GiveEInterp(Absby)
endgiveeinterp=time()
startq=time()
Q = tpc.GiveQ(Absorbed, eta = 1)
endq=time()
startrr0=time()
rr0 = tpc.RR0(eta,Absorbed,TcellT)
endrro=time()
startgenerated =time()
tpc.Generated(eta,Absorbed)
endgenerated=time()
startpmp=time()
pmp = tpc.Give_Pmp(eta, Absorbed, Rs, Rsh, TcellT, n = 1, Ns = 1)
endpmp=time()
starttcell=time()
TcellNew = TcellCalc(spectra['As'], eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
endtcell=time()
startneovlt=time()
neoVLT = cvs.getVLT(spectra['Ts'],lams)
endneovlt=time()
print('Time to calculate in sec: ','layers',endlayer-startlayer,'spectra',endspectra-startspectra,'vlt',endvlt-startvlt,'neovlt',endneovlt-startneovlt,'giveeinterp',endgiveeinterp-startgiveeinterp,'Q',endq-startq,'rr0',endrro-startrr0,'generated',endgenerated-startgenerated,'pmp',endpmp-startpmp,'tcell',endtcell-starttcell)
'''
#VLT can vary due to the way it is calculated.
#This part caluclates VLT 5 times and lets you see how different they are.
VLT1 = cvs.getVLT(spectrat['Ts'],lams)#VLT(layers)
VLT2 =  cvs.getVLT(spectrat['Ts'],lams)#VLT(layers)
VLT3 =  cvs.getVLT(spectrat['Ts'],lams)#VLT(layers)
VLT4 =  cvs.getVLT(spectrat['Ts'],lams)#VLT(layers)
VLT5 =  cvs.getVLT(spectrat['Ts'],lams)#VLT(layers)
print('preop VLT',VLTtest,VLT1,VLT2,VLT3,VLT4,VLT5)
print('Standard deviaiton of VLT is',stdev([VLTtest,VLT1,VLT2,VLT3,VLT4,VLT5]))

# = tpc.GiveMinMaxVLTFromMaterials(Materials,AbsorberLayer, AbsorberBoundary)
#print('VLT range and thicknesses of absorber are',minmax)

#This gives useful information before optimizing.Good for comparing the effects of the optimizer on PV performance
PreOpInfo = GiveImportantInfo(Thickness, Materials,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,0)
print('PreOp Calculations give',PreOpInfo)




#_________________________________Optimization happens here__________________________

#Here the optimization is calculated. It is timed and useful info is printed out.
start1 = time()
print('Sim PCE for Optimization =', MediumOptimize(Thickness))
end1 = time()
TimePCE = (end1-start1)
print('Time to calculate PCE in sec', TimePCE, 'PCE + VLT time is', TimePCE+(endvlt-startvlt))

start2 = time()
'''Optimization is actually performed here'''
WERT = dotheoptimize(Thickness)
end2 = time()


print(WERT)
#TimePCE = (end1-start1)
TimeOptimize = (end2 - start2)
print('time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
#print('VLT = ',VLTconstraint(WERT['x'])+.5)
#WERT2 = TotalOptimize(eta, Thickness, Materials, AbsorberLayer, Boundary, Ti = 300, To = 300, Ui = 8.3, Uo = 17, Rs = .02, Rsh = 100, n = 1, Ns = 1)
#print(WERT2)

WERTinfo = GiveImportantInfo(WERT['x'], Materials, eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,0)
#print('WERTinfo = ', WERTinfo)

import CalculateVLTFromSpectrum as cvs
WERTlayers = GiveLayers(WERT['x'],Materials)
WERTspectra = Spectra(WERTlayers,AbsorberLayer)
Ts = WERTspectra['Ts']
WERTVLT1 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT2 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT3 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT4 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT5 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT6 =  cvs.getVLT(Ts,lams)#VLT(layers)
WERTVLT7 =  cvs.getVLT(Ts,lams)#VLT(layers)
print(WERTVLT1,WERTVLT2,WERTVLT3,WERTVLT4,WERTVLT5,WERTVLT6,WERTVLT7)
print(sum((WERTVLT1,WERTVLT2,WERTVLT3,WERTVLT4,WERTVLT5,WERTVLT6,WERTVLT7))/7)

WERTVLT01 = VLT(WERTlayers)
WERTVLT02 = VLT(WERTlayers)
WERTVLT03 = VLT(WERTlayers)
WERTVLT04 = VLT(WERTlayers)
WERTVLT05 = VLT(WERTlayers)
WERTVLT06 = VLT(WERTlayers)
print(WERTVLT01,WERTVLT02,WERTVLT03,WERTVLT04,WERTVLT05,WERTVLT06)
print(sum((WERTVLT01,WERTVLT02,WERTVLT03,WERTVLT04,WERTVLT05,WERTVLT06))/6)


#___________________________Global Optimization happens here___________________________________________#
'''
#This stuff is for the global optimization routine using differential evolution.
#The optimization step is timed and useful info is printed out.
start3 = time()
GlobalWERT = GlobalOptimize(Thickness)
end3 = time()
GlobalTime = (end3-start3)
print(GlobalWERT)
print('Time to optimize globally in minutes = ',GlobalTime/60)
#Time to optimize globally in minutes =  161.08489356835682,nfev: 1452,fun: -0.083127734,x: array([5.99999484e+03, 1.74480480e-01, 2.50000000e-02, 9.29653260e-01,2.84263201e-02, 2.03697117e-01, 2.99967137e+03, 5.99991107e+03,7.00000000e-02, 1.50100000e-02, 7.00000000e-02])
GlobalWERTinfo = GiveImportantInfo(GlobalWERT['x'], Materials)
print(GlobalWERTinfo)
'''
'''
#This runs forever. Do not use--> Messed up VLT range.
#If VLT Range cannot be reached with absorber type then it will not solve.
start4 = time()
GlobalWERT2 = GlobalOptimize2(Thickness)
end4 = time()
GlobalTime2 = (end4-start4)
print(GlobalWERT2)
print('Time to optimize globally in minutes = ',GlobalTime2/60)
GlobalWERTinfo2 = tpc.GiveImportantInfo(GlobalWERT2['x'], Materials)
print('GlobalWERTinfo2 = ', GlobalWERTinfo2)
'''






'''
#calculated using regular minimization
BERT = {'x':[6.00000000e+03, 1.00000000e-01, 8.87880248e-02, 7.81564055e-01]}
BERTinfo = tpc.GiveImportantInfo(BERT, Materials,eta,Ti,To,Ui,Uo,Rs,Rsh)
#calculated using global optimizer differential evolution
YERT = {'x':[5.99992328e+03, 1.00000000e-01, 3.44122617e-02, 9.00000000e-01]}
YERTinfo = tpc.GiveImportantInfo(YERT, Materials,eta,Ti,To,Ui,Uo,Rs,Rsh)
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

#Materials = [FTO,TiO2]
x = np.linspace(.02, 1,num =100)
y = np.linspace(.2, 1,num =100)
Thick = [x,y]
Materials = [FTO,TiO2]
xgrid, ygrid = np.meshgrid(x, y)
xy = np.stack([xgrid, ygrid])

layerss = tpc.GiveLayers(Thick, Materials)
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


#Optimized   [6000, 0.175, 0.025, 1,0.0353, 0.1, 3000, 6000,0.04, 0.013, 0.02]
#Unoptimized [6000, 0.25, 0.05, 0.5, 0.05, 0.2, 3000, 6000, 0.03, 0.015, 0.03]



#These are reports and other things saved for later
#PCE =  0.08324261113434049 VLT =  0.4560079972898518 SHGC =  0.5032511317369591 Tcell =  316.0491481068552