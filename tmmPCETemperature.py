# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:40:15 2021

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
import tmmPCECalc as tpc
from tmmPCECalc import Glass,TiO2, FTO, MAPI,AZO,ITO,ITOlowE,SnO2,SnO2lowE,SnO2lowEfat,SiO2,NiO,Ag,TiO2lowE,TiO2lowEfat,Bleach,ClAlPc,C60,IR,MAPBr,EVA
import time

#Designing more advanced temperature calculation
#Thicknesses
degree = np.pi/180
inc_angle = 0.*degree   
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams)

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

GlassTh = 6000
TiO2Th = .050
FTOTh = .250
MAPITh = .130  #.800 
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

Thickness = [GlassTh,FTOTh,TiO2Th,MAPBrTh]#,NiOTh,ITOTh,EVATh,GlassTh,TiO2lowETh,AgTh,TiO2lowETh]
#Thickness = [FTOTh,TiO2Th]#,NiOTh,ITOTh,EVATh,GlassTh,TiO2lowETh,AgTh,TiO2lowETh]
LayersMaterials = [Glass,FTO,TiO2,MAPBr]#,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#LayersMaterials = [FTO,TiO2]#,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
Boundary = [GlassBound,FTOBound,TiO2Bound,MAPBrBound]#,NiOBound,ITOBound,EVABound,GlassBound,TiO2lowEBound,AgBound,TiO2lowEBound]
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
h = 12
k = 1
TcellAbsorber = 309
TList=[300,300,300,309]
#L = layer thickness
#Qinc = solar_constant
#solar_constant = sum(Qlayers)

#x = len(Thickness)
#        for i in range(x):
def GiveLayerAbs(Thickness,LayersMaterials):
    layers = tpc.GiveLayers(Thickness,LayersMaterials)

    thicks = [tmm.inf]
    iorcs = ['i']
    for layer in layers:
        thicks.append(layer.d)
        iorcs.append(layer.i_or_c)
    thicks.append(tmm.inf)
    iorcs.append('i')
    #thicks_bw = thicks[::-1]
    #iorcs_bw = iorcs[::-1]
    
    x = len(Thickness)
    AbsLayer = []#np.zeros((x,0))
    for i in range(x):
        AbsLayer.append([])
        for lam in lams:
            
            nks = [1]
            for layer in layers:
                 nks.append(layer.nk(lam))
            nks.append(1)
            #nks_bw = nks[::-1]
            front_spol = tmm.inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
            front_ppol = tmm.inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)
            #back_spol = tmm.inc_tmm('s',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
            #back_ppol = tmm.inc_tmm('p',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
        
            AbsLayer_spol = tmm.inc_absorp_in_each_layer(front_spol)[i+1]
            AbsLayer_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[i+1]
            AbsLayer[-1].append( (AbsLayer_spol + AbsLayer_ppol) / 2. )
            
    AbsLayer = np.array(AbsLayer)
    return AbsLayer


def qdot(LayersMaterials, Thickness, eta, TList, AbsorberLayer):  
    layerSpectra = GiveLayerAbs(Thickness,LayersMaterials)
    #layerSpectraInterp = tpc.GiveEInterp(layerSpectra)
    #Need to change this temperature eventaully to integral of temeprature distributtion across absorber layer only.
    #Something like Integral(T(x)) from L1 to L2. Ts1 = AbsorberLayer-2 and Ts2 = absorberlayer-1
    TAbsorber = TList[AbsorberLayer-1]
    x = len(Thickness)
    ThickSum = 0
    for j in range(x):
        ThickSum = ThickSum + Thickness[j]
    qList = []
    for i in range(x):
        if i == AbsorberLayer-1: #done --> need to fix this to subtract electrical energy from total energy absorbed by absorber to give heat loss only
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            Qabs = tpc.GiveQ(layerSpectraInterp)#/Thickness[i]
            Qconv = tpc.Give_Pmp(eta, layerSpectraInterp, Rs, Rsh, TAbsorber, n = 1, Ns = 1)
            qList.append((Qabs-Qconv)/ThickSum)#Thickness[i])
            #qList[-1].append((Qabs-Qconv)/Thickness[i])
        else:
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            qList.append(tpc.GiveQ(layerSpectraInterp)/ThickSum)#/Thickness[i])
            #qList[-1].append(tpc.GiveQ(layerSpectraInterp,1)/Thickness[i])
    return qList

def qdotSum(LayersMaterials, Thickness, eta, TcellAbsorber, AbsorberLayer):  
    layerSpectra = GiveLayerAbs(Thickness,LayersMaterials)
    #layerSpectraInterp = tpc.GiveEInterp(layerSpectra)
    x = len(Thickness)
    QList = 0
    ThickSum = 0
    for i in range(x):
        if i == AbsorberLayer-1: #need to fix this to subtract electrical energy from total energy absorbed by absorber to give heat loss only
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            Qabs = tpc.GiveQ(layerSpectraInterp)#/Thickness[i]
            Qconv = tpc.Give_Pmp(eta, layerSpectraInterp, Rs, Rsh, TcellAbsorber, n = 1, Ns = 1)
            QList = QList+(Qabs-Qconv)
            #qList[-1].append((Qabs-Qconv)/Thickness[i])
            ThickSum = ThickSum + Thickness[i]
            
        else:
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            QList = QList+(tpc.GiveQ(layerSpectraInterp))
            #qList[-1].append(tpc.GiveQ(layerSpectraInterp,1)/Thickness[i])
            ThickSum = ThickSum + Thickness[i]
    return QList/ThickSum

'''
#Maybe add a funciton where all the Qs added together = solar_constant?
def Qdot(LayersMaterials, Thickness):  
    layerSpectra = GiveLayerAbs(Thickness,LayersMaterials)
    #layerSpectraInterp = tpc.GiveEInterp(layerSpectra)
    x = len(Thickness)
    #QList = []
    QLayerTot = 0
    for i in range(x):
        #QList.append([])
        layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
        QLayerTot = QLayerTot + (tpc.GiveQ(layerSpectraInterp))
    return QLayerTot
'''


'''
#solve for surface temps
T(0) = qdot*L^^2/(2*k)+(Ts1+Ts2)/2


T0 = Ts = To + qdot,0 * L / h #boundary condition
Ts1 = qdot,1 *L**2/(2*k)+(Ts1+Ts2)/2
Ts2 = qdot,2 *L**2/(2*k)+(Ts2+Ts3)/2
Ts3 = qdot,3 *L**2/(2*k)+(Ts3+Ts4)/2
Tf = Tsf = Ti + qdot,f * L / h
'''
#Describes surface temperature facing in or out based on Tinf
#def SurfaceTemp(Tinf, Thick, h,qdot):
#def SurfaceTemp(layer0, Tinf, Thick, h, TcellAbsorber,AbsorberLayer):
    #Tease out h here from tmm data or something
    #return Tinf + (qdot(layer0, Thick,eta,TcellAbsorber,AbsorberLayer)*Thick)/h
    #return Tinf + (qdot*Thick)/h

def SurfaceTemp(Tinf, Thick, qdot, h):
    #Tease out h here from tmm data or something
    return Tinf + (qdot*Thick)/h

#Describes temperature of Ts1 between two layers.
def InterfaceTemp(Thick, qdot, Ts1, Ts2, k):
    #Figure out k somehwere here
    #return qdot(layer, Thick,eta,Tcell,AbsorberLayer)*Thick**2/(2*k)+(Ts1+Ts2)/2
    return qdot*Thick**2/(2*k)+(Ts1+Ts2)/2
   
    
#ST = SurfaceTemp([Glass], 300, [6000],300, 5)
#print(ST)

layerSpectra = GiveLayerAbs(Thickness,LayersMaterials)
layer1 = tpc.GiveEInterp(layerSpectra[0])
layer2 = tpc.GiveEInterp(layerSpectra[1])
layer3 = tpc.GiveEInterp(layerSpectra[2])
layer4 = tpc.GiveEInterp(layerSpectra[3])
q1 = tpc.GiveQ(layer1,1)/Thickness[0]
q2 = tpc.GiveQ(layer2,1)/Thickness[1]
q3 = tpc.GiveQ(layer3,1)/Thickness[2]
q4 = tpc.GiveQ(layer4,1)/Thickness[3]
qSum = (q1+q2+q3+q4)

Q1 = tpc.GiveQ(layer1,1)
Q2 = tpc.GiveQ(layer2,1)
Q3 = tpc.GiveQ(layer3,1)
Q4 = tpc.GiveQ(layer4,1)
qSum2 = (Q1+Q2+Q3+Q4)/(Thickness[0]+Thickness[1]+Thickness[2]+Thickness[3])

qdotsum = qdotSum(LayersMaterials,Thickness,eta, 309, AbsorberLayer)
print('small q layer sum', qSum)
print('small q sum from big Q', qSum2)
print('small q sum from sum function',qdotsum)
Fullspectra = tpc.Spectra(tpc.GiveLayers(Thickness, LayersMaterials),AbsorberLayer)
As = Fullspectra['As']
AbsByAbsorbers = Fullspectra['AbsByAbsorbers']
qAs = tpc.GiveQ(tpc.GiveEInterp(As))/(GlassTh+FTOTh+TiO2Th+MAPBrTh)


AbsorbanceTest = GiveLayerAbs(Thickness,LayersMaterials)
dotq = qdot(LayersMaterials,Thickness,eta, TList, AbsorberLayer)
print('small q sum from qdot function',(dotq[0]+dotq[1]+dotq[2]+dotq[3]))
print('using qdot function',dotq)
#Qdotsum = Qdot(LayersMaterials, Thickness)
#print(Qdotsum) #sum of independent layer absorbances
print('q for full stack',qAs) #Absorbance of all layers together


'''
plt.figure()
plt.plot(tpc.Ephoton, layer1(tpc.Ephoton),color='gray', label = 'layer1')
plt.plot(tpc.Ephoton, layer2(tpc.Ephoton),color='orange', label = 'layer2')
plt.plot(tpc.Ephoton, layer3(tpc.Ephoton),color='green', label = 'layer3')
plt.plot(tpc.Ephoton, layer4(tpc.Ephoton),color='blue', label = 'layer4')
plt.plot(tpc.Ephoton, layer3(tpc.Ephoton)+layer2(tpc.Ephoton)+layer1(tpc.Ephoton)+layer4(tpc.Ephoton),color='black', label = 'sum of layers')

#plt.plot(tpc.Ephoton, layer1(tpc.Ephoton)+layer2(tpc.Ephoton), label = 'layer1+layer2')
plt.plot(tpc.Ephoton, As,color='red', label = 'Full Stack')
plt.plot(tpc.Ephoton, AbsByAbsorbers,color='red', label = 'AbsByAbsorbers')
plt.ylabel("Intensity")
plt.title("Absorbance")
plt.legend(loc = 'upper left')
plt.show()
'''

ST1 = SurfaceTemp(Ti,q1, Thickness[0], 12)
ST2 = SurfaceTemp(Ti,q2, Thickness[1], 12)
ST3 = SurfaceTemp(Ti,q3, Thickness[2], 12)
ST4 = SurfaceTemp(Ti,q4, Thickness[3], 12)

print('surface temp test calcs',ST1,ST2,ST3,ST4)
#IFT = InterfaceTemp(LayersMaterials, Thickness, Ts1, Ts2, Tcell,k)

#Moopsbrgd

    
    
    
#Need special equation for absorber layer in terms of qdot.
    
#Might need to have a GiveLayers somewhere if i can separate the thickness portion later

#Series of Equations #check that the x-2 and array specifiying parts are correct
def GiveTSList4Solv(eta, Thickness, LayersMaterials, qdots, Ti, To, TList,h,k):#List of surface temp calcs for solving
    x = len(LayersMaterials)
    if x == len(Thickness):
        #Figure out h or soemthing
        TDistList = [(SurfaceTemp(To, Thickness[0], qdots[0], h))- TList[0]]#LayersMaterials[0], To, Thickness[0], h) - TList[0])]
        for i in range(x-2+1): #Add 1 to account for final interface. -2 to account for boundaries
            TDistList.append(InterfaceTemp(Thickness[i+1], qdots[i+1], TList[i+1], TList[i+2],k)- TList[i+1])#eta, LayersMaterials[i], Thickness[i], TList[i], TList[i+1]) - TList[i])
        TDistList.append(SurfaceTemp(Ti, Thickness[-1], qdots[-1], h)- TList[-1])#LayersMaterials[x], Ti, Thickness[x], h) - TList[x])
        return TDistList
    else:  
        raise ValueError ('layers and Thickness lengths do not match')
#Ideas for equations that might need to be added to this system:
#1.1 Include a qdot funciton so it only needs to be performed once and the other equations work better
#1.2 qdot has to remain a function since it depends on Tcell which is the final result
#2.1 Include a funciton giving Tcell average in order to calculate qdot correctly
#2.2 Tcell average is something like the average of the sums of integrals of layer temp distributions.    
#2.3 This will require another set of loops to generate the list of integrals

#Initial Guess
def GiveTListGuess(LayersMaterials):
    x = len(LayersMaterials)+1
    TList = [330]*x
    #for i in range(x):
    #    TList.append(300) #Generate a list of variables T1,T2,...,Ti
    return TList

TList = GiveTListGuess(LayersMaterials)

Ts1 = SurfaceTemp(To, Thickness[0], dotq[0], 12)
Ts2 = InterfaceTemp(Thickness[1], dotq[1], TList[1], TList[2], 1)
Ts3 = InterfaceTemp(Thickness[2], dotq[2], TList[2], TList[3], 1)
Ts4 = InterfaceTemp(Thickness[3], dotq[3], TList[3], TList[4], 1)
Ts5 = SurfaceTemp(Ti, Thickness[3], dotq[3], 12)

print('Test of manual GiveTSList4Solv', Ts1,Ts2,Ts3,Ts4,Ts5)

ListOTs = GiveTSList4Solv(eta, Thickness, LayersMaterials, dotq, Ti, To, TList,12,1)
print('List of T calcs',ListOTs)

#MoopsNSons


#Function to solve
#COuld add qdot in here like the medium optimize function. qdot = ... and referenced by tht ereturned function
def FunctionForSolving(TList):
    qdots = qdot(LayersMaterials, Thickness, eta, TList, AbsorberLayer)
    return GiveTSList4Solv(eta, Thickness, LayersMaterials, qdots, Ti, To, TList,h,k)

#Function that does the solving for T at all surfaces
def GiveMeTheFormuoli(eta, Thickness, LayersMaterials, Ti, To, AbsorberLayer):
    eta = eta
    Thickness = Thickness
    LayersMaterials = LayersMaterials
    Ti = Ti
    To = To
    AbsorberLayer = AbsorberLayer
    TListGuess = GiveTListGuess(LayersMaterials)
    return scipy.optimize.fsolve(FunctionForSolving, TListGuess)

start1 = time.time()
Result = GiveMeTheFormuoli(eta, Thickness, LayersMaterials, Ti, To, AbsorberLayer)
end1 = time.time()

print('result is',Result)
print('Time to solve in seconds = ',end1-start1)
print('Time to solve in minutes = ',(end1-start1)/60)

