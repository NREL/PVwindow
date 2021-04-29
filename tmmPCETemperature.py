# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:40:15 2021

@author: aduell
"""
'''
_____________________Look here for notes________________________________________________
TAbsorber in the qdot equation is throwing things off.  GIvepmp returns exp overflow when attempting to fsolve
Need to figure out what TAbsorber actually is
solver is not set up correctly. idk what is wrong with it
interfacetemp has a transformation applied to it. Verify that it works
Thickness as a varibale is annoying. Need two versions in m and um. Eg thickness1/0
um goes to GiveLayerAbs and m goes to everything else

'''

import numpy as np
#import tmm
from tmm import inc_tmm, inc_absorp_in_each_layer, inf
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
DictTh={'GlassTh':GlassTh,'TiO2Th':TiO2Th,'FTOTh':FTOTh,'MAPITh':MAPITh,'AZOTh':AZOTh,'ITOTh':ITOTh,'ITOlowETh':ITOlowETh,'SnO2Th':SnO2Th,'SnO2lowETh':SnO2lowETh,'SnO2lowEfatTh':SnO2lowEfatTh,'SiO2Th':SiO2Th,'NiOTh':NiOTh,'AgTh':AgTh,'TiO2lowETh':TiO2lowETh,'TiO2lowEfatTh':TiO2lowEfatTh,'BleachTh':BleachTh,'ClAlPcTh':ClAlPcTh,'C60Th':C60Th,'IRTh':IRTh,'MAPBrTh':MAPBrTh,'EVATh':EVATh}

Materials = [Glass,FTO,TiO2,MAPBr]#,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
Thickness = np.array(tpc.GiveThicks(Materials, DictTh))
Thickness1=np.array(Thickness)
#Thickness*=1e-6
Thickness1*=1e-6 #m
Thickness0=np.array(Thickness) #um
AbsorberLayer = 4
Rs = .002 #* ohm #series resistance
Rsh = 1000 #* ohm #shunt resistance
eta = 0.6
n = 1
Ns = 1
Ti = 304 #K
To = 316 #K
Ui = 8.3 #W/(m**2 *K) 
Uo = 17 #W/(m**2 *K)
Rtot = 1/Ui
h = 100 #W/m**2
k = .01 #W/(m*K)
TList=[300,300,300,300,300,300,300,300,300,300,300,300]
#L = layer thickness
#Qinc = solar_constant
#solar_constant = sum(Qlayers)

'''I return an abosrptance spectrum for each layer
Input units for thickness must be um'''
def GiveLayerAbs(Thickness,Materials):
    layers = tpc.GiveLayers(Thickness,Materials)

    thicks = [inf]
    iorcs = ['i']
    for layer in layers:
        thicks.append(layer.d)
        iorcs.append(layer.i_or_c)
    thicks.append(inf)
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
            front_spol = inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
            front_ppol = inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)
            #back_spol = inc_tmm('s',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
            #back_ppol = inc_tmm('p',nks_bw,thicks_bw,iorcs_bw,inc_angle,lam)
        
            AbsLayer_spol = inc_absorp_in_each_layer(front_spol)[i+1]
            AbsLayer_ppol = inc_absorp_in_each_layer(front_ppol)[i+1]
            AbsLayer[-1].append( (AbsLayer_spol + AbsLayer_ppol) / 2. )
            
    AbsLayer = np.array(AbsLayer)
    return AbsLayer



#Energy statement for finding qloss
#d/dx(k*dT/dx)=qabs-qconv=qloss
#boundaries T(xL)=TL, T(xR)=TR

'''I return a list of qdots, 1 for each layer. Units are W/m^3.
The part with TAbsorber needs to be changed. TAbsorber is the temperature of the absorbing layer.
currently it is just a number but it needs to be the actual calculated temperature of the absorbing layer.
It will ultimatley come out of FunctionForSolving or GiveSurfaceTemps
layerspectra requires thickness units of um while all other uses of thickness require m
Basic calculation of qdot is Q/thickness'''
def qdot(Materials, Thickness, eta, TList, AbsorberLayer):
    #Thickness0 = np.array(Thickness)
    #Thickness0 *= 1e6
    '''Units for thickness in layerspectra must be um'''
    layerSpectra = GiveLayerAbs(Thickness0,Materials)
    #layerSpectraInterp = tpc.GiveEInterp(layerSpectra)
    #Need to change this temperature eventaully to integral of temeprature distributtion across absorber layer only.
    #Something like Integral(T(x)) from L1 to L2. Ts1 = AbsorberLayer-2 and Ts2 = absorberlayer-1
    #TAbsorber = TList[AbsorberLayer-1]
    TAbsorber=316
    x = len(Thickness)
    
    qList = []
    for i in range(x):
        if i == AbsorberLayer-1: #done --> need to fix this to subtract electrical energy from total energy absorbed by absorber to give heat loss only
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            Qabs = tpc.GiveQ(layerSpectraInterp)#/Thickness[i]
            Qconv = (tpc.Give_Pmp(eta, layerSpectraInterp, Rs, Rsh, TAbsorber, n = 1, Ns = 1))#/Thickness[i])
            qList.append((Qabs-Qconv)/Thickness[i])
            #qList[-1].append((Qabs-Qconv)/Thickness[i])
        else:
            #qList.append([])
            layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
            qList.append(tpc.GiveQ(layerSpectraInterp)/Thickness[i])
            #qList[-1].append(tpc.GiveQ(layerSpectraInterp,1)/Thickness[i])
    return qList
    #Units W/m^3


#qq = qdot(Materials, Thickness1, eta, TList, AbsorberLayer)
#print(qq)
#MOopeort


'''I add up all the Qs and return a qdot for the whole cell'''
'''
def qdotSum(Materials, Thickness, eta, TcellAbsorber, AbsorberLayer):  
    layerSpectra = GiveLayerAbs(Thickness0,Materials)
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

'''I describe surface temperature facing in or out based on ambient temp, Tinf.
Equation is eq 3.46 from fundamentals of heat and mass transfer incropera'''
def SurfaceTemp(Tinf, Thick, qdot, h):
    #Tease out h here from tmm data or something
    return Tinf + (qdot*Thick)/h

'''I describe temperature of Ts1 or Ts2 between two layers.
To solve for Ts1 set x = 0, to solve for Ts2 set x = Thick
To find edge temperature, set Ts1 or Ts2 equal to SurfaceTemp
Equation is a transformed version of eq 3.41 from fundamentals of heat and mass transfer incropera
In incropera, x ranges from -1/2*Thick to +1/2*Thick. in the transformed version here, x ranges from 0 to Thick
See Incropera for diagram'''
def InterfaceTemp(Thick, qdot, Ts1, Ts2, k, x=0):
    #x=0
    LTD = (qdot*(1/2*Thick)**2/(2*k))*(1-((x-1/2*Thick)**2/(1/2*Thick)**2)) + ((Ts2-Ts1)/2)*(x-1/2*Thick)/(1/2*Thick) + (Ts1+Ts2)/2
    #LTD = ((Ts2-Ts1)/2) + (Ts1+Ts2)/2 = Ts2
    return LTD





'''This stuff is all for verifying how GiveLayerAbs and qdot are working'''
layerSpectra = GiveLayerAbs(Thickness,Materials)
layer1 = tpc.GiveEInterp(layerSpectra[0])
layer2 = tpc.GiveEInterp(layerSpectra[1])
layer3 = tpc.GiveEInterp(layerSpectra[2])
layer4 = tpc.GiveEInterp(layerSpectra[3])

Fullspectra = tpc.Spectra(tpc.GiveLayers(Thickness, Materials),AbsorberLayer)
As = Fullspectra['As']
AbsByAbsorbers = Fullspectra['AbsByAbsorbers']
#AbsorbanceTest = GiveLayerAbs(Thickness,Materials)

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

dotq = qdot(Materials,Thickness1,eta, TList, AbsorberLayer)
print('qdots test are',dotq)
#Moopsbrgd


 
'''Series of Equations
This returns a set of equations that are used to calculate the surface and interface temperatures of the solar cell
Uses a loop to generate the equations so it can handle an arbitrary number of layers (At least 2)

This function produces a temperature distribution list (TDistList) that will give surface temps as follows-->
The first equation gives the temp of the outside surface based on the first layer (Ts1 in InterfaceTemp). 
The loop gives the surface temp on the right side of each layer (Ts2) starting with the first layer.
The final equation gives the surface temperature on the inside of the building 
The inside and outside surface temp equations are not in the loop becuase they must include the SurfaceTemp function
in place of either Ts1 or Ts2'''
def GiveTSList4Solv(eta, Thickness, Materials, qdots, Ti, To, TList,h,k):#List of surface temp calcs for solving
    x = len(Materials)
    if x == len(Thickness):
        #Figure out h or soemthing
        TDistList = [InterfaceTemp(Thickness[0], qdots[0], (SurfaceTemp(To, Thickness[0], qdots[0], h)), TList[1],k, x=0)- TList[0]]
        #TDistList.append(InterfaceTemp(Thickness[0], qdots[0], (SurfaceTemp(To, Thickness[0], qdots[0], h)), TList[1],k, x=Thickness[0])- TList[1])
        for i in range(x-2+1): #Add 1 to account for final interface. -2 to account for boundaries
            #TDistList.append(InterfaceTemp(Thickness[i+1], qdots[i+1], TList[i+1], TList[i+2],k)- TList[i+1])
            TDistList.append(InterfaceTemp(Thickness[i], qdots[i], TList[i], TList[i+1],k,x=Thickness[i])- TList[i+1])
        TDistList.append(InterfaceTemp(Thickness[-1], qdots[-1],TList[-2],(SurfaceTemp(Ti, Thickness[-1], qdots[-1], h)),k,x=Thickness[-1])- TList[-1]) #x is set equal to Thickness to give S2 instead of S1
        return TDistList
        #TDistList = [SurfaceTemp(To, Thickness1[0], qdots[0], h)-TList[0]]
        #for i in range(x-2+1): #Add 1 to account for final interface. -2 to account for boundaries
        #    TDistList.append(InterfaceTemp(Thickness1[i+1], qdots[i+1], TList[i+1], TList[i+2],k)- TList[i+1])#eta, Materials[i], Thickness[i], TList[i], TList[i+1]) - TList[i])
        #TDistList.append(SurfaceTemp(Ti, Thickness1[-1], qdots[-1], h)- TList[-1])#Materials[x], Ti, Thickness[x], h) - TList[x])
        #return TDistList
    else:  
        raise ValueError ('layers and Thickness lengths do not match')
        
       
        
#Ideas for equations that might need to be added to this system:
#1.1 Include a qdot funciton so it only needs to be performed once and the other equations work better
#1.2 qdot has to remain a function since it depends on Tcell which is the final result
#2.1 Include a funciton giving Tcell average in order to calculate qdot correctly
#2.2 Tcell average is something like the average of the sums of integrals of layer temp distributions.    
#2.3 This will require another set of loops to generate the list of integrals

'''Initial Guess. This is used to start the numerical solver. To change the starting value,
adjust the numbers in the for loop, or switch to the commented out TList=[310]*x and adjust the temperature there'''
def GiveTListGuess(Materials):
    x = len(Materials)+1
    TList=[]
    for i in range(x):
        TList.append(316-i**2)
    #TList=[310]*x
    return TList

TList = GiveTListGuess(Materials)
print(TList)
#moopbrgd
#Ts1 = SurfaceTemp(To, Thickness1[0], dotq[0], 12)#-TList[0]
Ts2 = InterfaceTemp(Thickness1[0], dotq[0], SurfaceTemp(To, Thickness1[0], dotq[0], 12), TList[1], 1)
Ts3 = InterfaceTemp(Thickness1[1], dotq[1], TList[1], TList[2], 1)
Ts4 = InterfaceTemp(Thickness1[2], dotq[2], TList[2], TList[3], 1)
Ts5 = InterfaceTemp(Thickness1[3], dotq[3], TList[3], TList[4], 1)
'''Ts6 = InterfaceTemp(Thickness1[4], dotq[4], TList[4], TList[5], 1)
Ts7 = InterfaceTemp(Thickness1[5], dotq[5], TList[5], TList[6], 1)
Ts8 = InterfaceTemp(Thickness1[6], dotq[6], TList[6], TList[7], 1)
Ts9 = InterfaceTemp(Thickness1[7], dotq[7], TList[7], TList[8], 1)
Ts10 = InterfaceTemp(Thickness1[8], dotq[8], TList[8], TList[9], 1)
Ts11 = InterfaceTemp(Thickness1[9], dotq[9], TList[9], TList[10], 1)
Ts12 = InterfaceTemp(Thickness1[10], dotq[10],SurfaceTemp(Ti, Thickness1[10], dotq[10], 12), TList[10], 1)
#Ts13 = SurfaceTemp(Ti, Thickness1[10], dotq[10], 12)#-TList[11]'''

'''
Ts1 = SurfaceTemp(To, Thickness1[0], dotq[0], 12)#-TList[0]
Ts2 = InterfaceTemp(Thickness1[1], dotq[1], TList[1], TList[2], 1)
Ts3 = InterfaceTemp(Thickness1[2], dotq[2], TList[2], TList[3], 1)
Ts4 = InterfaceTemp(Thickness1[3], dotq[3], TList[3], TList[4], 1)
Ts5 = InterfaceTemp(Thickness1[4], dotq[4], TList[4], TList[5], 1)
Ts6 = InterfaceTemp(Thickness1[5], dotq[5], TList[5], TList[6], 1)
Ts7 = InterfaceTemp(Thickness1[6], dotq[6], TList[6], TList[7], 1)
Ts8 = InterfaceTemp(Thickness1[7], dotq[7], TList[7], TList[8], 1)
Ts9 = InterfaceTemp(Thickness1[8], dotq[8], TList[8], TList[9], 1)
Ts10 = InterfaceTemp(Thickness1[9], dotq[9], TList[9], TList[10], 1)
Ts11 = InterfaceTemp(Thickness1[10], dotq[10], TList[10], TList[11], 1)
Ts12 = SurfaceTemp(Ti, Thickness1[10], dotq[10], 12)
'''
print('Test of manual GiveTSList4Solv', Ts2,Ts3,Ts4,Ts5)#,Ts6,Ts7,Ts8,Ts9,Ts10, Ts11,Ts12)#,Ts13)
#print('Test of manual GiveTSList4Solv', Ts6,Ts7,Ts8,Ts9,Ts10, Ts11,Ts12)

ListOTs = GiveTSList4Solv(eta, Thickness1, Materials, dotq, Ti, To, TList,12,1)
print('List of T calcs',ListOTs)

#MoopsNSons


'''I get solved to give surface temperatures'''
def FunctionForSolving(TList):
    qdots = qdot(Materials, Thickness1, eta, TList, AbsorberLayer)
    return GiveTSList4Solv(eta, Thickness1, Materials, qdots, Ti, To, TList,h,k)

'''I solve for Temperature at all surfaces'''
def GiveSurfaceTemps(Materials):
    TListGuess = GiveTListGuess(Materials)
    return scipy.optimize.fsolve(FunctionForSolving, TListGuess)

#___________________________Calculating Temps____________________________________________

start1 = time.time()
Result = GiveSurfaceTemps(Materials)
end1 = time.time()

print('result is',Result)
print('Time to solve in seconds = ',end1-start1,'and minutes =',(end1-start1)/60)













#Plotting

#__________________________________Plotting________________________


'''I repackage the InterfaceTemp function so that it feeds easily into the following functions'''
def LayerTempDist(Thick, x, qdot, Ts1, Ts2, k):
    #LTD = (qdot*(1/2*Thick)**2/(2*k))*(1-((x-1/2*Thick)**2/(1/2*Thick)**2)) + ((Ts2-Ts1)/2)*(x-1/2*Thick)/(1/2*Thick) + (Ts1+Ts2)/2##
    #LTD = (qdot*Thick**2/(2*k))*(1-((x-1/2*Thick)**2/(1/2*Thick)**2))# + ((Ts2-Ts1)/2)*(x-1/2*Thick)/(1/2*Thick) + (Ts1+Ts2)/2##
    #return LTD
    return InterfaceTemp(Thick,qdot,Ts1,Ts2,k,x = x)

'''I give a list of arrays for calculating and plotting temperature distribution'''
'''Each array goes from 0 to the thickness of each layer'''
def GiveZ(Thickness):
     x = len(Thickness)
     z=[]
     
     for j in range(x):
        z.append([])
        z[-1].append(np.linspace(0,Thickness[j], num = 100))
     return z

'''I Use Z to calculate the temperature distribution for each layer'''
def GiveTSList4Plot(eta, Thickness, Materials, Ti, To, TList, AbsorberLayer, h,k):#List of surface temp calcs for solving
    W = len(Materials)
    if W == len(Thickness):
        #Figure out h or soemthing
        TDistList = []
        #Thickness = np.array(Thickness)
        #Thickness *= 1e-6
        z=GiveZ(Thickness)
        qdots = qdot(Materials,Thickness,eta,TList,AbsorberLayer)
        '''
        ThickSum=0
        for Y in range(W):
            ThickSum = ThickSum + Thickness[Y]
        '''
        #for j in range(x):
        #    z.append(np.linspace(Thickness[j],Thickness[j-1],num = (Thickness[j-1]-Thickness[j])*100))
        for i in range(W): 
            #Thickness[i] = Thickness[i]/Thickness[i]
            #Distance = Distance+Thickness[i]
            TDistList.append([])
            TDistList[-1].append(LayerTempDist(Thickness[i], z[i][0], qdots[i], TList[i], TList[i+1],k))
            #Thickness[i] = Thickness[i]+Thickness[i+1]
        return TDistList
    else:  
        raise ValueError ('layers and Thickness lengths do not match')


NeoTList = Result #TList#[300,303,306,309,309,320,317,308,307,304,300,300]#
start2 = time.time()
PlotResult = GiveTSList4Plot(eta, Thickness1, Materials, Ti, To, NeoTList, AbsorberLayer, h,k)
end2 = time.time()  
z = GiveZ(Thickness)

z1 = z[0]
z2=z1[0]
z3=z2[1]
x=z

#ThickSum = Thickness[0]+Thickness[1]+Thickness[2]+Thickness[3]
#FixME =LayerTempDist(Thickness1[1], z[1][0], dotq[1], 300, 300, 1)
#print('LayerTempDIst sample',FixME)




#print('Temp Dist is',PlotResult)
print('Time to solve in seconds = ',end2-start2)
print('Time to solve in minutes = ',(end2-start2)/60)
#print('z is',z)



plt.figure()
plt.plot(z[0][0],PlotResult[0][0],color='magenta',marker=None,label="$Layer1$")
plt.plot(z[1][0]+Thickness[0],PlotResult[1][0],color='gold',marker=None,label="$Layer2$")
plt.plot(z[2][0]+Thickness[1]+Thickness[0],PlotResult[2][0],color='red',marker=None,label="$Layer3$")
plt.plot(z[3][0]+Thickness[1]+Thickness[0]+Thickness[2],PlotResult[3][0],color='blue',marker=None,label="$Layer4$")
plt.xlabel('Thickness, $\mu$m')
plt.ylabel('Temperature, K')
plt.legend(loc = 'upper right')
plt.show()



#print((dotq[1]*(1/2*Thickness1[1])**2/(2*k))*(1-((x[1][0]*(1e-6)-1/2*Thickness1[1])**2/(1/2*Thickness1[1])**2)))
#print(  (1-((x[1][0]*(1e-6)-1/2*Thickness1[1])**2/(1/2*Thickness1[1])**2)) )
#print( (dotq[1]*(1/2*Thickness1[1])**2/(2*k))  )

#print(((300-300)/2)*(x-1/2*Thickness1[1])/(1/2*Thickness1[1]))
#print((300+300)/2)

#print((dotq[1]*(1/2*Thickness[1])**2/(2*k))*(1-((x[1][0]-1/2*Thickness[1])**2/(1/2*Thickness[1])**2))+((300-300)/2)*(x[1][0]-1/2*Thickness[1])/(1/2*Thickness[1])+(300+300)/2)
#    LTD = (qdot*(1/2*Thick)**2/(2*k))*(1-((x-1/2*Thick)**2/(1/2*Thick)**2)) + ((Ts2-Ts1)/2)*(x-1/2*Thick)/(1/2*Thick) + (Ts1+Ts2)/2##

SpecialPlotting
xval = np.linspace(0,1,num=100)
fudgefactor = lambda x: -(x-0.5)**2
#PlotResult[3][0]*=fudgefactor(x-1)

plt.figure()
plt.xticks([0,1,2,3,4])
plt.plot(xval,PlotResult[0][0],color='magenta',marker=None,label="$Layer1$")
plt.plot(xval+1,PlotResult[1][0],color='black',marker=None,label="$Layer2$")
plt.plot(xval+2,PlotResult[2][0],color='red',marker=None,label="$Layer3$")
plt.plot(xval+3,PlotResult[3][0],color='blue',marker=None,label="$Layer4$")
#plt.xlim(0,4)
#plt.plot(xval+4,PlotResult[4][0],color='gold',marker=None,label="$Layer2$")
#plt.plot(xval+5,PlotResult[5][0],color='red',marker=None,label="$Layer3$")
#plt.plot(xval+6,PlotResult[6][0],color='blue',marker=None,label="$Layer4$")
#plt.plot(xval+7,PlotResult[7][0],color='red',marker=None,label="$Layer3$")
#plt.plot(xval+8,PlotResult[8][0],color='blue',marker=None,label="$Layer4$")
#plt.plot(xval+9,PlotResult[9][0],color='gold',marker=None,label="$Layer2$")
#plt.plot(xval+10,PlotResult[10][0],color='red',marker=None,label="$Layer3$")
#plt.xlabel('Thickness, $\mu$m')
plt.xlabel('Layer')
plt.ylabel('Temperature, K')
plt.legend(loc = 'upper right')
plt.show()


