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
import PVColor as pvc
import tmmPCECalc as tpc

#Designing more advanced temperature calculation


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
        
            AbsLayer_spol = tmm.inc_absorp_in_each_layer(front_spol)[i]
            AbsLayer_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[i]
            AbsLayer[-1].append( (AbsLayer_spol + AbsLayer_ppol) / 2. )
        


    AbsLayer = np.array(AbsLayer)
    return AbsLayer

def qdot(LayersMaterials, Thickness):  
    layerSpectra = GiveLayerAbs(Thickness,LayersMaterials)
    #layerSpectraInterp = tpc.GiveEInterp(layerSpectra)
    x = len(Thickness)
    qList = []
    for i in range(x):
        qList.append([])
        layerSpectraInterp = tpc.GiveEInterp(layerSpectra[i])
        qList[-1].append(tpc.GiveQ(layerSpectraInterp)/Thickness[i])
    return qList

#Maybe add a funciton where all the Qs added together = solar_constant?

#solve for surface temps
T(0) = qdot*L^^2/(2*k)+(Ts1+Ts2)/2


T0 = Ts = To + qdot,0 * L / h #boundary condition
T1 = qdot,1 *L**2/(2*k)+(Ts1+Ts2)/2
T2 = qdot,2 *L**2/(2*k)+(Ts2+Ts3)/2
T3 = qdot,3 *L**2/(2*k)+(Ts3+Ts4)/2
Tf = Tsf = Ti + qdot,f * L / h

#Describes surface temperature facing in or out based on Tinf
def SurfaceTemp(layer0, Tinf, Thick, h):
    #Tease out h here from tmm data or something
    return Tinf + (qdot(layer0, Thick)*Thick)/h

#Describes temperature of Ts1 between two layers.
def InterfaceTemp(layer, Thick, Ts1, Ts2):
    #Figure out k somehwere here
    return qdot(eta, layer, Thick)*Thick**2/(2*k)+(Ts1+Ts2)/2
    
    
#Need special equation for absorber layer in terms of qdot.
    
#Might need to have a GiveLayers somewhere if i can separate the thickness portion later

#Series of Equations #check that the x-2 and array specifiying parts are correct
def GiveTSList4Solv(eta, Thickness, LayersMaterials, Ti, To, Tlist):#List of surface temp calcs for solving
    x = len(LayersMaterials)
    if x == len(Thickness):
        #Figure out h or soemthing
        TDistList = [(SurfaceTemp(LayersMaterials[0], To, Thickness[0], h) - TList[0])]
        for i in range(x-2):
            TDistList.append(InterfaceTemp(eta, LayersMaterials[i], Thickness[i], TList[i], TList[i+1]) - TList[i]))
        return TDistList
        TDistList.append(SurfaceTemp(LayersMaterials[x], Ti, Thickness[x], h) - TList[x])
    
    else:  
        raise ValueError ('layers and Thickness lengths do not match')

#Initial Guess
def GiveTListGuess(LayersMaterials):
    x = len(LayersMaterials)
    TList = [300]*x
    #for i in range(x):
    #    TList.append(300) #Generate a list of variables T1,T2,...,Ti
    return TList

#Function to solve
def FunctionForSolving(Tlist):
    return GiveTSList4Solv(eta, Thickness, LayersMaterials, Ti, To, Tlist)

#Function that does the solving for T at all surfaces
def GiveMeTheFormuoli(eta, Thickness, LayersMaterials, Ti, To):
    eta = eta
    Thickness = Thickness
    LayersMaterials = LayersMaterials
    Ti = Ti
    To = To
    TListGuess = GiveTListGuess(LayersMaterials)
    return scipy.optimize.fsolve(FunctionForSolving, TListGuess)






