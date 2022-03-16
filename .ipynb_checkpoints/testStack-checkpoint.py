import numpy as np
from wpv import Layer,Stack

#import pandas as pd

# This whole thing uses microns for length

degree = np.pi/180
#inc_angle = 10.*degree
inc_angle = 0.*degree
        
#num_lams = 500
#lams = np.linspace(0.3,2.5,num=num_lams)
lamrange = [0.3,2.5]

Glass = Layer(4000,'nkLowFeGlass','i')
TiO2 = Layer(0.05,'nkTiO2','c')
FTO = Layer(0.3,'nkFTO','c')
MAPI = Layer(0.005,'nkMAPI','c')
ITO = Layer(0.4,'nkITO','c')
SnO2 = Layer(0.5,'nkSnO2','c')
NiO = Layer(0.05,'nkNiO','c')
Ag = Layer(0.01,'nkAg','c')
TiO2lowE = Layer(0.02,'nkTiO2','c')
Bleach = Layer(0.5,'nkTiO2','c')
EVA = Layer(1500,'nkEVA','i')


#MAPI.plotnk(lams)


#layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
layers = [MAPI]

stack = Stack(layers)

A_sol = stack.get_solar_weighted_absorption(lamrange,inc_angle)

print("A_sol = " + str(A_sol))

VLT = stack.get_visible_light_transmission(lamrange,inc_angle)

print("VLT = " + str(VLT))
