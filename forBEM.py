import numpy as np
import tmm
import pandas as pd
#import tmm as tmm
import matplotlib.pyplot as plt
from wpv import Layer, Stack
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
# This whole thing uses microns for length


#defines angle of light hitting cell
degree = np.pi/180
#inc_angle = 10.*degree
inc_angle = 0.*degree


#range of wavelengths        
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams)


#defines materials with thickness and data file
Glass = Layer(6000,'nkLowFeGlass','i')
TiO2 = Layer(0.050,'nkTiO2','c')
FTO = Layer(0.250,'nkFTO','c')
MAPI = Layer(0.130,'nkMAPI','c')
AZO = Layer(0.200,'nkAZO','c')
ITO = Layer(0.200,'nkITO','c')
ITOlowE = Layer(0.075,'nkITO','c')
SnO2 = Layer(0.05,'nkSnO2','c')
SnO2lowE = Layer(0.030,'nkSnO2','c')
SnO2lowEfat = Layer(0.050,'nkSnO2','c')
SiO2 = Layer(0.024,'nkSiO2','c')
NiO = Layer(0.050,'nkNiO','c')
Ag = Layer(0.015,'nkAg','c')
TiO2lowE = Layer(0.030,'nkTiO2','c')
TiO2lowEfat = Layer(0.060,'nkTiO2','c')
Bleach = Layer(0.370,'nkBleach','c')
ClAlPc = Layer(0.300,'nkClAlPc','c')
C60 = Layer(0.200,'nkC60','c')
IR = Layer(0.060,'nkPTB7_ThIEICO_4F','c')
MAPBr = Layer(0.500,'nkMAPbBr3','c')
EVA = Layer(3000,'nkEVA','i')



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
layers = [Glass,FTO,TiO2,C60,ClAlPc,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]


# 50% VLT with non-wavelength-selective absorber, MAPbBr3 = 500 nm
#layers = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

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
EQEs = []
EQEs2 = []
IREQEs = []

#layerchoice = 4
layerchoice = 4
layerchoice2 = 5

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
    
    EQE_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice]
    EQE_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
    EQEs.append( (EQE_spol + EQE_ppol) / 2. )
    
    EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
    EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
    EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
    Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
    Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
    Ts.append( (front_spol['T']+front_ppol['T']) / 2. )



Ts = np.array(Ts)
Rfs = np.array(Rfs)
Rbs = np.array(Rbs)
As = 1-Ts-Rfs
sanities = Ts+Rfs+As

EQEs = np.array(EQEs)
EQEs2 = np.array(EQEs2)
IREQEs=EQEs+EQEs2

# Here I calculate VLT and spit it out to the screen
VLTstack=Stack(layers)
VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
print("VLT =",VLT)
#

X = np.transpose([lams,EQEs])
np.savetxt('./Output/EQE.txt',X,delimiter=',',header="wavelength [micron], EQE [1]")

Y = np.transpose([lams,Ts,Rfs,Rbs])
np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")

# This is for when there are 2 layers contributing to the EQE:
#Z = np.transpose([lams,IREQEs])
#np.savetxt('./Output/IREQE.txt',Z,delimiter=',',header="wavelength [micron], EQE [1]")


#all the different plots put together in one graph.
plt.figure()
plt.plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
plt.plot(lams,Ts,color='green',marker=None,label="$T$")
plt.plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(lams,As,color='black',marker=None,label="A")
plt.plot(lams,EQEs,color='black',linestyle='--',marker=None,label="EQE")
# This is for when there are 2 layers contributing to the EQE:
plt.plot(lams,IREQEs,color='gray',linestyle='--',marker=None,label="EQE")
plt.plot(lams,sanities,color='gold',marker=None,label="R+A+T")
# This is the photopic eye response
plt.plot(lams,VLTstack.cieplf(lams),color='red',marker=None,label="photopic")
# This is the solar spectrum
 #plt.plot(lams,VLTstack.Is(lams)/max(VLTstack.Is(lams)),color='gray',marker=None,label="AM1.5")
plt.xlabel('wavelength, $\mu$m')
plt.legend(loc = 'upper right')
plt.show()

### Added by AD


#Determine color of reflected and transmitted light
T = Y[:,1]
Rf = Y[:,2]
lamsnm = lams*1000 #convert um to nm
spectrumT = np.vstack((lamsnm, T)).T
spectrumRf = np.vstack((lamsnm, Rf)).T

plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()
plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
plt.show()

def give_xy(spectrum):
    xyz = ciexyz.xyz_from_spectrum(spectrum)
    xyz1 = colormodels.xyz_normalize (xyz)
    xyz2 = xyz1[0:2]
    #print(xyz2)
    return(xyz2)
 
def plot_xy_on_fin(Tspectrum, Rfspectrum):
    xyT = give_xy(Tspectrum)
    xyRf = give_xy(Rfspectrum)
    plots.shark_fin_plot()
    #print(xyT[0])
    plt.plot(xyT[0], xyT[1], 'ro',color = 'black', label = 'T')
    plt.plot(xyRf[0], xyRf[1], 'ro',color = 'black', label = 'T')
    plt.annotate('T', (xyT[0], xyT[1]))
    plt.annotate('Rf', (xyRf[0], xyRf[1]))
    plt.show()
 
xyT = give_xy(spectrumT)
xyRf = give_xy(spectrumRf)
plot_xy_on_fin(spectrumT, spectrumRf)