"""
The purpose of this script is to give a basic demonstration of PVwindow
Author: vmwheeler
"""

import numpy as np
import wpv
import matplotlib.pyplot as plt
import pandas as pd
from colour import SpectralDistribution, XYZ_to_sRGB, convert, XYZ_to_RGB, cctf_decoding, cctf_encoding
from colour.colorimetry import sd_to_XYZ, sd_to_XYZ_integration
from colour.utilities import numpy_print_options
from colour.notation import RGB_to_HEX
from colour.plotting import plot_single_colour_swatch, ColourSwatch, plot_multi_colour_swatches

Glass = wpv.Layer(4000,'nkLowFeGlass.csv','i')
FTO = wpv.Layer(0.03,'nkFTO.csv','c')
MAPI = wpv.Layer(0.6,'nkMAPI.csv','c',isPV=True)
Ag = wpv.Layer(0.01,'nkAg.csv','c')
TiO2lowE = wpv.Layer(0.002,'nkTiO2.csv','c')
EVA = wpv.Layer(1500,'nkEVA.csv','i')


#layers = [Glass,FTO,MAPI,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
#layers = [Glass,FTO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
layers = [Glass,FTO,MAPI]

stack = wpv.Stack(layers)


num_lams = 200
lams = np.linspace(0.3,2.5,num=num_lams)


Rfs = []
As = []
Ts = []

iang = 0.

pvabs = stack.get_specular_PV_abs(lams, iang)

[Rfs,As,Ts] = stack.get_specular_RAT(lams,iang)


plt.figure()
plt.plot(lams,Rfs,label=r"$R$")
plt.plot(lams,As,label=r"$A$")
plt.plot(lams,Ts,label=r"$T$")
plt.plot(lams,pvabs,label=r"$A_{PV}$")
plt.plot(lams,Rfs+As+Ts,label=r"$R+A+T$")
plt.xlabel(r"$\lambda$, micron")
plt.ylabel(r"R, A, or T")
plt.legend(loc='upper right')
#plt.show()


eta = 1 #electron-hole pair extraction efficienc
Ti = 298 #inside temperature
To = 311 #outside temperature
Ui = 8.3 #overall heat transfer coefficient of layers inside active layer
Uo = 17 #same for outside
Rs = 0 #series resistence
Rsh = 1e5 #shunt resistence
AbsorberLayer = 3 #which layer is PV absorber layer


#stuff = wpv.get_performance_characteristics_old(stack,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,iang)

stuff = wpv.get_performance_characteristics(stack,Ti,To,Ui,Uo,Rs,Rsh,iang)


print(stuff)

'''
nmlams = lams*1000
#d = {'wavelength':nmlams,'Ts':Ts}
df_spec = pd.Series(data=Ts,index=nmlams)

sd_spec = SpectralDistribution(df_spec)

df_ill = pd.Series(data=stack.Is(lams),index=nmlams)

sd_ill = SpectralDistribution(df_ill)

#with numpy_print_options(suppress=True):
#print(sd)

XYZ = sd_to_XYZ_integration(sd_spec,illuminant=sd_ill)
#https://colour.readthedocs.io/en/develop/tutorial.html?highlight=colourswatch#convert-to-display-colours
#see above for why 100 is below
#RGB = XYZ_to_sRGB(xyz_small)
sRGB = XYZ_to_sRGB(XYZ/100.)
print(XYZ)
#print(sd)
print(sRGB)

# remove gamma transfer function to get chromatricity by scaling
# see https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_.28CIE_xyY_or_CIE_XYZ_to_sRGB.29
RGB_degamma = cctf_decoding(np.clip(sRGB,0,1),'GAMMA 2.2')
print('degamma: ' + str(RGB_degamma))

RGB_scale = RGB_degamma/np.max(RGB_degamma)
RGB_chrom = cctf_encoding(RGB_scale,'GAMMA 2.2')
print('RGB_chrome: ' + str(RGB_chrom))

#RGB = np.clip(RGB/np.max(RGB),0,1)
sRGB = np.clip(sRGB,0,1)

print(sRGB)

RGB_badchrom = sRGB/np.max(sRGB)

print('***1**err?')
HEX = RGB_to_HEX(sRGB)
print('***2**err?')
print(HEX)

HEX_chrom = RGB_to_HEX(RGB_chrom)
print(HEX_chrom)

#Hexcheck = convert(RGB,'sRGB','HexaDecimal')
#print(Hexcheck)

plot_multi_colour_swatches(
    [ColourSwatch(sRGB,'sRGB'),ColourSwatch(RGB_chrom,'chromatricity'),ColourSwatch(RGB_badchrom,'~chromatricity')],
    text_kwargs={'size': 'x-large'})


hidden = stack.get_transmitted_color(lams,iang)

print(hidden)
#hexcol = convert(sd,'Spectral Distribution','Hexadecimal')
#print(hexcol)
'''