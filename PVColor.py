# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:08:09 2021

@author: aduell
"""
#import tmm as tmm
import matplotlib.pyplot as plt
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
# This whole thing uses microns for length

def GiveColorSwatch(spectrumRf, spectrumT):
    plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    return

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
    return
'''
xyT = give_xy(spectrumT)
xyRf = give_xy(spectrumRf)
plot_xy_on_fin(spectrumT, spectrumRf)
'''