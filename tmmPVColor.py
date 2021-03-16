# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:08:09 2021

@author: aduell
"""
#import tmm as tmm
import numpy as np
import matplotlib.pyplot as plt
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
import tmmPCECalc as pce
# This whole thing uses microns for length

'''Gives color swatches of reflected and trsnmitted light based on input of intensity of spectrum'''
def GiveColorSwatch(T, Rf):
    lamsnm = np.array(pce.lams)
    lamsnm*=1000
    spectrumT = np.vstack((lamsnm, T)).T
    spectrumRf = np.vstack((lamsnm, Rf)).T
    plots.spectrum_plot (spectrumRf, 'Rf', 'Rf_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    plots.spectrum_plot (spectrumT, 'T', 'T_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    return

'''Converts a spectrum (nm vs intenisty) into ciexyz coordinates'''
def give_xy(spectrum):
    xyz = ciexyz.xyz_from_spectrum(spectrum)
    xyz1 = colormodels.xyz_normalize (xyz)
    xyz2 = xyz1[0:2]
    #print(xyz2)
    return(xyz2)
 
'''Plots color as points on the CIE colorchart based intensity of spectra'''
def plot_xy_on_fin(T, Rf):
    lamsnm = np.array(pce.lams)
    lamsnm*=1000
    Tspectrum = np.vstack((lamsnm, T)).T
    Rfspectrum = np.vstack((lamsnm, Rf)).T
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