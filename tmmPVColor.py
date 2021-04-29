# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:08:09 2021

@author: aduell
"""
from numpy import linspace, array,vstack
import matplotlib.pyplot as plt
from colorpy import plots, ciexyz, colormodels #need to install colorpy to call all packages at once
import tmmPCECalc as tpc
# This whole thing uses microns for length

'''This code is all based on the ColorPy package. Check the online documentation
for more info and more functions if needed'''


'''Gives color swatches of reflected and transmitted light based on input of intensity of spectrum'''
def GiveColorSwatch(T, Rf):
    lamsnm = linspace(0.3,2.5,num=tpc.num_lams) #um
    lamsnm*=1000
    spectrumT = vstack((lamsnm, T)).T
    spectrumRf = vstack((lamsnm, Rf)).T
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
 
'''Plots color as points on the CIE colorchart based on intensity of spectra
T and Rf are arrays of intensity for transmission and front reflection'''
def plot_xy_on_fin(T, Rf):
    lamsnm = array(tpc.lams)
    lamsnm*=1000
    Tspectrum = vstack((lamsnm, T)).T
    Rfspectrum = vstack((lamsnm, Rf)).T
    xyT = give_xy(Tspectrum)
    xyRf = give_xy(Rfspectrum)
    plots.shark_fin_plot()
    #print(xyT[0])
    plt.plot(xyT[0], xyT[1], 'ro',color = 'black', label = 'T')
    plt.plot(xyRf[0], xyRf[1], 'ro',color = 'black', label = 'Rf')
    plt.annotate('T', (xyT[0], xyT[1]))
    plt.annotate('Rf', (xyRf[0], xyRf[1]))
    plt.show()
    return

'''I give a color swatch for any spectrum. Spectrum is a list of intensities and wavelength
 is the corresponding list of wavelgnths in nm. name is the title of the sample and must be a string'''
def GiveSingleColorSwatch(spectrum, wavelength, name):
    spectrumT = vstack((wavelength, spectrum)).T
    plots.spectrum_plot (spectrumT, name, 'T_Color', 'Wavelength ($nm$)', 'Intensity')
    plt.show()
    return
