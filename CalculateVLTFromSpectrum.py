# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:44:40 2021

@author: aduell
"""

'''
This file takes a data file of transmission spectra and gives a new file
containing the VLT for each spectrum

Datafile format: .xlsx or.xls, First row is header, 1st column is wavlength in um
All other rows are intensities of transmission ranging from 0-100
If the datafile has trans intensity from 0-1, go to getVLT, Transmission, 
and change TS/100 to TS

'''


'''_____________________Input Things Here________________________________'''

#DataFileLocation='./YourDataFileHere.xlsx'
#OutputTitle='./YourVLTDataHere'












from numpy import array
import pandas as pd
from scipy.interpolate import interp1d
import sys
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
import vegas




'''We are constants and help control units'''
q = 1.602176634e-19 #coulombs. elementary charge  
c0 = 299792458 #m/s #Speed of light
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  


'''I give the AM1.5G aka the solar spectrum and Photopic Eye Response data'''
def AM15GandPER():
    alldata = pd.read_excel('./Data/ASTMG173.xls',header=1)
    
    '''I provide the AM1.5G data'''
    Intensities = array(alldata['Direct+circumsolar W*m-2*nm-1'])
    wavelengths = array(alldata['Wvlgth nm'].values)
    AM15G = interp1d(wavelengths/1000.,Intensities*1000)
    
    '''I provide the photopic eye response'''
    ciedata = pd.read_csv('./Data/CIEPhotopicLuminosity.csv',names=['lams','phis'])
    cieplf = interp1d(array(ciedata['lams'])/1000.,array(ciedata['phis']),bounds_error=False,fill_value=(0.0,0.0))
    AMPER = [AM15G,cieplf]
    
    return AMPER

'''
AMPER = AM15GandPER()
AM15G = AMPER[0]
cieplf = AMPER[1]

#I perform the VLT calculation when given a spectrum and associated wavelength range
def getVLT(TS,lamrange):
    integ = vegas.Integrator([lamrange])
    Transmission = interp1d(lamrange,TS/100)
    numerator = integ(lambda lam: AM15G(lam)*cieplf(lam)*Transmission(lam), nitn=10, neval=200)[0]
    denominator = integ(lambda lam: AM15G(lam)*cieplf(lam), nitn=10, neval=200)[0]
    VLT = numerator/denominator
    return VLT.mean
'''


'''I perform the VLT calculation when given:
    a spectrum and associated wavelength range, cieplf, and AM1.5G data.'''
def getVLTComplete(TS,lamrange,AM15G,cieplf):
    integ = vegas.Integrator([lamrange])
    Transmission = interp1d(lamrange,TS/100)
    numerator = integ(lambda lam: AM15G(lam)*cieplf(lam)*Transmission(lam), nitn=10, neval=200)[0]
    denominator = integ(lambda lam: AM15G(lam)*cieplf(lam), nitn=10, neval=200)[0]
    VLT = numerator/denominator
    return VLT.mean

'''I compile a lit of VLTs given an array of transmission data'''

def BuildVLTList(Data): 
    Tlams = Data[:,0]
    Tlams*=1/1000
    AMPER = AM15GandPER()
    AM15G = AMPER[0]
    cieplf = AMPER[1]
    x = Data.shape[1]
    VLTList = []
    for i in range(x-1):
        VLTList.append(getVLTComplete(Data[:,i+1],Tlams,AM15G,cieplf))
        print('...working...')
    VLTList= array(VLTList)
    return VLTList


'''I import a data file given a location and file name, extract the column names, and neatly package it'''
def ImportDataFile(DataFileLocation):
    Data = pd.read_excel(DataFileLocation)
    ColName = list(Data.columns)
    ColName.pop(0)
    Data = array(Data)
    return {'Data':Data,'ColName':ColName}


'''I export the VLT data as a .csv file. Title must be a string eg: 'Title'   '''
def ExportVLTData(Title,ColName,VLTList):
    dict = {'Sample': ColName, 'VLT': VLTList}  
    df = pd.DataFrame(dict) 
    df.to_csv(Title+'.csv') 
    VLTcsv =pd.read_csv(Title+'.csv',header=1)
    print(VLTcsv)

'''
I produce a file containing all the calculated VLTs.
Inputs are the location and name of a data file
and the desired title of the output file
Both inputs must be strings

See top of file for Datafile Format
'''
def GiveMeVLTFile(DataFileLocation, OutputTitle):
    Data = ImportDataFile(DataFileLocation)
    VLTList = BuildVLTList(Data['Data'])
    ExportVLTData(OutputTitle,Data['ColName'],VLTList)
    return

#GiveMeVLTFile(DataFileLocation,OutputTitle)
#GiveMeVLTFile('./Data/TransmissionSpectra.xlsx','./TestFile')


