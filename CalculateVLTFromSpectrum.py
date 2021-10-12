# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:44:40 2021

@author: aduell
"""

'''
This file takes a data file of transmission spectra and gives a new file
containing the VLT for each spectrum

Datafile format: .xlsx or.xls, First row is header, 1st column is wavlength in nm
All other rows are intensities of transmission ranging from 0-100
If the datafile has trans intensity from 0-1, go to getVLTComplete--> Transmission, 
and change TS/100 to TS. Additionally go to BuildVltList, find the line
SpectraList.append(vstack((Tlamsnm, Data[:,i+1]/100)).T) and remove the '/100'
Preferably just multiply all the intensities by 100 so they are on a scale from 0-100

Any negative vlaues in the transmission spectra will be converted to 0

when calculating VLT, neval can be increased to increase the precision of the calculation
at the cost of speed

Data file must have a wavelength range within 280-4000 nm 
Wavelength range will be cropped from 370 to 770 nm

getVLT is easier to import and use in other programs but is not currently used for 
producing a data file here

All color calculations use ColorPy. Check the online documentation for more info
and additional functions. http://markkness.net/colorpy/ColorPy.html

'''


'''_____________________Input Things Here________________________________'''

#DataFileLocation='./YourDataFileHere.xlsx'
#OutputTitle='./YourVLTDataHere'












from numpy import array, where, vstack
import pandas as pd
from scipy.interpolate import interp1d
import sys
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
import vegas
#import tmmPVColor as pvc
from colorpy import ciexyz, colormodels #need to install colorpy to call all packages at once. 
#If colorpy doesnt work, go to the file and change import 'blank' to from . import 'blank'
#from wpv import Layer, Stack





'''We are constants and help control units'''
q = 1.602176634e-19 #coulombs. elementary charge  
c0 = 299792458 #m/s #Speed of light
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  




'''I give the AM1.5G aka the solar spectrum and Photopic Eye Response data'''
def AM15GandPER():
    alldata = pd.read_excel('./Data/ASTMG173.xls',header=1)
    
    '''I provide the AM1.5G data (AM)'''
    Intensities = array(alldata['Direct+circumsolar W*m-2*nm-1'])
    wavelengths = array(alldata['Wvlgth nm'].values)
    AM15G = interp1d(wavelengths/1000.,Intensities*1000)
    
    '''I provide the photopic eye response (PER)'''
    ciedata = pd.read_csv('./Data/CIEPhotopicLuminosity.csv',names=['lams','phis'])
    #global lamrange
    #lamrange = array(ciedata['lams'])/1000.
    cieplf = interp1d(array(ciedata['lams'])/1000.,array(ciedata['phis']),bounds_error=False,fill_value=(0.0,0.0))
    AMPER = [AM15G,cieplf]
    
    return AMPER


AMPER = AM15GandPER()
AM15G = AMPER[0]
cieplf = AMPER[1]


'''I perform the VLT calculation when given a spectrum and associated wavelength range.
Adapted from wpv code by Vince Wheeler
TSScale tells what the transmission values can span. eg 0-100 set as 100, 0-1 set as 1'''
def getVLT(TS,lamrange, TSScale = 1):
    integ = vegas.Integrator([lamrange])
    Transmission = interp1d(lamrange,TS/TSScale)
    numerator = integ(lambda lam: AM15G(lam)*cieplf(lam)*Transmission(lam), nitn=10, neval=300)[0]
    denominator = integ(lambda lam: AM15G(lam)*cieplf(lam), nitn=10, neval=300)[0]
    VLT = numerator/denominator
    return VLT.mean


'''I perform the VLT calculation when given:
    a spectrum and associated wavelength range, cieplf, and AM1.5G data.
    I get used by BuildVLTList'''
def getVLTComplete(TS,lamrange, AM15G,cieplf):
    integ = vegas.Integrator([lamrange])
    Transmission = interp1d(lamrange,TS/100)
    numerator = integ(lambda lam: AM15G(lam)*cieplf(lam)*Transmission(lam), nitn=10, neval=700)[0]
    denominator = integ(lambda lam: AM15G(lam)*cieplf(lam), nitn=10, neval=700)[0]
    VLT = numerator/denominator
    return VLT.mean


'''I give the CIE xyz values for a spectrum that is nm vs intensity
I get used in BuildVLTList'''
def give_xyz(spectrum):
    xyz = ciexyz.xyz_from_spectrum(spectrum)
    xyz1 = colormodels.xyz_normalize (xyz)
    #xyz2 = xyz1[0:2]
    #print(xyz2)
    return(xyz1)


'''I take a list of arrays for CIE data and separate them into 3 lists
1 for x values, 1 for y values, and 1 for z values
I get used in BuildVLTList to process the CIE data from give_xyz'''
def SeparateCIEData(CIEData):
    w = len(CIEData)
    x=[]
    y=[]
    z=[]
    for i in range(w):
        x.append(CIEData[i][0])
    for i in range(w):
        y.append(CIEData[i][1])
    for i in range(w):
        z.append(CIEData[i][2])
    return [x,y,z]


'''I compile a lit of VLTs given an array of transmission data'''
def BuildVLTList(Data): 
    #global Tlams
    Tlamsnm = Data[:,0]
    Tlams = array(Tlamsnm)
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
    
    '''Color calcs happen here'''
    SpectraList = []
    '''I assemble the raw data into spectra usable by give_xy'''
    for i in range(x-1):
        SpectraList.append(vstack((Tlamsnm, Data[:,i+1]/100)).T)
    SpectraList=array(SpectraList)
    CIEColorList = []
    '''I calcualte the CIE xy values'''
    for i in range(x-1):
        CIEColorList.append(give_xyz(SpectraList[i]))
    CIEColorList = array(CIEColorList)
    CIEData = SeparateCIEData(CIEColorList)
    
    return [VLTList,CIEData]

    
'''I import a data file given a location and file name, extract the column names, and neatly package it'''
'''I also convert negative transmission values to 0''' 
def ImportDataFile(DataFileLocation):
    Data = pd.read_excel(DataFileLocation)
    ColName = list(Data.columns)
    ColName.pop(0)
    Data = array(Data)
    Data = where(Data < 0, 0, Data)
    return {'Data':Data,'ColName':ColName}


'''I export the VLT data as a .csv file. Title must be a string eg: 'Title'   '''
def ExportVLTData(Title,ColName,VLTList):
    dict = {'Sample': ColName, 'VLT': VLTList[0], 'CIE x':VLTList[1][0], 'CIE y':VLTList[1][1], 'CIE z':VLTList[1][2]}  
    df = pd.DataFrame(dict) 
    df.to_csv(Title+'.csv') 
    VLTcsv =pd.read_csv(Title+'.csv',header=1)
    print(VLTcsv)

'''
I produce a file containing all the calculated VLTs. I also contain cie xyz data
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

#Data = ImportDataFile('./VLT_Data/ColorBalanceSpectra.xlsx')
#Data = pd.read_excel('./VLT_Data/Transparent_KS_device_IZO_10.xlsx')
#Data = Data['Data']

#GiveMeVLTFile(DataFileLocation,OutputTitle)
#GiveMeVLTFile('./VLT_Data/BAR2115_1_MeOH_TransmissionData.xlsx','./VLT_Data/VLT_Output/BAR2115_1_MeOH_VLTData')
#GiveMeVLTFile('./VLT_Data/BAR2115_2_H2O_TransmissionData.xlsx','./VLT_Data/VLT_Output/BAR2115_2_H2O_VLTData')
#GiveMeVLTFile('./VLT_Data/Transparent_KS_device_IZO_10.xlsx','./VLT_Data/VLT_Output/Transparent_KS_device_IZO_10_VLTData')
#GiveMeVLTFile('./VLT_Data/DeviceTransmissionData.xlsx','./VLT_Data/VLT_Output/DeviceVLTData')
#GiveMeVLTFile('./VLT_Data/Static_PV_Semitransparant_TransmissionData_Kevin.xlsx','./VLT_Data/VLT_Output/Static_PV_Semitransparant_VLTData_Kevin')

#GiveMeVLTFile('./VLT_Data/ColorBalanceSpectra.xlsx', './VLT_Data/VLT_Output/ColorBalanceVLTandColorData')
    