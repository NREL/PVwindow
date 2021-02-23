import csv
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tmm
import vegas

class Layer:
    """ 
    I am a layer class for organizing data for each layer. I should make constructing stacks easier in the future and reduce possible mistakes
    """
    
    def __init__(self, thickness, fname_root, i_or_c,**kwargs):
        self.d = thickness
        self.i_or_c = i_or_c
        #self.nk = 1.0
        self.datasource = fname_root
        
        if kwargs.get('onecol'):
            print('dumb data')
            self.get_dumb_data()
        else:
            self.get_sensible_data()
        
    
    
    def get_dumb_data(self):
        
        matfilename = 'Data/' + self.datasource + '.csv'
        lct = 0
        bothdat = []
        with open(matfilename, newline='') as csvfile:
            rawdat = csv.reader(csvfile, delimiter=' ')
            for row in rawdat:
                lct += 1
                if row:
                    bothdat.append(row[0])
                    if 'k' in row[0]:
                        kstart = lct
            lct = 1
            nlams = []
            ns = []
            klams = []
            ks = []
            for line in bothdat:
                if lct < kstart-1:
                    if 'n' not in line:
                        nlams.append(float(line.split(',')[0]))
                        ns.append(float(line.split(',')[1]))
                else:
                    if 'k' not in line:
                        #print(line)
                        klams.append(float(line.split(',')[0]))
                        ks.append(float(line.split(',')[1]))
                lct += 1
        nlams = np.array(nlams)
        ns = np.array(ns)
        #print(nlams)
        klams = np.array(klams)
        ks = np.array(ks)
        #print(klams)
        
        self.n = interp1d(nlams,ns,fill_value="extrapolate")
        self.k = interp1d(klams,ks,fill_value="extrapolate")
        

        
    def get_sensible_data(self):
        """
                next we will unpack n and k data from a csv file and turn it into a callable interpolation function
        """
        
        matfilename = 'Data/' + self.datasource + '.csv'
        testdat = np.genfromtxt(matfilename,delimiter=',',skip_header=1)
        
        nlams = testdat[:,0]
        ns = testdat[:,1]
        ks = testdat[:,2]
        
        self.n = interp1d(nlams,ns,fill_value="extrapolate")
        self.k = interp1d(nlams,ks,fill_value="extrapolate")
        
        
    def nk(self,lam):
        return complex(self.n(lam),self.k(lam))
    
    def plotnk(self,lams):
        
        plt.figure()
        plt.plot(lams, self.n(lams),label='n')
        plt.plot(lams, self.k(lams),label='k')
        plt.title(self.datasource)
        plt.legend()
        plt.show()
                 
  
class Stack:
    """
    I organize layers, interface with tmm, 
    and calculate interesting things like color, VLT, etc.
    """
    def __init__(self, layers,**kwargs):
        
        self.layers = layers

        #import data from NIST solar spectrum
        alldata = pd.read_excel('./Data/ASTMG173.xls',header=1)

        Intensities = np.array(alldata['Direct+circumsolar W*m-2*nm-1'])
        wavelengths = np.array(alldata['Wvlgth nm'].values)
        
        self.Is = interp1d(wavelengths/1000.,Intensities*1000)

        ciedata = pd.read_csv('./Data/CIEPhotopicLuminosity.csv',names=['lams','phis'])

        self.cieplf = interp1d(np.array(ciedata['lams'])/1000.,np.array(ciedata['phis']),bounds_error=False,fill_value=(0.0,0.0))
        
        '''
        plt.figure()
        plt.plot(wavelengths/1000,self.cieplf(wavelengths/1000))
        plt.show()
        '''
            
    
    def get_solar_weighted_absorption(self,lamrange,inc_angle):
                
        
        integ = vegas.Integrator([lamrange])
        
        Asol = integ(lambda lam: self.Is(lam)*self.get_RAT(lam,inc_angle)[1], nitn=10, neval=100)[0]
        Asol /= integ(self.Is, nitn=10, neval=1000)[0]
        
        #print(type(Asol.mean))
        
        return Asol.mean
    
    def get_visible_light_transmission(self,lamrange,inc_angle):
        
        integ = vegas.Integrator([lamrange])
        
        numerator = integ(lambda lam: self.Is(lam)*self.cieplf(lam)*self.get_RAT(lam,inc_angle)[2], nitn=10, neval=100)[0]
        denominator = integ(lambda lam: self.Is(lam)*self.cieplf(lam), nitn=10, neval=100)[0]
        VLT = numerator/denominator
        
        #print(type(Asol.mean))
        
        return VLT.mean
        
        
    
    def get_RAT(self,lam,inc_angle):
        
        thicks = [tmm.inf]
        iorcs = ['i']
        nks = [1]
        for layer in self.layers:
            nks.append(layer.nk(lam))
            thicks.append(layer.d)
            iorcs.append(layer.i_or_c)
        thicks.append(tmm.inf)
        iorcs.append('i')
        nks.append(1)
        
        front_spol = tmm.inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
        front_ppol = tmm.inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)
        
        R = (front_spol['R']+front_ppol['R']) / 2.
        T = (front_spol['T']+front_ppol['T']) / 2. 
        
        return [R,1-R-T,T]
        
    def reverse(self):
        
        flippedstack = Stack(self.layers[::-1])
        return flippedstack
        
