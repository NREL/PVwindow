import csv
import numpy as np
from numpy import array, linspace, exp, pi, inf, vstack
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.integrate import quad,trapz
import pandas as pd
import tmm
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,figure,xlabel,ylabel,show,ylim,legend
import sys
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
from pvlib.pvsystem import singlediode

#import tmmPVColor as pvc



'''We determine the size and scaling of the photon wavelength scale. Units are um'''   
num_lams = 500
lams = linspace(0.3,2.5,num=num_lams) #um

'''We are constants and help control units'''
q = 1.602176634e-19 #coulombs. elementary charge  
c0 = 299792458 #m/s #Speed of light
hPlanck = 6.62607015e-34 #J*s   4.135667516e-15 #eV*s               
kB = 1.380649e-23 #J/K    8.61733034e-5 #eV/K  

class Layer:
    """ 
    I am a layer class for organizing data for each layer. I should make constructing stacks easier in the future and reduce possible mistakes
    """
    
    def __init__(self, thickness, fname_root, i_or_c='c', isPV=False, **kwargs):
        
        if thickness:
            self.d = thickness
        else:
            self.d = None
        
        self.i_or_c = i_or_c
        #self.nk = 1.0
        
        if fname_root:
            self.datasource = fname_root
            if kwargs.get('onecol'):
                print('dumb data')
                self.get_dumb_data()
            else:
                self.get_sensible_data()
        else:
            self.datasource = None
        
        
        self.isPV = isPV
        self.abs = 0
        
        
        
    
    
    def get_dumb_data(self):
        
        matfilename = 'Data/Materials/' + self.datasource + '.csv'
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
        
        matfilename = 'Data/Materials/' + self.datasource# + '.csv'
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
        
    def self_summary(self):
        
        print('    Material: ' + str(self.datasource) )
        print('    Thickness: ' + str(self.d) )
        print('    PV?: ' + str(self.isPV))
                 
  
class Stack:
    """
    I organize layers, interface with tmm, 
    and calculate interesting things like color, VLT, etc.
    """
    def __init__(self, layers,**kwargs):
        
        self.layers = layers

        #import data from NIST solar spectrum
        alldata = pd.read_excel('./Data/Spectra/ASTMG173.xls',header=1)

        Intensities = np.array(alldata['Direct+circumsolar W*m-2*nm-1'])
        wavelengths = np.array(alldata['Wvlgth nm'].values)
        
        self.Is = interp1d(wavelengths/1000.,Intensities*1000)

        ciedata = pd.read_csv('./Data/Spectra/CIEPhotopicLuminosity.csv',names=['lams','phis'])

        self.cieplf = interp1d(np.array(ciedata['lams'])/1000.,np.array(ciedata['phis']),bounds_error=False,fill_value=(0.0,0.0))

    def get_visible_light_transmission(self,lams,inc_angle):
        
        
        numerator = trapz(self.Is(lams)*self.cieplf(lams)*self.get_specular_RAT(lams,inc_angle)[2],lams)
        denominator = trapz(self.Is(lams)*self.cieplf(lams),lams)
        VLT = numerator/denominator
        
        #print(type(Asol.mean))
        
        return VLT
        
    
    def get_specular_RAT(self,lams,inc_angle):
        
        Ts = []
        Rs = []
        for lam in lams:
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
            Rs.append(R)
            T = (front_spol['T']+front_ppol['T']) / 2. 
            Ts.append(T)
        
        Rs = np.array(Rs)
        Ts = np.array(Ts)
        As = 1.-Rs-Ts
        return [Rs,As,Ts]
        
    def reverse(self):
        
        flippedstack = Stack(self.layers[::-1])
        return flippedstack
    
    
    def get_specular_PV_abs(self, lams, inc_angle):
        
        '''
        note from Byrnes:     
        Assumes the final layer eventually absorbs all transmitted light.
        Assumes the initial layer eventually absorbs all reflected light.
        '''
        
        
        thicks = [inf]
        iorcs = ['i']
        lnum = 0
        pvlayer = 0
        
        pvs = []
        for layer in self.layers:
            thicks.append(layer.d)
            iorcs.append(layer.i_or_c)
            pvs.append(layer.isPV)
            if layer.isPV:
                pvlayer = lnum+1 #+1 because of how tmm is written: always a layer above and below stack
            lnum += 1
        
        print('pvlayer: ' + str(pvlayer))
        print('lnum: ' +str(lnum)) 
        print(any(pvs))
        print(np.invert(any(pvs)))
        if np.invert(any(pvs)):
            print('no PV')
            return np.zeros(np.shape(lams))
        
        thicks.append(inf)
        iorcs.append('i')


        thicks_bw = thicks[::-1]
        iorcs_bw = iorcs[::-1]
        
        pvabs = []
        
        for lam in lams:

            nks = [1]
            for layer in self.layers:
                nks.append(layer.nk(lam))
            nks.append(1)

            #note the plus one because of the assumed before and after layers
            front_spol = tmm.inc_tmm('s',nks,thicks,iorcs,inc_angle,lam)
            front_ppol = tmm.inc_tmm('p',nks,thicks,iorcs,inc_angle,lam)

            pvabs_s = tmm.inc_absorp_in_each_layer(front_spol)[pvlayer]
            pvabs_p = tmm.inc_absorp_in_each_layer(front_ppol)[pvlayer]


            pvabs.append( (pvabs_s + pvabs_p) / 2. )

            
        #print(allabs)
            
        return pvabs 
    
    
    def update_from_dash(self,dashdata):
        
        ct = 0
        layers = []
        for entry in dashdata:
            
            if entry['Thickness [μm]']:
            
                if float(entry['Thickness [μm]'])>100:
                    ic = 'i'
                else:
                    ic = 'c'
            
            layer = Layer(entry['Thickness [μm]'], 
                          entry['Material'], 
                          i_or_c=ic, 
                          isPV=entry['PV'])
            layers.append(layer)
            
        self.layers = layers
       
        return False
    
    
    def self_summary(self):
        
        print('=======================================')
        print('I am a stack with the following layers:')
        print('=======================================')

        ct = 0
        for layer in self.layers:
            ct+=1
            print('  Layer ' + str(ct))
            layer.self_summary()
            print('')
    
        '''
        plt.figure()
        plt.plot(wavelengths/1000,self.cieplf(wavelengths/1000))
        plt.show()
        '''
            
    '''
    def get_solar_weighted_absorption(self,lamrange,inc_angle):
                
        
        integ = vegas.Integrator([lamrange])
        
        Asol = integ(lambda lam: self.Is(lam)*self.get_RAT(lam,inc_angle)[1], nitn=10, neval=100)[0]
        Asol /= integ(self.Is, nitn=10, neval=1000)[0]
        
        #print(type(Asol.mean))
        
        return Asol.mean
   
    
    def get_visible_light_transmission_OLD(self,lamrange,inc_angle):
        
        integ = vegas.Integrator([lamrange])
        
        numerator = integ(lambda lam: self.Is(lam)*self.cieplf(lam)*self.get_RAT(lam,inc_angle)[2], nitn=10, neval=150)[0]
        denominator = integ(lambda lam: self.Is(lam)*self.cieplf(lam), nitn=10, neval=150)[0]
        VLT = numerator/denominator
        
        #print(type(Asol.mean))
        
        return VLT.mean
    '''
   

# STUFF FROM ADAM

## wierd stuff to fix

'''Gives a spectrum of VLT. Used for plotting'''
def VLTSpectrum(layers):
    return Stack(layers)

'''Calculates Spectra Based on the layers of the cell
AbsorberLayer is an integer giving the position of the PV layer in the stack. Currently supports 1 PV layer'''
def Spectra(layers, AbsorberLayer, inc_angle):
    thicks = [inf]
    iorcs = ['i']
    for layer in layers:
        thicks.append(layer.d)
        iorcs.append(layer.i_or_c)
    thicks.append(inf)
    iorcs.append('i')
    
       
    thicks_bw = thicks[::-1]
    iorcs_bw = iorcs[::-1]

    Ts = []
    Rfs = []
    Rbs = []
    AbsByAbsorbers = []
    #EQEs2 = []
    #IREQEs = []


    layerchoice = AbsorberLayer 
    #layerchoice2 = 5

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
    
        AbsByAbsorber_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice]
        AbsByAbsorber_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
        AbsByAbsorbers.append( (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2. )
    
        # EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
        # EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
        # EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
        Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
        Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
        Ts.append( (front_spol['T']+front_ppol['T']) / 2. )


    Ts = array(Ts)
    Rfs = array(Rfs)
    Rbs = array(Rbs)
    As = 1-Ts-Rfs
    sanities = Ts+Rfs+As

    AbsByAbsorbers = array(AbsByAbsorbers)
    Spectra = {'AbsByAbsorbers':AbsByAbsorbers, 'Ts':Ts,'Rfs':Rfs,'Rbs':Rbs,'As':As,'Total':sanities}
    return Spectra


## other stuff on color


'''This code is all based on the ColorPy package. Check the online documentation
for more info and more functions if needed'''


'''Gives color swatches of reflected and transmitted light based on input of intensity of spectrum'''
def GiveColorSwatch(T, Rf):
    lamsnm = linspace(0.3,2.5,num=num_lams) #um
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
    lamsnm = array(lams)
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




## stuff on PCE
# ******************** Here I add PCE calculation *********************#
            
'''This stuff imports a spreadsheet of the solar spectrum'''
worksheet = pd.read_excel('./Data/Spectra/ASTMG173.xls')#('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
downloaded_array = array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)

# Interpolate to get a continuous function which I will be able to do integrals on:
'''Interpolated solar spectrum
when using, inputs must be within 300-2500 nm'''
AM15interp = interp1d(AM15[:,0]/1000, AM15[:,1])


# Here’s the plot, it looks correct:
'''Plot of the solar spectrum for verification'''
'''
y_values = np.array([AM15interp(x) for x in lams])
figure()
plot(lams , y_values)
xlabel("Wavelength (nm)")
ylabel("Spectral intensity (W/m$^2$/nm)")
title("Light from the sun");
show()
'''


'''I convert wavelength to energy. E_min and max are used for integration limits '''
Ephoton = hPlanck * c0 / lams *1e6 #J
E_min = min(Ephoton) #J   energy units from hPlanck
E_max = max(Ephoton) #J   energy units from hPlanck


'''I give the number of photons per......'''
def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton *1e6  #um
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2) * 1e9
'''I give the power for each......'''
def PowerPerTEA(Ephoton):
    return Ephoton * SPhotonsPerTEA(Ephoton)
'''I give the solar constant which is the W/m*2 emitted by the sun. Should be ~1000'''
def Solar_Constant(Ephoton):
    #PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    return quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
'''This is the solar constant value. It is called by optimization and used in a variety of functions here
Should always be ~1000'''
solar_constant = Solar_Constant(Ephoton)

'''I return an interpolated function of a spectrum relative to photon wavelength. Used for plotting'''
def GivelamsInterp(Parameter):
    Curve = Parameter.round(8)
    return interp1d(lams, Curve)

'''I return an interpolated function of a spectrum relative to photon energy'''
def GiveEInterp(Parameter):
    Curve = Parameter.round(8)
    return interp1d(Ephoton, Curve)

'''I give Q based on a given spectrum. Units are W/m^2
Input is a spectrum interpolated with respect to energy, E
eta should only be used if looking at a PV layer. Otherwise it is set to 1'''
def GiveQ(Spectra, eta = 1):#Spectra must be an interpolated function
        def integrand(E):
            return eta * Spectra(E) * PowerPerTEA(E)
        return quad(integrand, E_min, E_max, full_output=1)[0]       
    
'''
#trapz calcs
def GiveQ(Spectra, eta = 1):#Spectra must be an array
        integrand = eta*Spectra*PowerPerTEA(Ephoton)
        return -np.trapz(integrand, Ephoton)     
'''

'''
def GivePhotons(Spectra, eta):#Spectra must be an interpolated function
        def integrand(E):
            return eta * Spectra(E) * SPhotonsPerTEA(E)
        return quad(integrand, E_min, E_max)[0]        
'''
# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

'''I give the rate of recombination for the solar cell, Units are photons/(s*m**2)'''
def RR0(eta,Absorbed,Tcell):
    integrand = lambda E : eta * Absorbed(E) * (E)**2 / (exp(E / (kB * Tcell)) - 1)
    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return ((2 * pi) / (c0**2 * hPlanck**3)) * integral# / 1.60218e-19 #J/eV
#units = photons/(s*m**2)

'''I give the amount of energy converted to electricity in terms of photons, units are photons(s/m**2)'''
def Generated(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(E) * SPhotonsPerTEA(E)
#    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return quad(integrand, E_min, E_max, full_output=1)[0]
#units photons/(s*m**2)
'''
#Using trapezoidal rule for integration instaed of quad
#AbsByAbsorbers is an aray of intensities, not an interpolated function.
def RR0(eta,Absorbed,Tcell):
    AbsByAbsorbers = AbsByAbsorbers.round(8)
    integrand = eta * AbsByAbsorbers * (Ephoton)**2 / (np.exp(Ephoton / (kB * Tcell)) - 1)
    integral = trapz(integrand, Ephoton)
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

def Generated(eta,Absorbed):
    Absorbed = Absorbed.round(8)
    integrand = eta * Absorbed * SPhotonsPerTEA(Ephoton)
#    integral = quad(integrand, E_min, E_max, full_output=1)[0]
    return np.trapz(integrand, Ephoton)
'''

'''I use the single diode equation to return the max power of the cell in watts
Check PVlib documentation for details'''
def Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell, n = 1, Ns = 1):
    data = singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed,Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)
    return data['p_mp']

'''I calculate equilibrium tmperature of the cell assuming the cell is infinitely thin
TotalAbs is the full absorptance of the stack as an array of intensities, uninterpolated. 
Absorbed is PV layer absorptance interpolated
Temperature calculation is implicit so the numerical solver fsolve is used.
This equation is derived from Wheeler and Wheeler Detailed Balance Analysis of Photovoltaic Windows'''
def TcellCalc(TotalAbs, eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh):
    AbsTotal = GiveEInterp(TotalAbs)
    Qabs = GiveQ(AbsTotal)
    Temp = lambda Tcell: (Qabs - Give_Pmp(eta,Absorbed,Rs,Rsh, Tcell) + Ui*Ti + Uo*To)/(Ui + Uo)-Tcell
    return fsolve(Temp, 300)[0]


'''I use the single diode equation to produce an IV curve and power plot
I also return related values such as Voc, Isc, and Pmp in units volts, amps, and watts
See pvlib singlediode equation for more information'''
def GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1):
    data = singlediode(Generated(eta, Absorbed)*q, RR0(eta, Absorbed, Tcell)*q, Rs, Rsh, n*Ns*kB*Tcell/q, ivcurve_pnts = 500)


    Vvalues = array(data['v'])
    Ivalues = array(data['i'])
    #print('Isc = ', Isc, ', Voc = ', Voc, ', Imp = ', Imp, ', Vmp = ', Vmp, ', Pmp =', Pmp)

    figure()
    plot(Vvalues,Ivalues, label = 'IV')
    xlabel('Voltage, (V)')
    ylabel('Current (A) or Power (W/m^2)')
    ylabel('Power (W/m^2)')
    P_values = array([Ivalues * Vvalues])
    plot(Vvalues , P_values.T, label = 'Power')
    ylim(-1, 150)
    legend(loc = 'upper right')
    show()
    return data



'''I give the solar heat gain coefficient. unitless numebr between 0 and 1
Ts is the transmission spectra. Must be a list of intensities, not an interpolated function
This equation comes form a combination of Wheeler and Wheeler Detailed Balance Analysis of Photovoltaic Windows
and equation 3.18 from Fundamentals of Heat and Mass Transfer 6ed Incropera'''
def SHGC(Ts, Ti, To, Tcell, Ui):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    Rtot = 1/Ui #This is approximate because Ui is assumed
    #Included in GiveQ for simplicity but should not be used for calculating SHGC
    TransTotal = GiveEInterp(Ts)
    Qtrans = GiveQ(TransTotal,1)
    return (Qtrans + Ui*(Tcell-Ti) - ((To-Ti)/Rtot))/solar_constant

'''I give max efficiency also called PCE'''
'''Absorbed must be an interpolated function of the absorption spectrum of the PV layer'''
def max_efficiency(eta,Absorbed,Tcell, Rs, Rsh):
    #Tcell = TcellCalc(As,Ti,To,eta,Absorbed)
    return Give_Pmp(eta, Absorbed, Rs, Rsh, Tcell) / solar_constant

'''We determine the incident angle of the sun shining on the cell. Input is in degrees'''
def giveincangle(angle):
    degree = pi/180
    return angle*degree

'''I assemble a list of layer objects using Thicknesses and Materials''' 
def GiveLayers(Thickness,Materials):
    x = len(Materials)
    if x == len(Thickness):
        Layers = []
        for i in range(x):
            Layers.append(Materials[i](Thickness[i]))
        return Layers
    else:  
        raise ValueError ('layers and Thickness lengths do not match')

'''I give important info about a solar cell such as PCE, SHGC, Temperature, etc'''
def GiveImportantInfo(Thickness, Materials,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,Angle=0):
    global inc_angle
    inc_angle = giveincangle(Angle)
    
    layers = GiveLayers(Thickness,Materials)
    
    spectra = Spectra(layers ,AbsorberLayer)
    AbsByAbsorbers = spectra['AbsByAbsorbers']
    Ts = spectra['Ts']
    Rfs = spectra['Rfs']
    Rbs = spectra['Rbs']
    As = spectra['As']
    sanities = spectra['Total']
    Absorbed = GiveEInterp(AbsByAbsorbers)
    VLTcalc =  cvs.getVLT(Ts,lams)#VLT(layers)
    Tcell = TcellCalc(As,eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
    #Absorbed = tpc.GiveEInterp(tpc.Spectra(tpc.GiveLayers(Thickness, Materials),4)['AbsByAbsorbers'])
    data = GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1)
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    SHGCcalc = SHGC(Ts, Ti, To, Tcell, Ui)
    PCE = max_efficiency(eta,Absorbed,Tcell, Rs, Rsh)


    #Spectral Curves
    figure()
    plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
    plot(lams,Ts,color='green',marker=None,label="$T$")
    plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
    plot(lams,As,color='black',marker=None,label="A")
    plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
    plot(lams,sanities,color='gold',marker=None,label="R+A+T")
    plot(lams,VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    xlabel('wavelength, $\mu$m')
    ylabel('Intensity')
    legend(loc = 'upper right')
    show()
    
    EphotoneV = Ephoton*6.241509e+18 
    figure()
    plot(EphotoneV, Ts, color='magenta',marker=None,label="$T$")
    plot(EphotoneV, Rfs,color='green',marker=None,label="$R_f$")
    plot(EphotoneV, Rbs,color='orange',marker=None,label="$R_b$")
    plot(EphotoneV, AbsByAbsorbers,color='black',marker=None,label="Abs")
    #plot(Ephoton,tpc.VLTSpectrum(layers).cieplf(lams),color='red',marker=None,label="photopic")
    legend(loc = 'upper right')
    xlabel('Energy, eV')
    ylabel('Intensity')
    show()

    GiveColorSwatch(Ts, Rfs)
    plot_xy_on_fin(Ts, Rfs)

    print('PCE = ',PCE,'VLT = ', VLTcalc, 'SHGC = ',SHGCcalc, 'Tcell = ',Tcell)#,'time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
    return {'PCE':PCE, 'VLT':VLTcalc, 'SHGC':SHGCcalc, 'Tcell':Tcell,'Isc':Isc, 'Voc': Voc, 'Imp': Imp, 'Vmp': Vmp,'Pmp': Pmp}

# Vince hacks some stuff from Adam

def get_performance_characteristics(stack,eta,Ti,To,Ui,Uo,Rs,Rsh,AbsorberLayer,Angle):
    
    layers = stack.layers
    
    spectra = Spectra(layers ,AbsorberLayer,Angle)
    AbsByAbsorbers = spectra['AbsByAbsorbers']
    Ts = spectra['Ts']
    Rfs = spectra['Rfs']
    Rbs = spectra['Rbs']
    As = spectra['As']
    sanities = spectra['Total']
    Absorbed = GiveEInterp(AbsByAbsorbers)
    VLTcalc =  stack.get_visible_light_transmission(lams,Angle) #cvs.getVLT(Ts,lams)#VLT(layers)
    Tcell = TcellCalc(As,eta, Ti,To, Absorbed, Ui, Uo, Rs, Rsh)
    #Absorbed = tpc.GiveEInterp(tpc.Spectra(tpc.GiveLayers(Thickness, Materials),4)['AbsByAbsorbers'])
    data = GiveIVData(eta, Absorbed, Rs, Rsh,Tcell, n = 1, Ns = 1)
    Isc = data['i_sc']
    Voc = data['v_oc']
    Imp = data['i_mp']
    Vmp = data['v_mp']
    Pmp = data['p_mp']
    SHGCcalc = SHGC(Ts, Ti, To, Tcell, Ui)
    PCE = max_efficiency(eta,Absorbed,Tcell, Rs, Rsh)

    #print('PCE = ',PCE,'VLT = ', VLTcalc, 'SHGC = ',SHGCcalc, 'Tcell = ',Tcell)#,'time to calculate PCE from scratch in seconds = ', TimePCE, 'Time to run optimizer in minutes = ',TimeOptimize/60)
    return {'PCE':PCE, 'VLT':VLTcalc, 'SHGC':SHGCcalc, 'Tcell':Tcell,'Isc':Isc, 'Voc': Voc, 'Imp': Imp, 'Vmp': Vmp,'Pmp': Pmp}

