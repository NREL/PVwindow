import numpy as np
import tmm
import pandas as pd
#import tmm_vw as tmm
import matplotlib.pyplot as plt
from wpv import Layer, Stack
import scipy.interpolate, scipy.integrate, pandas, sys
from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e
assert sys.version_info >= (3,6), 'Requires Python 3.6+'


# This whole thing uses microns for length

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
ClAlPc = Layer(0.200,'nkClAlPc','c')
C60 = Layer(0.100,'nkC60','c')
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
#layers = [Glass,FTO,TiO2,C60,ClAlPc,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

# 50% VLT with non-wavelength-selective absorber, MAPbBr3 = 500 nm
layers = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

# Different thicknesses of MAPI: 50% VLT = 40 nm, 25% VLT = 130 nm, 5% VLT = 370 nm, 0.5% VLT = 775 nm
#layers = [Glass,FTO,TiO2,MAPI,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
# Here's the corresponding bleached layers for 5 and 0.5%
layers = [Glass,FTO,TiO2,Bleach,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

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


# ******************** Here I add PCE calculation *********************#

# Solar cell temperature is 300 kelvin:

Tcell = 300 * K

# worksheet = pandas.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
#worksheet = pandas.read_excel('/Users/lwheeler/Code/pv-window-bem/Data/astmg173.xls')
worksheet = pandas.read_excel('./Data/ASTMG173.xls')
downloaded_array = np.array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)


# Tack on the appropriate units:
AM15[:,0] *= nm
AM15[:,1] *= W / m**2 / nm
    

λ_min = 300 * nm
λ_max = 2500 * nm
#λ_min = 280 * nm
#λ_max = 4000 * nm
E_min = hPlanck * c0 / λ_max
E_max = hPlanck * c0 / λ_min

# Interpolate to get a continuous function which I will be able to do integrals on:

AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])

# Here’s the plot, it looks correct:

λs = np.linspace(λ_min, λ_max, num=500)
y_values = np.array([AM15interp(x) for x in λs])
plt.figure()
plt.plot(λs / nm , y_values / (W/m**2/nm))
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral intensity (W/m$^2$/nm)")
plt.title("Light from the sun");
plt.show()

def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2)

PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
#print(solar_constant / (W/m**2))

#def AbsVsE(Ephoton):
#    λ = hPlanck * c0 / Ephoton 
#    return AbsInterp(λ)
#return AbsInterp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2)
    
# def EQE(eta,AbsInterp):
#    lam = hPlanck * c0 / Ephoton 
#    return eta * AbsInterp(lam)

# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

def RR0(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(hPlanck * c0 / E) * E**2 / (np.exp(E / (kB * Tcell)) - 1)
    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

def Generated(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(hPlanck * c0 / E) * SPhotonsPerTEA(E)
#    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]

def current_density(voltage, eta,Absorbed):
    return e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * voltage / (kB * Tcell)))

# ## Open-circuit voltage, short-circuit current

def JSC(eta,Absorbed):
    return current_density(0, eta,Absorbed)
def VOC(eta,Absorbed):
    return (kB * Tcell / e) * np.log(Generated(eta,Absorbed) / RR0(eta,Absorbed))

# SciPy only comes with minimization, not maximization. Let's fix that...
from scipy.optimize import fmin

def fmax(func_to_maximize, initial_guess=0.5*V):
    """return the x that maximizes func_to_maximize(x)"""
    func_to_minimize = lambda x : -func_to_maximize(x)
    return fmin(func_to_minimize, initial_guess, disp=False)[0]

def V_mpp(eta,Absorbed):
    """ voltage at max power point """
    return fmax(lambda voltage : voltage * current_density(voltage, eta,Absorbed))

def J_mpp(eta,Absorbed):
    """ current at max power point """
    return current_density(V_mpp(eta,Absorbed), eta,Absorbed)

def Power(voltage, eta,Absorbed):
    return voltage * current_density(voltage, eta,Absorbed)

def max_power(eta,Absorbed):
    voltage = V_mpp(eta,Absorbed)
    return voltage * current_density(voltage, eta,Absorbed)

def max_efficiency(eta,Absorbed):
    return max_power(eta,Absorbed) / solar_constant

#Vs = np.linspace(0.1 * V, 2 * V, num=500)
#y_values = np.array([Power(x,0.9,AbsInterp) for x in Vs])
#plt.figure()
#plt.plot(Vs / V , y_values / (W/m**2))
#plt.xlabel("Voltage (V)")
#plt.ylabel("Power (W m-2)")
#plt.ylim(-1, 150)
#plt.title("Light from the sun");
#plt.show()


# array of angles
degree = np.pi/180
#inc_angle = 0.*degree
num_angles = 50
angles = np.linspace(0,89.999*degree,num=num_angles)
   
# array of wavelengths     
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams)

# declare arrays
#Ts = []
#Rfs = []
#Rbs = []
#AbsByAbsorbers = []
#AbsByAbsorber = np.array(num_lams)
#EQEs2 = []
#IREQEs = []
PCEs = []

# This is the electron-hole pair extraction efficiency. I could also loop this.
eta = 0.8

# This parameter defines which part of the stack of materials is the absorber 
# layer for PCE calculation:

layerchoice = 4
# This is the second layer for a selective absorber stack
layerchoice2 = 5

thicks = [tmm.inf]
iorcs = ['i']
for layer in layers:
    thicks.append(layer.d)
    iorcs.append(layer.i_or_c)
thicks.append(tmm.inf)
iorcs.append('i')

thicks_bw = thicks[::-1]
iorcs_bw = iorcs[::-1]

for angle in angles:
    # Start arrays fresh
    Ts = []
    Rfs = []
    Rbs = []
    AbsByAbsorbers = []
    AbsByAbsorbers2 = []
    
    for lam in lams:

        nks = [1]
        for layer in layers:
            nks.append(layer.nk(lam))
        nks.append(1)

        nks_bw = nks[::-1]
    
        front_spol = tmm.inc_tmm('s',nks,thicks,iorcs,angle,lam)
        front_ppol = tmm.inc_tmm('p',nks,thicks,iorcs,angle,lam)
        back_spol = tmm.inc_tmm('s',nks_bw,thicks_bw,iorcs_bw,angle,lam)
        back_ppol = tmm.inc_tmm('p',nks_bw,thicks_bw,iorcs_bw,angle,lam)
    
        # EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
        # EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
        # EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
        Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
        Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
        Ts.append( (front_spol['T']+front_ppol['T']) / 2. )
           
        AbsByAbsorber_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice]
        AbsByAbsorber_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
        AbsByAbsorbers.append( (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2. )
        
        # I'll add the absorption of a second layer for when I use the two layers
        # for selective absorption
        AbsByAbsorber2_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
        AbsByAbsorber2_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
        AbsByAbsorbers2.append( (AbsByAbsorber2_spol + AbsByAbsorber2_ppol) / 2. )
        #np.append(AbsByAbsorbers, (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2.)
    
    # Comment in/out AbsByAbsrober2 on the next few lines to include a second
    # abosrber layer in the calculation
    AbsByAbsorber = np.array(AbsByAbsorbers)
    #AbsByAbsorber2 = np.array(AbsByAbsorbers2)
    #AbsByAbsorber = AbsByAbsorber + AbsByAbsorber2

    AbsByAbsorber = AbsByAbsorber.round(8)
    AbsInterp = scipy.interpolate.interp1d(lams *1000 * nm, AbsByAbsorber)
        
    PCEs.append(max_efficiency(eta,AbsInterp))
    
    print("angle =",angle)
    print("PCE =",max_efficiency(eta, AbsInterp))
    
    # Here I calculate VLT and spit it out to the screen
    VLTstack=Stack(layers)
    VLT=VLTstack.get_visible_light_transmission(lams,angle)
    print("VLT =",VLT)

Vs = np.linspace(0.1 * V, 2 * V, num=500)
y_values = np.array([Power(x,eta,AbsInterp) for x in Vs])
plt.figure()
plt.plot(Vs / V , y_values / (W/m**2))
plt.xlabel("Voltage (V)")
plt.ylabel("Power (W m-2)")
plt.ylim(-1, 150)
#plt.title("Light from the sun");
plt.show()


Ts = np.array(Ts)
Rfs = np.array(Rfs)
Rbs = np.array(Rbs)
PCEs = np.array(PCEs)
As = 1-Ts-Rfs
sanities = Ts+Rfs+As

AbsByAbsorbers = np.array(AbsByAbsorbers)
#EQEs2 = np.array(EQEs2)
#IREQEs=EQEs+EQEs2

# I need to get a continuous function of the fraction of the photons absorbed
# in the absorber layer. Here I tack on units and interpolate:
    
# X[:,0] *= 1000 * nm

# AbsInterp = scipy.interpolate.interp1d(X[:,0], X[:,1])

# Tack on units
lams *= 1000 * nm

# Round AbsByAbsorber to make really small numbers equal to zero
AbsByAbsorbers = AbsByAbsorbers.round(8)
AbsInterp = scipy.interpolate.interp1d(lams, AbsByAbsorbers)

#λs = np.linspace(λ_min, λ_max, num=500)
#Abs_values = np.array([AbsInterp(x) for x in λs])
Abs_values = np.array([AbsInterp(x) for x in lams])
plt.figure()
#plt.plot(λs / nm , Abs_values )
plt.plot(lams / nm , Abs_values )
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorptance")
plt.title("Absorbed in absorber layer");
plt.show()

#

X = np.transpose([lams,AbsByAbsorbers])
np.savetxt('./Output/AbsByAbsorber.txt',X,delimiter=',',header="wavelength [micron], AbsByAbsorber [1]")

Y = np.transpose([lams,Ts,Rfs,Rbs])
np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")

PCEangleData = np.transpose([angles/degree,PCEs])
np.savetxt('./Output/PCE_v_angle.txt',PCEangleData,delimiter=',',header="Angle [degrees], PCE")

# This is for when there are 2 layers contributing to the AbsByAbsorber:
#Z = np.transpose([lams,IREQEs])
#np.savetxt('./Output/IREQE.txt',Z,delimiter=',',header="wavelength [micron], EQE [1]")

plt.figure()
plt.plot(lams,Rfs,color='magenta',marker=None,label="$R_f$")
plt.plot(lams,Ts,color='green',marker=None,label="$T$")
plt.plot(lams,Rbs,color='purple',marker=None,label="$R_b$")
plt.plot(lams,As,color='black',marker=None,label="A")
plt.plot(lams,AbsByAbsorbers,color='black',linestyle='--',marker=None,label="AbsByAbsorber")
# This is for when there are 2 layers contributing to the EQE:
#plt.plot(lams,IREQEs,color='gray',linestyle='--',marker=None,label="EQE")
plt.plot(lams,sanities,color='gold',marker=None,label="R+A+T")
# This is the photopic eye response
plt.plot(lams,VLTstack.cieplf(lams),color='red',marker=None,label="photopic")
# This is the solar spectrum
 #plt.plot(lams,VLTstack.Is(lams)/max(VLTstack.Is(lams)),color='gray',marker=None,label="AM1.5")
plt.xlabel('wavelength, $\mu$m')
plt.legend(loc = 'upper right')
plt.show()


#print("PCE =",max_efficiency(0.9,AbsInterp))
