'''
Adam's version of the tmmPCE file
trying to incorporate shunt resistance and series resistance
Testing on how to solve a complex implicit equation
'''

import numpy as np
import tmm
import pandas as pd
#import tmm_vw as tmm
import matplotlib.pyplot as plt
from wpv import Layer, Stack
import scipy.interpolate, scipy.integrate, pandas, sys
import scipy
from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e, A, C, ohm, J, um
import sympy
import sympy.solvers.solvers
assert sys.version_info >= (3,6), 'Requires Python 3.6+'

#import numericalunits
#help(numericalunits)

# This whole thing uses microns for length

degree = np.pi/180
#inc_angle = 10.*degree
inc_angle = 0.*degree
        
num_lams = 500
lams = np.linspace(0.3,2.5,num=num_lams)

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
#layers = [Glass,FTO,TiO2,C60,ClAlPc,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

# 50% VLT with non-wavelength-selective absorber, MAPbBr3 = 500 nm
layers = [Glass,FTO,TiO2,MAPBr,NiO,ITO,EVA,Glass,TiO2lowE,Ag,TiO2lowE]

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
AbsByAbsorbers = []
#EQEs2 = []
#IREQEs = []

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
    
    AbsByAbsorber_spol = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice]
    AbsByAbsorber_ppol = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice]
    
    AbsByAbsorbers.append( (AbsByAbsorber_spol + AbsByAbsorber_ppol) / 2. )
    
   # EQE_spol2 = tmm.inc_absorp_in_each_layer(front_spol)[layerchoice2]
   # EQE_ppol2 = tmm.inc_absorp_in_each_layer(front_ppol)[layerchoice2]
    
   # EQEs2.append( (EQE_spol2 + EQE_ppol2) / 2. )
    
    Rfs.append( (front_spol['R']+front_ppol['R']) / 2.)
    Rbs.append( (back_spol['R']+back_ppol['R']) / 2.)
    Ts.append( (front_spol['T']+front_ppol['T']) / 2. )


Ts = np.array(Ts)
Rfs = np.array(Rfs)
Rbs = np.array(Rbs)
As = 1-Ts-Rfs
sanities = Ts+Rfs+As

AbsByAbsorbers = np.array(AbsByAbsorbers)
#EQEs2 = np.array(EQEs2)
#IREQEs=EQEs+EQEs2

# Here I calculate VLT and spit it out to the screen
VLTstack=Stack(layers)
VLT=VLTstack.get_visible_light_transmission(lams,inc_angle)
print("VLT =",VLT)
#

X = np.transpose([lams,AbsByAbsorbers])
np.savetxt('./Output/AbsByAbsorber.txt',X,delimiter=',',header="wavelength [micron], AbsByAbsorber [1]")

Y = np.transpose([lams,Ts,Rfs,Rbs])
np.savetxt('./Output/TRfRb.txt',Y,delimiter=',',header="wavelength [micron], T [1], R_f [1], R_b [1]")

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

# ******************** Here I add PCE calculation *********************#

# Solar cell temperature is 300 kelvin:

Tcell = 300 * K

worksheet = pandas.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
#worksheet = pandas.read_excel('/Users/lwheeler/Code/pv-window-bem/Data/astmg173.xls')
downloaded_array = np.array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
# print(AM15)


# Tack on units
lams *= 1000 * nm

# Tack on the appropriate units:
AM15[:,0] *= nm
AM15[:,1] *= W / m**2 / nm
    
Ephoton = (hPlanck * c0 / lams)
#λ_min = 300 * nm
#λ_max = 2500 * nm
#λ_min = 280 * nm
#λ_max = 4000 * nm
E_min = min(Ephoton) #hPlanck * c0 / λ_max
E_max = max(Ephoton) #hPlanck * c0 / λ_min

# Interpolate to get a continuous function which I will be able to do integrals on:
AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])

# Here’s the plot, it looks correct:

#λs = np.linspace(λ_min, λ_max, num=500)
y_values = np.array([AM15interp(x) for x in lams])
#y_values = np.array([AM15interp(x) for x in λs])
plt.figure()
plt.plot(lams / nm , y_values / (W/m**2/nm))
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
print('Solar Constant = ',solar_constant / (W/m**2))

# I need to get a continuous function of the fraction of the photons absorbed
# in the absorber layer. Here I tack on units and interpolate:
    
# X[:,0] *= 1000 * nm

# AbsInterp = scipy.interpolate.interp1d(X[:,0], X[:,1])

# Tack on units
#lams *= 1000 * nm

# Round AbsByAbsorber to make really small numbers equal to zero
AbsByAbsorbers = AbsByAbsorbers.round(8)
AbsInterp = scipy.interpolate.interp1d(lams, AbsByAbsorbers)
EInterp = scipy.interpolate.interp1d(Ephoton, AbsByAbsorbers)

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

#AbsByAbsorbers

#def AbsVsE(Ephoton):
#    λ = hPlanck * c0 / Ephoton 
#    return AbsInterp(λ)
#return AbsInterp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2)
    
# def EQE(eta,AbsInterp):
#    lam = hPlanck * c0 / Ephoton 
#    return eta * AbsInterp(lam)





Rs = .02 * ohm #lumped series resistance parameter
Rsh = 10 * ohm #shunt resistance
#.01,10000
#0,.2

# Here I input the spectrum of photons absorbed by the absorber material (Absorbed)
# and the electron-hole pair extraction efficiency (eta). EQE = eta * Absorbed

def RR0(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(E) * E**2 / (np.exp(E / (kB * Tcell)) - 1)
    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

def Generated(eta,Absorbed):
    integrand = lambda E : eta * Absorbed(E) * SPhotonsPerTEA(E)
#    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
    return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]


def current_density_ideal(voltage, eta,Absorbed):
    return e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * voltage / (kB * Tcell)))



def current_density(voltage, eta,Absorbed):
    Current_dens = lambda current:(e*Generated(eta,Absorbed) - e*RR0(eta,Absorbed) * np.exp(e * (voltage  + current/(1/m**2) * Rs) / (kB * Tcell))) - ((voltage*(1/m**2) + current * Rs)/Rsh) - (current)
    Solution = scipy.optimize.fsolve(Current_dens, current_density_ideal(voltage, eta, Absorbed), maxfev = 400)
    return Solution

#def RR0(eta,Absorbed):
 #   integrand = lambda E : eta * Absorbed(hPlanck * c0 / E) * E**2 / (np.exp(E / (kB * Tcell)) - 1)
  #  integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
   # return ((2 * np.pi) / (c0**2 * hPlanck**3)) * integral

#def Generated(eta,Absorbed):
 #   integrand = lambda E : eta * Absorbed(hPlanck * c0 / E) * SPhotonsPerTEA(E)
#    integral = scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]
  #  return scipy.integrate.quad(integrand, E_min, E_max, full_output=1)[0]


#def current_density_ideal(voltage, eta,Absorbed):
 #   return e * (Generated(eta,Absorbed) - RR0(eta,Absorbed) * np.exp(e * voltage / (kB * Tcell)))



#def current_density(voltage, eta,Absorbed):
 #   Current_dens = lambda current:(e*Generated(eta,Absorbed) - e*RR0(eta,Absorbed) * np.exp(e * (voltage  + current/(1/m**2) * Rs) / (kB * Tcell))) - ((voltage*(1/m**2) + current * Rs)/Rsh) - (current)
  #  Solution = scipy.optimize.fsolve(Current_dens, current_density_ideal(voltage, eta, Absorbed), maxfev = 400)
   # return Solution
    #Units W/m^2 / V = C/(s*m**2) = A/m**2
                    #        C *    [    V / m**2    +   C *          (A/m**2)  *  ohm) = V/m**2            / J ]             V*C = J
print('Exponent Factor = ', (e * (1 * V  + current_density_ideal(1*V, 0.6, EInterp)/(1/m**2)* Rs) / (kB * Tcell)))# / (V*C/(J)))
#units cancel out
                            # e is in coulombs. This equation^^ results in C*V/(*J) = J/J 
#RR0test = RR0(0.6, AbsInterp)
#print('RR0test =', RR0test / (1/(s*m**2)))
#current = 4 * (A/(m**2))

#print('I = ', current/(A/(m**2)))
print('Shunt resistance factor = ', ((1*V*(1/m**2)  + (current_density(1 * V, 0.6, EInterp)) * Rs)/Rsh) / (A/m**2))
print('RR0 factor = ', e*RR0(0.6,EInterp) * np.exp(e * (.1*V + (current_density(1 * V, 0.6, EInterp))/(1/m**2) * Rs) / (kB * Tcell)) / (A/(m**2)))
#print( 'Generated factor = ', (e*Generated(0.6, AbsInterp) / (A/(m**2))))
print ('RR0 = ',RR0(0.6, EInterp) / (1/(s*m**2)))
print('Gen = ;', Generated(0.6, EInterp) / (1/(s*m**2)))
print('current without resistance (A/m**2) = ', (current_density_ideal(1 * V, 0.6, EInterp)) / (A/(m**2)))

#print('Testnum = ', ((6*V + 6*ohm*A/m**2 / (1/m**2))/(6*ohm)) * (1/m**2) / (A/m**2))
#print('Tesnum2 = ', (e*(V+A*ohm) / (kB * Tcell)))# / (V*C/(J))))
#print('e is units of C?', e/C)
#themoops
#print('current with resistance NEO = ', (current_density_NEO(.1 * V, 0.6, AbsInterp)) / (C/(s*m**2)))
print('current with resistance = ', (current_density(1 * V, 0.6, EInterp)) / (A/m**2))
#moops



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
    #Units V

def J_mpp(eta,Absorbed):
    """ current at max power point """
    return current_density(V_mpp(eta,Absorbed), eta,Absorbed)

def Power(voltage, eta,Absorbed):
    return voltage * current_density(voltage, eta,Absorbed)


Vs = np.linspace(0.1 * V, 2 * V, num=500)
y_values = np.array([Power(x,0.6,EInterp) for x in Vs])
j_values = np.array([current_density(x, 0.6,EInterp) for x in Vs])
plt.figure()
plt.plot(Vs / V , y_values / (W/m**2))
plt.plot(Vs / V , j_values / (A/m**2))
plt.xlabel("Voltage (V)")
plt.ylabel("Power (W m-2)")
plt.ylim(-1, 150)
#plt.title("Light from the sun");
plt.show()

def max_power(eta,Absorbed):
    voltage = V_mpp(eta,Absorbed)
    return voltage * current_density(voltage, eta,Absorbed)
    #Units W/m^2

def max_efficiency(eta,Absorbed):
    return max_power(eta,Absorbed) / solar_constant
    #Units W/m^2 / W/m^2

print("PCE =",max_efficiency(0.6,EInterp))

#print('Vmp = ',V_mpp(0.6,EInterp)/V)
#0.039122219399482544/V == 1.434680175781252

### Stuff from AD
'''
to incorporate nonradiative recombination 
   divide I0 in the IV curve by eta_ext = the external fluorescence efficiency.
   I0 is ideal dark current
see equations 6 from https://aip.scitation.org/doi/10.1063/1.4905277
and 3rd from https://www.pveducation.org/pvcdrom/solar-cell-operation/iv-curve

added eta_ext to def current_density


to incorporate non zero resistance add factor (V-IRs) in place of V in def current_density
to shockley quessier equation as per eq 7 from https://aip.scitation.org/doi/10.1063/1.4905277 
Also add factor (V-IRs)/Rsh to the end. Rsh is shunt resistance,
Rs is lumped series resistance parameter.

Added shunt resistance and series resistance to def curent_density (series resistance :: (voltage = voltage - current*Rs))
Check pveducation shunt resistance, series resistance, for implications


Rch = Vmp/Imp ~ Voc/Isc
Rs = f * Rch
Where f is the fraction power loss from 0 to 1
Power loss is needed to calculate Rs and Rs is needed to calculate powerloss?
    see: pveducation series resistance and characteristic resistance


f=Rch/Rs
f=Rsh/Rch
f/f = rch/rs / rsh/rch
1 = rch^2/(rs*rch)
rch^2 = rs*rch

'''





