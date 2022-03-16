from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def fresnel_Ts(theta,n2):
    n1 = 1.
    num = n1*np.cos(theta) - n2*np.sqrt(1-(n1/n2*np.sin(theta))**2.)
    dem = n1*np.cos(theta) + n2*np.sqrt(1-(n1/n2*np.sin(theta))**2.)
    return 1 - np.abs(num/dem)**2.

def fresnel_Tp(theta,n2):
    n1 = 1.
    num = n1*np.sqrt(1-(n1/n2*np.sin(theta))**2.) - n2*np.cos(theta)
    dem = n1*np.sqrt(1-(n1/n2*np.sin(theta))**2.) + n2*np.cos(theta)
    return 1 - np.abs(num/dem)**2.

def fresnel_T(theta,n2):
    Ts = fresnel_Ts(theta,n2)
    Tp = fresnel_Tp(theta,n2)
    return 0.5*(Ts+Tp)

def fresfit(theta,n2,s):
    return s*fresnel_T(theta,n2)

def schlick(theta,R0,s):
    R = R0 + (1.-R0)*(1. - np.cos(theta))**5.
    return s*R

def dothefits(filename):

    rawdat = np.loadtxt(filename,delimiter=",")

    angles = rawdat[:,0]*np.pi/180
    PCEs = rawdat[:,1]

    popt_s,pcov_s = curve_fit(schlick,angles,PCEs)


    popt_f,pcov_f = curve_fit(fresfit,angles,PCEs)


    return [angles,PCEs,popt_s,pcov_s,popt_f,pcov_f]

fnames = ["PCE_v_angle_50pct_S_IQE0pt8.txt",
          "PCE_v_angle_50pct_S_IQE0pt6.txt",
          "PCE_v_angle_50pct_S_IQE0pt4.txt"]

plt.figure()
for fname in fnames:
    thisname = './Output/' + fname
    [angles,PCEs,popt_s,pcov_s,popt_f,pcov_f] = dothefits(thisname)
    
    print('file: ' + fname)
    print('scale factor (avg PCE kinda): ' + str(popt_f[1]))
    print('effective ref index: ' + str(popt_f[0]))
    print('')
    
    plt.plot(angles,PCEs,linestyle='None',marker='s',label="PCE, "+fname)
    # plt.plot(angles,schlick(angles,*popt_s),label="fit, Schlick")
    plt.plot(angles,fresfit(angles,*popt_f),label="fit, Fresnel")

plt.xlabel("angle, rad")
plt.ylabel("PCE")
plt.legend()
plt.show()

'''
Ts_ex = fresnel_Ts(angles,1.5)
Tp_ex = fresnel_Tp(angles,1.5)
T_ex = fresnel_T(angles,1.5)

plt.figure()
plt.plot(angles,Ts_ex,label="Ts")
plt.plot(angles,Tp_ex,label="Tp")
plt.plot(angles,T_ex,label = "T")
plt.legend()
plt.show()
'''
