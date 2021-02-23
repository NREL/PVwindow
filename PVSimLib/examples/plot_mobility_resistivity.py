import photovoltaic.semiconductors as pv
import numpy as np
import matplotlib.pyplot as plt


N_D = np.logspace(12, 21, 200)  # sweep the doping in the base (cm-3)

# plot the resisitivity
plt.figure()
plt.plot(N_D, pv.resistivity_Si_n(N_D), label='n-type')
plt.plot(N_D, pv.resistivity_Si_p(N_D), label='p-type')
plt.loglog()
plt.legend()
plt.title('resisitivity in silicon as a function of doping')
plt.xlabel('doping (cm^{-3})')
plt.ylabel('resistivity (ohm cm)')
plt.xlim(1e12, 1e21)
plt.ylim(1e-4, 1e4)
plt.show()
delta_n = pv.implied_carrier(0.7, 1e15)
print(pv.impliedV(1e15, 6.7e15))
G = pv.current2gen(0.038)
tau = delta_n / G
print(delta_n, G, tau)

# plot the mobility models
N_D = np.logspace(14, 22)
mmobility = pv.mob_masetti_phos(N_D)
PC1Dmobility = pv.mob_thurber(N_D, False)
plt.plot(N_D, mmobility, label='Masetti')
plt.plot(N_D, PC1Dmobility, label='Thurber')
plt.xlabel('doping (/cm³)')
plt.ylabel('Mobility (cm²/V·s)')
plt.loglog()
plt.ylim(10, 3000)
plt.xlim(1e14, 1e22)
plt.legend(loc=0)
plt.title('electron mobility in silicon')
plt.show()