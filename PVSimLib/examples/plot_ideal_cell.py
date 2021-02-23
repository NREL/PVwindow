import matplotlib.pyplot as plt
import numpy as np
from photovoltaic import *
T = 300
# Emitter parameters
We = 0.25e-4        # emitter thickness (cm)
Le = 3.87e-4        # emitter diffusion length (cm)
De = 3              # emitter diffusivity (cm2/s)
Se = 1000           # front surface recomb (cm/s)

# Base Parameters
Wb = 250e-4         # Base thickness (cm)
Lb = 158.11e-4      # Base diffusion length (cm)
Db = 5              # Base diffusivity (cm2/s)
Sb = 2000           # Back surface recomb (cm/s)

# Concentration parameters
ni = semiconductors.ni_Si(T)   # Intrinsic carrier concentration in silicon at 300K (cm^-3)
Nd = 2e19           # Majority carrier concentration on N side
Na = 3e16           # Majority carrier concentration on P side

IL = 0.03642
I0 = semiconductors.J0_layer(We, Nd, De, Le, Se, ni) + semiconductors.J0_layer(Wb, Na, Db, Lb, Sb, ni)

V = np.linspace(0, 0.7, 100)
I = solarcells.I_cell(V, IL, I0)
plt.plot(V, I, color="blue")
plt.xlim(0, 0.8)
plt.ylim(0, 0.04)
plt.xlabel('voltage (V)')
plt.ylabel('current density (A cm$^{-2}$)')
plt.grid(True)
plt.plot(V, I, label="Ideal Cell")
plt.legend(loc='lower left')
plt.savefig('plot_ideal_cell.png')
plt.show()
