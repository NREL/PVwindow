import photovoltaic.common
from photovoltaic import *
import matplotlib.pyplot as plt
T = 300                     # Temperature (K)
E_Si = 11.7                 # Relative electric permittivity of Silicon (No dimension)

# Emitter parameters
We = 0.25e-4        # emitter thickness (cm)
Le = 3.87e-4        # emitter diffusion length (cm)
De = 3              # emitter diffusivity (cm2/s)
Se = 1000           # front surface recomb (cm/s)
thickness_e = 0.25  # emitter thickness (um)

# Base Parameters
Wb = 250e-4         # Base thickness (cm)
Lb = 158.11e-4      # Base diffusion length (cm)
Db = 5              # Base diffusivity (cm2/s)
Sb = 2000           # Back surface recomb (cm/s)
thickness_b = 250   # Base thickness (um)

# Concentration parameters
ni = semiconductors.ni_Si(T)   # Intrinsic carrier concentration in silicon at 300K (cm^-3)
Nd = 2e19           # Majority carrier concentration on N side
Na = 3e16           # Majority carrier concentration on P side

# Reading solar spectra file
wavelength, AM0, AM15G, AM15D = common.solar_spectra()

# Calculate depletion region and built-in voltage
V_0 = semiconductors.built_in_voltage(Nd, Na, ni, T)
print("Built-in voltage: " + str(V_0) + " V")

Wd = semiconductors.depletion_region_width(Nd, Na, E_Si, V_0)
print("Depletion region width: " + str(Wd*10000) + " um")

wavelength_2, abs_coeff = common.read_abs_coefficient()
abs_coeff = photovoltaic.common.adapt_vector_to_wavelength(wavelength, wavelength_2, abs_coeff)
Wd = semiconductors.depletion_region_width(Nd, Na, E_Si, V_0)
QE, EQE, BQE, DQE = semiconductors.quantum_efficiencies(abs_coeff, Wd, Se, Le, De, We, Sb, Wb, Lb, Db)

plt.figure('Quantum Efficiency by Regions')
plt.plot(wavelength, EQE, label="Emitter")
plt.plot(wavelength, DQE, label="Depletion region")
plt.plot(wavelength, BQE, label="Base")
plt.plot(wavelength, QE, label="Total")
plt.legend(loc='upper right')
plt.ylim(0, 1.05)
plt.xlim(200, 1400)
plt.grid(True)
plt.xlabel('wavelength (nm)')       #  add axis labels and plot title
plt.ylabel('Quantum Efficiency')
plt.title('Quantum Efficiency by Regions')
plt.savefig('Quantum_Efficiency_by_Region.png')

plt.show()

wv2, absorp = common.read_abs_coefficient('absorption_wv.txt')

QE, EQE, BQE, DQE = semiconductors.quantum_efficiencies(absorp, Wd, Se, Le, De, We, Sb, Wb, Lb, Db)

plt.figure('Quantum Efficiency by Regions')
plt.plot(wavelength, EQE, label="Emitter")
plt.plot(wavelength, DQE, label="Depletion region")
plt.plot(wavelength, BQE, label="Base")
plt.plot(wavelength, QE, label="Total")
plt.legend(loc='upper right')
plt.ylim(0, 1.05)
plt.xlim(200, 1400)
plt.grid(True)
plt.xlabel('wavelength (nm)')       #  add axis labels and plot title
plt.ylabel('Quantum Efficiency')
plt.title('Quantum Efficiency by Regions')
plt.savefig('Quantum_Efficiency_by_Region_2.png')
plt.show()