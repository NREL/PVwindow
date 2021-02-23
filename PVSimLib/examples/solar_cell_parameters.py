import photovoltaic.common
from photovoltaic import *
import numpy as np
# General Parameters
T = 300                     # Temperature (K)
E_Si = 11.7                 # Relative electric permittivity of Silicon (No dimension)
Eg = semiconductors.Eg0_paessler(T)    # Band-gap of Silicon (eV)
W = 15.6                    # Solar cell width (cm)
Z = 15.6                    # Solar cell height (cm)

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
# Reading solar spectra file
wavelength, AM0, AM15G, AM15D = common.solar_spectra()

# Obtaining the current in band
current_in_band = semiconductors.current_in_band(wavelength, AM15G)

# Calculate depletion region and built-in voltage
V_0 = semiconductors.built_in_voltage(Nd, Na, ni, T)
print("Built-in voltage:  {:.2f} V".format(V_0))

Wd = semiconductors.depletion_region_width(Nd, Na, E_Si, V_0)
print("Depletion region width: {:.2f} um".format(Wd*10000))

# Reading absorption coefficient and calculation QE for all the regions
wavelength_2, abs_coeff = common.read_abs_coefficient()
abs_coeff = photovoltaic.common.adapt_vector_to_wavelength(wavelength, wavelength_2, abs_coeff)

QE, EQE, BQE, DQE = semiconductors.quantum_efficiencies(abs_coeff, Wd, Se, Le, De, We, Sb, Wb, Lb, Db)

# Total J_sc
print('\n')
print("Max photocurrent: {:.2f} mA/cm^2".format(semiconductors.photocurrent(wavelength, current_in_band, Eg)))

J_SC = solarcells.short_circuit_current(current_in_band, QE)
print("J_SC: {:.2f} mA/cm^2".format(J_SC))


# Total J_0
J_0_E = semiconductors.J0_layer(We, Nd, De, Le, Se, ni)*1000
J_0_B = semiconductors.J0_layer(Wb, Na, Db, Lb, Sb, ni)*1000
J_0 = J_0_E + J_0_B
print('\n')
print("J_0(Narrow Base): {:.2e} mA/cm^2".format(J_0))

# By regions
print("By regions: Emitter and Base")
print("J_0 Emitter: {:.2e} mA/cm^2".format(J_0_E))
print("J_0 Base: {:.2e} mA/cm^2".format(J_0_B))

# Voc
V_OC = solarcells.Voc(J_SC, J_0, 1, T)
print('\n')
print("Voc: {:.2f} V".format(V_OC))

# Fill factor
FF = solarcells.FF_ideal(V_OC, 1, T)
print("FF: {:.2f}".format(FF))

# Jmp and Vmp
Jmp = J_SC*np.sqrt(FF)
Vmp = V_OC*np.sqrt(FF)
print("Jmp: {:.2f} mA/cm^2".format(Jmp))
print('Vmp: {:.2f} V'.format(Vmp))

# Efficiency
efficiency = solarcells.efficiency(V_OC, J_SC/1000, FF, 1)
print('Efficiency: {:.2f} %'.format(efficiency))
