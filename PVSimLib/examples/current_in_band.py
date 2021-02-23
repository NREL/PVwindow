from photovoltaic import *
import matplotlib.pyplot as plt

# Reading solar spectra file
wavelength, AM0, AM15G, AM15D = common.solar_spectra()

# Obtaining the current in band
current_in_band = semiconductors.current_in_band(wavelength, AM15G)

plt.figure('Current in band')
plt.plot(wavelength,current_in_band)
# plt.ylim(0, 0.08)
# plt.xlim(200, 1400)
plt.grid(True)
plt.xlabel('wavelength (nm)')       #  add axis labels and plot title
plt.ylabel('Current density in band (mA cm$^{-2}$ nm$^{-1}$)')
plt.title('Current density in band')
plt.savefig('Current_density_in_band.png')
plt.show()