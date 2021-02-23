import numpy as np

q = 1.60217662e-19  # (coulombs) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
k_eV = 8.6173303e-05  # (eV K^-1)
Wien = 2.898e-3  # (m K)
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
π = np.pi  # yes, I use unicode
pi = np.pi  # compatibility with class
h = 6.62607004e-34  # (J.s)
hbar = 6.62607004e-34 / (2 * π)  # usable
c = 299792458.0  # (m s^-1)
hc_q = h * c / q  # 1.2398419745831506e-06
ε = 8.85418782e-12  # Vacuum permittivity (s^4 A^2 m^-3 kg^-1)
ε_Si = 11.7
