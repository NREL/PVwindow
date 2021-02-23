import numpy as np
import os
import photovoltaic.constants as const


def sind(angle):
    """
    Return the sine of an angle
    :param angle: Angle(Degrees)
    :return: Sine of the angle
    """
    return np.sin(np.radians(angle))


def cosd(angle):
    """
    Return the cosine of an angle
    :param angle: Angle(Degrees)
    :return: Cosine of the angle
    """
    return np.cos(np.radians(angle))


def tand(angle):
    """
    Return the tangent of an angle
    :param angle: Angle(Degrees)
    :return: Cosine of the angle
    """
    return np.tan(np.radians(angle))


def arcsind(x):
    """
    Return the arcsine of an angle
    :param x:
    :return: arcsine of the angle (Degrees)
    """
    return np.degrees(np.arcsin(x))


def arccosd(x):
    """
    Return the arccosine of an angle
    :param x:
    :return: arcsine of the angle (Degrees)
    """
    return np.degrees(np.arccos(x))


def nm2eV(x):
    """
    Given wavelength (nm) of a photon return the energy (eV)
    :param x: Wavelength of the photon (nm)
    :return: Energy (eV)
    """
    return const.hc_q * 1e9 / x


def eV2nm(x):
    """
    Given energy (eV) of a photon return the wavelength (nm)
    :param x: Energy (eV)
    :return: Wavelength of the photon (nm)
    """
    return const.hc_q * 1e9 / x


def nm2joule(x):
    """
    Given wavelength (nm) of a photon return the energy (J)
    :param x: Wavelength of the photon (nm)
    :return: Energy (J)
    """
    return const.h * const.c * 1e9 / x


def Vt(T=298.15):
    """
    Return thermal voltage (volts) at given temperature, T(Kelvin).
    The default temperature is 298.15 K, which is equal to 25 °C
    :param T: Temperature (K)
    :return: Thermal voltage (V)
    """
    return const.k * T / const.q


def probability_fermi_dirac(E, Ef, T):
    """
    Return the fermi dirac function (units) where E is the energy (), Ef is the fermi
    given the energies in electron volts
    :param E: Energy (eV)
    :param Ef: Fermi level (eV)
    :param T: Temperature (K)
    :return: Fermi-Dirac function.
    """
    return 1 / (np.exp((E - Ef) / Vt(T)) + 1.0)


def probability_maxwell_boltzmann(E, Ef, T):
    """
    Return the Maxwell-Blotzmann function where E is the energy (), Ef is the fermi
    given the energies in electron volts
    :param E: Energy (eV)
    :param Ef: Fermi level (eV)
    :param T: Temperature (K)
    :return: Maxwell-Blotzmann function
    """
    return 1 / (np.exp((E - Ef) / Vt(T)))


def probability_bose_einstein(E, Ef, T):
    """
    Return the Bose-Einstein function where E is the energy (), Ef is the fermi
    given the energies in electron volts
    :param E: Energy (eV)
    :param Ef: Fermi level (eV)
    :param T: Temperature (K)
    :return: Bose-Einstein function
    """
    return 1 / (np.exp((E - Ef) / Vt(T)) - 1.0)

def read_abs_coefficient(fname=None):
    """
    Returns an array with the absorbtion coefficient of a material
    :param fname: Name of the file (String). Eg: 'files/Si_abs_coeff.txt'
    :return: wavelngth (nm), absorption coefficient (cm-1)
    """
    if fname is None:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/Si_abs_coeff.txt')
    else:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/' + fname)
    wavelength, abs_coeff= np.loadtxt(fname, skiprows=1, unpack=True)
    return wavelength, abs_coeff


def read_nk(fname=None):
    """
    Return refractive indexes.
    :param fname: Name of the file (String). Eg: 'files/refractive_index.csv'
    :return: real refractive index, imaginary refractive index
    """

    if fname is None:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/Si_refractive_index.csv')
    else:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/' + fname)
    wv, nd, kd = np.loadtxt(fname, skiprows=1, delimiter=',', unpack=True)
    return wv, nd, kd


def solar_spectra(fname=None):
    """Return wavelength (nm) and AM0, AM15G, AM15D (W/m²/nm)
    of the standard spectrum.
    reference XXX, DOI XXX"""
    if fname is None:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/ASTMG173.txt')
    else:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'files/' + fname)
    wavelength, AM0, AM15G, AM15D = np.genfromtxt(fname, skip_header=2, unpack=True)
    return wavelength, AM0, AM15G, AM15D

def station_info(fname):
    """
    Read TMY data file, import by station number
    :param fname:
    :return: city, station, GMT_offset, latitude, longitude, altitude, ETR(Wm^-2), GHI(Wm^-2), DNI(Wm^-2), DHI(Wm^-2),
     ambient_temperature(C)
    """

    package_path = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(package_path, 'files/'+fname)
    city = np.genfromtxt(fname, dtype='U', max_rows=1, delimiter=",", usecols=1)
    station, GMT_offset, latitude, longitude, altitude = np.genfromtxt(fname, max_rows=1, delimiter=",",
                                                                       usecols=(0, 3, 4, 5, 6))
    ETR, GHI, DNI, DHI, ambient_temperature = np.genfromtxt(fname, skip_header=2, delimiter=",",
                                                            usecols=(2, 4, 7, 10, 31), unpack=True)
    return city, station, GMT_offset, latitude, longitude, altitude, ETR, GHI, DNI, DHI, ambient_temperature


def adapt_vector_to_wavelength(x1, x2, y):
    """Adapts vector y to the length and step of vector x1. vector x2 is the current key of vector y"""
    i = 0
    j = 0
    z = 0
    final_coeff = np.zeros(x1.shape[0])
    last = y[0]
    temp_vec = []
    while (i <= (x1.shape[0]-1)) and (j <= x2.shape[0]-1):
        if x1[i] == x2[j]:
            if z > 0 :
                increment = (y[j] - last)/z
                while z > 0 :
                    final_coeff[i-z] = y[j] - z*increment
                    z -= 1
            final_coeff[i] = y[j]
            last = y[j]
            i += 1
            j += 1
        elif x1[i] > x2[j]:
            j += 1
        elif x1[i] < x2[j]:
            z += 1
            i += 1
    return final_coeff
