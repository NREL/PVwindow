import photovoltaic.common as common
import photovoltaic.constants as const
import numpy as np

from scipy import integrate


def am_intensity(airmass):
    """Return radiation intensity (W/m**2) given airmass (units) """
    It = 1.353 * (0.7 ** (airmass ** 0.678))
    Ig = It * 1.1
    return It, Ig


def air_mass(angle):
    """Return air mass (units) where *angle* is the zenith angle (degrees) """
    return 1 / common.cosd(angle)


def am_shadow(s, h):
    """Return the air mass (units) where h is the height of a pole and s is the length of its shadow.
     s and h are the same length units, i.e. both in m, ft, cm etc."""
    am = np.sqrt(1 + (s / h) ** 2)
    return am


def blackbody_peak(T):
    """Return the peak wavelength (nm) of a black body spectrum where T is the temperature (K)"""
    return 1e9 * const.Wien / T


def blackbody_integrated(T):
    """Return integrated irradiance (W/m2/steradian) from a blackbody where T is the temperature (K)"""
    return const.Stefan_Boltzmann * T ** 4


def blackbody_spectrum(wavelength, T=6000):
    """Return the blackbody irradaiance (W/nm) at a given wavelength (nm) and temperature, T (K). """
    wavelength = wavelength * 1e-9
    F = 2 * const.π * const.h * const.c ** 2 / ((wavelength ** 5) * (np.exp(const.h * const.c /
                                                                            (wavelength * T * const.k)) - 1))
    return F * 1e-9  # convert to W/m to W/nm


def equal_spacing(x, y, x_min, x_max, x_step):
    """Returns spectra with equal spacking and truncation (W/m2) (NOT W/m2/nm)
    given a spectrum (W/m2/nm) as a function of wavelength (nm)
    wavelength minimum (nm) and wavlength maximum (nm)
    Note: would actually work for any x y data"""
    y_integrated = integrate.cumtrapz(y, x, initial=0)
    x_midpoint = np.arange(x_min, x_max + x_step, x_step)
    x_extents = np.arange(x_min - x_step / 2, x_max + 3 / 2 * x_step, x_step)
    y_spaced = np.diff(np.interp(x_extents, x, y_integrated))
    return x_midpoint, y_spaced


def space_solar_power(x):
    """Return the radiant power density (W/m²) where x is the distance from  the sun (m)"""
    return 2.8942e25 / (x ** 2)


def etr_earth(day_no):
    """Return extraterrestrial radiation at earth (W/m**2) where day_no is the day of the year (day).
     January 1 has a day_no of 1. There is no correction for leap year."""
    return (1 + .033 * (common.cosd((day_no - 2.0) * (360.0 / 365.0))) * 1353)


def declination(day_no):
    """Return declination angle of sun (degrees) where day_no is the day of the year (day).
    For Jan 1 day_no = 1, Dec 31 dayno = 365. There is no correction for leap years"""
    return 23.45 * common.sind((day_no - 81) * (360 / 365))


def equation_of_time(day_no):
    """Return the equation of time (minutes) where day_no is the day of the year (day). """
    B = 360.0 / 365.0 * (day_no - 81.0)
    EoT = 9.87 * common.sind(2 * B) - 7.53 * common.cosd(B) - 1.5 * common.sind(B)
    return EoT


def time_correction(EoT, longitude, GMTOffset):
    """Return the time correction (minutes) where EoT is the equation of time,
    and given location longitude (degrees) and the GMT offset (hours)"""
    LSTM = 15.0 * GMTOffset
    TimeCorrection = 4.0 * (longitude - LSTM) + EoT
    return TimeCorrection


def elevation(declination, latitude, local_solar_time):
    """Return the elevation angle of the sun (degrees)
    given declination (degrees), latitude (degrees) and local_solar_time (hours) """
    hra = 15.0 * (local_solar_time - 12.0)
    return common.arcsind(common.sind(declination) * common.sind(latitude) +
                          common.cosd(declination) * common.cosd(latitude) * common.cosd(hra))


def sun_rise_set(latitude, declination, time_correction):
    """Return the sunrise and sunset times in hours
    given the latitude (degrees) and the declination (degrees)
    """
    A = -1 * (common.sind(latitude) * common.sind(declination)) / (common.cosd(latitude) * common.cosd(declination))
    local_solar_time = common.arccosd(A) / 15.0
    sunrise = 12.0 - local_solar_time - (time_correction / 60.0)
    sunset = 12 + local_solar_time - (time_correction / 60.0)
    return sunrise, sunset


def elev_azi(declination, latitude, local_solar_time):
    """Return the elevation (degrees) and azimuth (degrees)"""
    hour_angle = 15.0 * (local_solar_time - 12.0)
    elevation = common.arcsind(common.sind(declination) * common.sind(latitude) +
                               common.cosd(declination) * common.cosd(latitude) * common.cosd(hour_angle))
    azimuth = common.arccosd((common.cosd(latitude) * common.sind(declination) -
                              common.cosd(declination) * common.sind(latitude) * common.cosd(hour_angle)) /
                             common.cosd(elevation))
    # the multiplication by 1.0 causes a single value return for single inputs, otherwise it returns an array of one
    # element
    azimuth = np.where(hour_angle > 0, 360.0 - azimuth, azimuth) * 1.0
    return elevation, azimuth


def module_direct(azimuth, elevation, module_azimuth, module_tilt):
    """Returns the faction of light on a arbtrarily tilted surface
     given sun azimuth (degrees) where north is zero and elevation
     module_azimuth and module_tilt, where """
    fraction = common.cosd(elevation) * common.sind(module_tilt) * common.cosd(azimuth - module_azimuth) + \
               common.sind(elevation) * common.cosd(module_tilt)
    return fraction


def sun_position(dayNo, latitude, longitude, GMTOffset, H, M):
    """Return the position of the sun as a elevation and azimuth given
    latitude, logitude and the GMTOffset, """
    EoT = equation_of_time(dayNo)
    TimeCorrection = time_correction(EoT, longitude, GMTOffset)
    local_solar_time = H + (TimeCorrection + M) / 60.0
    elevation, azimuth = elev_azi(declination(dayNo), latitude, local_solar_time)
    return elevation, azimuth


def spectrum_spacing(x, y, x_min, x_max, x_step):
    """Returns spectra in W/m2 (NOT W/m2/nm) with equal spacing
    given a spectrum (W/m2/nm) as a function of wavelength (nm)
    wavelength minimum (nm) and wavlength maximum (nm) """
    y_integrated = integrate.cumtrapz(y, x, initial=0)
    x_midpoint = np.arange(x_min, x_max + x_step, x_step)
    x_extents = np.arange(x_min - x_step / 2, x_max + 3 / 2 * x_step, x_step)
    y_spaced = np.diff(np.interp(x_extents, x, y_integrated))
    return x_midpoint, y_spaced


def photon_flux(power, wavelength):
    """Return the photon flux (/s) given the power of light (watts) and wavelength (nm)
    If power is in W/m2 then flux is in m-2s-1"""
    return power / common.nm2joule(wavelength)


def k2alpha(kd, wavelength):
    """Quick convert of extinction coefficient (units) to absorption coefficient(cm-1) given the wavelength (nm)"""
    return 1e7 * 4 * const.π * kd / wavelength


def refraction(n1, n2, θ1):
    """Return the refracted angle of light (degrees) where n1 is the refractive index of incident medium (units),
    n2 is the refractive index of the transmission medium (units) and θ1 is the incident angle to the normal"""
    θ2 = common.arcsind(n1 / n2 * common.sind(θ1))
    return θ2


def absorption_coeff(kd, wavelength):
    """absorption coefficient (cm-1) from extinction coefficient (units) and wavelength (nm)
    wavelength """
    return 1e7 * 4 * const.π * kd / wavelength


def transmittance(abs_coeff, thickness):
    """Return the fraction of light transmitted (units) where abs_coeff is the absorption coefficient (cm-1)
    and 'thickness' is the depth in the material (cm)
    """
    return np.exp(-abs_coeff * thickness)


def ARC_thick(wavelength, n1):
    """Return optimal anti-reflection coating thickness (nm) at a given wavelength (nm)
    where n1 is the refractive index of the antireflection coating.
    The returned unit is in the same as the wavelength unit."""
    return wavelength / (4 * n1)


def ARC_opt_n(n0, n2):
    """Return the optimal refractive index, typically denoted as n1, for an antireflection coating (units)
    where n0 is the refractive index of the incident medium and n2 is the refractive index of the object (units)"""
    return np.sqrt(n0 * n2)


def ARC_refl(wavelength, n0, n1, nSemi, thickness):
    """Return the reflectivity from a object (units) that has an anti-reflection coating.
    Where:
    n0 - the ambient refractive index units),
    n1 - refractive index of the dielectric layer (units)
    nSemi - refractive index of the semiconductor (units)
    The reflectivity is in the range of 0 to 1.
    """
    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - nSemi) / (n1 + nSemi)
    θ = (2 * const.π * n1 * thickness) / wavelength
    reflectivity = 100 * (r1 * r1 + r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ)) / (
        1 + r1 * r1 * r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ))
    return reflectivity


def DLARC_refl(wavelength, n0, n1, n2, nSemi, thickness1, thickness2):
    """Return the reflectivity from a object (units) that has a double layer anti-reflection coating, DLARC.
    Where:
    n0 - refractive index of the ambient (units)
    n1 - refractive index of the dielectric layer 1 (units)
    n2 - refractive index of the dielectric layer 2 (units)
    nSemi - refractive index of the semiconductor (units)
    wavelength, thickness1, thickness 2 all in same units (m) or (nm) etc.
    """
    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - n2) / (n1 + n2)
    r3 = (n2 - nSemi) / (n2 + nSemi)
    θ1 = (2 * const.π * n1 * thickness1) / wavelength
    θ2 = (2 * const.π * n2 * thickness2) / wavelength

    numerator = r1 * r1 + r2 * r2 + r3 * r3 + r1 * r1 * r2 * r2 * r3 * r3 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))
    denominator = 1 + r1 * r1 * r2 * r2 + r1 * r1 * r3 * r3 + r3 * r3 * r2 * r2 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))

    return numerator / denominator
