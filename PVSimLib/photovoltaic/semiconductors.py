#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Basic photovoltaic functions
Provides the formulas and equations typically used in an introductory photocoltaics textbook.

Typical solar units are used, NOT SI units. The units are denoted in parenthesis on the comment lines.
wavelength (nm)
Energy of  photon (eV)
semiconductor dimensions (cm)
degrees instead of radians.
Temperature of 298.15 K (25 degC) not 300 K

The first line on all input files is ignored to allow for column headers
# denotes a comment in input files and is ignored.

Contributions by: sgbowden, richter, heraimenka, jhul etc

Variables and acronyms
ARC - anti-reflection coating
"""
__version__ = '0.1.3'

import numpy as np
import photovoltaic.constants as const
import photovoltaic.common as common


def lifetime_f_length(L, D):
    """
    Return the lifetime (s) given the  diffusion length (cm) and diffusivity (cm2/s).
    :param L: Diffusion length (cm)
    :param D: Diffusivity (cm^2 s^-1)
    :return: Lifetime (s)
    """
    return L**2 / D


def diffusion_length(lifetime, D):
    """
    Return carrier diffusion length (cm) given carrier lifetime(s) and diffusivity (cm2/s).
    :param lifetime: Lifetime (s)
    :param D: Diffusivity (cm^2 s^-1)
    :return: Diffusion length (cm))
    """
    return np.sqrt(lifetime * D)


def tau_b__tau_eff(tau_eff, S, thickness):
    """
    Return the bulk lifetime (s) Given tau_eff (s) surface recombination (cm/s) thickness (cm)
    :param tau_eff: Effective lifetime (s)
    :param S: Surface recombination(cm s^-1)
    :param thickness: Thickness (cm)
    :return: bulk lifetime (s)
    """
    return tau_eff - thickness / (2 * S)


def diffusivity(mobility, T=298.15):
    """
    Return the diffusivity (cm²/s) given the mobility (cm²/Vs). This is also known as the Einstein relation.
    :param mobility: Carrier mobility (cm^2 V^-1 s^-1)
    :param T: Temperature (K)
    :return: Diffusivity (cm^2 s^-1)
    """
    return mobility * common.Vt(T)


def mobility(D, T=298.15):
    """
    Return the mobility of carriers (cm²/Vs) where D is the diffusivity (cm²/s). This is also known as the Einstein
    relation
    :param D: Diffusivity (cm^2 s^-1)
    :param T: Temperature (K)
    :return: Carrier mobility (cm^2 V^-1 s^-1)
    """
    return D / common.Vt(T)


def equilibrium_carrier(N, ni=8.6e9):
    """
    Return the majority and minority carrier concentrations (cm-3) of a semiconductor at equilibrium
    where N is the doping level (cm-3) and ni is the intrinsic carrier concentratoin (cm-3)
    Strictly N and ni just have to be in the same units but (cm-3 is almost always used.
    :param N: Doping level (cm-3)
    :param ni: Intrinsic carrier concentration (cm-3)
    :return: Majority and minority carrier concentrations (cm-3)
    """
    majority = N
    minority = N / (ni ** 2)
    return majority, minority


def conductivity(n, p, ue, uh):
    """
    Return the conductivity of a material(siemens)
    :param n: concentration of electrons (cm-3)
    :param p: concentration of holes (cm-3)
    :param ue: electron mobility (cm²/Vs)
    :param uh: hole mobility (cm²/Vs)
    :return: Conductivity (Siemens)
    """
    return const.q * ue * n + const.q * uh * p


def resistivity_Si_n(Ndonor):
    """
    Return the resistivity of n-type silicon (ohm cm) given the doping of donors(cm-3)
    :param Ndonor: Doping of donors(cm-3)
    :return: Resistivity of n-type silicon (ohm cm)
    """
    n_minority = ni_Si() ** 2 / Ndonor
    return 1 / ((const.q * mob_thurber(Ndonor, False) * Ndonor) +
                (const.q * mob_thurber(n_minority, False, False) * n_minority))


def resistivity_Si_p(Nacceptor):
    """
    Return the resistivity of p-type silicon (ohm cm) given the doping of acceptors(cm-3).
    :param Nacceptor: Doping of acceptors(cm-3)
    :return: Resistivity of p-type silicon (ohm cm)
    """
    n_minority = ni_Si() ** 2 / Nacceptor
    return 1 / ((const.q * mob_thurber(Nacceptor) * Nacceptor) +
                (const.q * mob_thurber(n_minority, True, False) * n_minority))


def mob_thurber(N, p_type=True, majority=True):
    """
    Return the mobility of carriers in silicon according to the model of Thurbur as a function of doping. P_type is
    True or 1 for p doped material and False or 0 for n-type. Majority is True or 1 for majority carriers and False or
     0 for minority carriers.
    :param N: Doping level (cm-3)
    :param p_type: P-type(True or 1) or N-type(False or 0)
    :param majority: Majority carriers(True or 1) or minority carriers(False or 0)
    :return: Mobility of carriers in silicon according to the model of Thurber as a function of doping(cm^2 V^-1 s^-1)
    """
    i = 2 * p_type + majority
    # n-type minority, n-type majority, p-type minority, p-type majority
    umax = [1417, 1417, 470, 470][i]
    umin = [160, 60, 155, 37.4][i]
    Nref = [5.6e16, 9.64E16, 1e17, 2.82E17][i]
    a = [0.647, 0.664, 0.9, 0.642][i]
    return umin + (umax - umin) / (1 + ((N / Nref) ** a))


def mob_masetti_phos(N):
    """
    Return the mobility of carriers (cm²/Vs) in phosphorus doped silicon accoording to the model of Masetti1983.
    Where N is the doping (cm-3). http://dx.doi.org/10.1109%2FT-ED.1983.21207
    :param N: Doping (cm^-3)
    :return: mobility of carriers (cm^-2 V^-1 s^-1)
    """
    µmax = 1414
    µmin = 68.5
    u1 = 56.1
    Cr = 9.20e16
    Cs = 3.41e20
    a = 0.711
    b = 1.98
    return µmin + (µmax - µmin) / (1 + ((N / Cr) ** a)) - u1 / (1 + ((Cs / N) ** b))


def mob_klassen(Nd, Na, Δn=1, T=298.16):
    """
    Return the mobility (cm2/Vs) according to Klassen model.
    :param Nd: Donor Doping (cm^-3)
    :param Na: Acceptors Doping (cm^-3)
    :param Δn: Excess of carriers (cm^-3)
    :param T: Temperature (K)
    :return: electrons and holes mobility (cm^2 V^-1 s^-1)
    """
    s1 = 0.89233
    s2 = 0.41372
    s3 = 0.19778
    s4 = 0.28227
    s5 = 0.005978
    s6 = 1.80618
    s7 = 0.72169
    r1 = 0.7643
    r2 = 2.2999
    r3 = 6.5502
    r4 = 2.367
    r5 = -0.01552
    r6 = 0.6478
    fCW = 2.459
    fBH = 3.828
    mh_me = 1.258
    me_m0 = 1

    T = 298.16
    n0, p0 = equilibrium_carrier(Nd)

    n_i = 8.31E+09

    cA = 0.5
    cD = 0.21
    Nref_A = 7.20E+20
    Nref_D = 4.00E+20

    p = p0 + Δn
    n = n0 + Δn
    cc = p + n

    Za_Na = 1 + 1 / (cA + (Nref_A / Na) ** 2)
    Zd_Nd = 1 + 1 / (cD + (Nref_D / Nd) ** 2)

    Na_h = Za_Na * Na
    Nd_h = Zd_Nd * Nd

    phosphorus = [1414, 68.5, 9.20E16, 0.711, 2.285]
    boron = [470.5, 44.9, 2.23E+17, 0.719, 2.247]

    boron_µmax = 470.5
    boron_µmin = 44.9
    boron_Nref_1 = 2.23E+17
    boron_α = 0.719
    boron_θ = 2.247

    phosphorus_µmax = 1414
    phosphorus_µmin = 68.5
    phosphorus_Nref_1 = 9.20E16
    phosphorus_α = 0.711
    phosphorus_θ = 2.285

    µ_eN = phosphorus_µmax ** 2 / (phosphorus_µmax - phosphorus_µmin) * (T / 300) ** (3 * phosphorus_α - 1.5)
    µ_hN = boron_µmax ** 2 / (boron_µmax - boron_µmin) * (T / 300) ** (3 * boron_α - 1.5)

    µ_ec = phosphorus_µmax * phosphorus_µmin / (phosphorus_µmax - phosphorus_µmin) * (300 / T) ** 0.5
    µ_hc = boron_µmax * boron_µmin / (boron_µmax - boron_µmin) * (300 / T) ** 0.5

    µ_eD = µ_eN * (phosphorus_Nref_1 / Nd) ** phosphorus_α + µ_ec * (cc / Nd)
    µ_hA = µ_hN * (boron_Nref_1 / Na) ** boron_α + µ_hc * (cc / Na)

    Ne_sc = Na_h + Nd_h + p
    Nh_sc = Na_h + Nd_h + n

    PBHe = 1.36e+20 / cc * me_m0 * (T / 300) ** 2
    PBHh = 1.36e+20 / cc * mh_me * (T / 300) ** 2

    PCWe = 3.97e+13 * (1 / (Zd_Nd ** 3 * (Nd_h + Na_h + p)) * ((T / 300) ** 3)) ** (2 / 3)
    PCWh = 3.97e+13 * (1 / (Za_Na ** 3 * (Nd_h + Na_h + n)) * ((T / 300) ** 3)) ** (2 / 3)

    Pe = 1 / (fCW / PCWe + fBH / PBHe)
    Ph = 1 / (fCW / PCWh + fBH / PBHh)

    G_Pe = 1 - s1 / ((s2 + (1 / me_m0 * 300 / T) ** s4 * Pe) ** s3) + s5 / (((me_m0 * 300 / T) ** s7 * Pe) ** s6)
    G_Ph = 1 - s1 / ((s2 + (1 / (me_m0 * mh_me) * T / 300) ** s4 * Ph) ** s3) + s5 / (
        ((me_m0 * mh_me * 300 / T) ** s7 * Ph) ** s6)

    F_Pe = (r1 * Pe ** r6 + r2 + r3 / mh_me) / (Pe ** r6 + r4 + r5 / mh_me)
    F_Ph = (r1 * Ph ** r6 + r2 + r3 * mh_me) / (Ph ** r6 + r4 + r5 * mh_me)

    Ne_sc_eff = Nd_h + G_Pe * Na_h + p / F_Pe
    Nh_sc_eff = Na_h + G_Ph * Nd_h + n / F_Ph

    # Lattice Scattering
    µ_eL = phosphorus_µmax * (300 / T) ** phosphorus_θ
    µ_hL = boron_µmax * (300 / T) ** boron_θ

    µe_Dah = µ_eN * Ne_sc / Ne_sc_eff * (phosphorus_Nref_1 / Ne_sc) ** phosphorus_α + µ_ec * ((p + n) / Ne_sc_eff)
    µh_Dae = µ_hN * Nh_sc / Nh_sc_eff * (boron_Nref_1 / Nh_sc) ** boron_α + µ_hc * ((p + n) / Nh_sc_eff)

    µe = 1 / (1 / µ_eL + 1 / µe_Dah)
    µh = 1 / (1 / µ_hL + 1 / µh_Dae)

    return µe, µh


def Eg0_paessler(T=298.15):
    """
    Return the bandgap of silicon (eV) according to Paessler2002, where T is the temperature (K).
    Code adapted from Richter Fraunhofer ISE. https://doi.org/10.1103/PhysRevB.66.085201
    :param T: Temperature (K)
    :return: Band-gap (eV)
    """
    # constants from Table I on page 085201-7
    α = 3.23 * 0.0001  # (eV/K)
    Θ = 446  # (K)
    Δ = 0.51
    Eg0_T0 = 1.17  # eV     band gap of Si at 0 K

    Tdelta = 2 * T / Θ
    wurzel = (1 + const.π ** 2 / (3 * (1 + Δ ** 2)) * Tdelta ** 2 + (
        3 * Δ ** 2 - 1) / 4 * Tdelta ** 3 + 8 / 3 * Tdelta ** 4 + Tdelta ** 6) ** (1 / 6)
    Eg0 = Eg0_T0 - α * Θ * ((1 - 3 * Δ ** 2) / (np.exp(Θ / T) - 1) + 3 / 2 * Δ ** 2 * (wurzel - 1))
    return Eg0


def ni_Si(T=298.15):
    """
    Return the intrinsic carrier concentration of silicon (cm**-3) according to Sproul94, where T is the
    temperature (K). http://dx.doi.org/10.1063/1.357521
    :param T: Temperature (K)
    :return:  intrinsic carrier concentration of silicon (cm^-3)
    """
    return 9.38e19 * (T / 300) * (T / 300) * np.exp(-6884 / T)


def ni_misiakos(T=298.15):
    """
    Return the intrinsic carrier concentration (cm-3) without band gap narrowing according to Misiakos,
    where T is the temperature (K). DOI http://dx.doi.org/10.1063/1.354551
    :param T: Temperature (K)
    :return:  intrinsic carrier concentration of silicon (cm^-3)
    """
    return 5.29E+19 * (T / 300) ** 2.54 * np.exp(-6726 / T)


def n_ieff(N_D, N_A, Δn, T=298.15):
    """
    Return effective ni (cm-3) given donor concentration N_D=n0 (1/cm³), acceptor concentration N_A=p0 (1/cm³),
    excess carrier density (1/cm³), temperature (K), calculation of the effective intrinsic concentration n_ieff
    including BGN according to Altermatt JAP 2003.
    :param N_D: Donor concentration (cm^-3)
    :param N_A: Donor concentration (cm^-3)
    :param Δn: Excess carrier density (cm^-3)
    :param T: Temperature (K)
    :return: Effective intrinsic carrier concentration (cm^-3)
    """
    # calculation of fundamental band gap according to Pässler2002
    Eg0 = Eg0_paessler(T)

    # n_i without BGN according to Misiakos93, parameterization fits very well
    # to value of Altermatt2003 at 300K
    ni0 = ni_misiakos(T)

    ni = ni0  # ni0 as starting value for n_ieff for calculation of n0 & p0

    n0 = np.where(N_D > N_A, N_D, N_A / ni ** 2)
    p0 = np.where(N_D > N_A, N_D / ni ** 2, N_A)

    # self-conistent iterative calculation of n_ieff

    for i in range(5):  # lazy programmer as it converges pretty fast anyway
        n = n0 + Δn
        p = p0 + Δn
        dEc, dEv = bandgap_schenk(n, p, N_A, N_D, Δn, T)
        ni = ni0 * np.exp(const.q * (dEc + dEv) / (2 * const.k * T))
        n0 = np.where(N_D > N_A, N_D, N_A / ni ** 2)
        p0 = np.where(N_D > N_A, N_D / ni ** 2, N_A)

    # print('iterations',ni)
    return ni


def bandgap_schenk(n_e, n_h, N_D, N_A, Δn, T=298.15):
    """
    Returns the band gap narrowing in silicon delta conduction band, delta valence band in eV.
    Band-gap narrowing after Schenk 1998, JAP 84(3689), model described in K. McIntosh IEEE PVSC 2010, model confirmed
    by Glunz2001 and Altermatt2003, nomenclatur and formula no. according to McIntosh2010,
    table no. according to Schenk1998.
    Code adapted from Richter at Fraunhofer ISE http://dx.doi.org/10.1063%2F1.368545
    :param n_e: Total electron density with Δn (cm^-3)
    :param n_h: Total hole density with Δn (cm^-3)
    :param N_D: Donor concentration (cm^-3)
    :param N_A: Acceptor concentration (cm^-3)
    :param Δn: Excess carrier density (cm^-3)
    :param T: Temperature (K)
    :return: conduction band and valence band narrowing in silicon (ev)
    """
    # Silicon material parameters (table 1)
    m_e_ = 0.321  # m_e/m_0 -> relative electron mass
    m_h_ = 0.346  # m_h/m_0 -> relative hole mass
    g_e = 12  # degeneracy factor for electrons
    g_h = 4  # degeneracy factor for holes
    m_ = 0.1665  # µ*/m_0 -> reduced effective mass / m_0
    alfa_e = 0.5187  # µ*/m_e
    alfa_h = 0.4813  # µ*/m_h
    Ry_ex = 0.01655  # eV    excitonic Rydberg constant
    alfa_ex = 0.0000003719  # cm     excitonic Bohr radius
    epsilon_s = 11.7  # static dielectric constant

    # Parameters for Pade-Approximation (tab. 2 & 3)
    b_e = 8
    b_h = 1
    c_e = 1.3346
    c_h = 1.2365
    d_e = 0.893
    d_h = 1.153
    p_e = 7 / 30
    p_h = 7 / 30
    h_e = 3.91
    h_h = 4.2
    j_e = 2.8585
    j_h = 2.9307
    k_e = 0.012
    k_h = 0.19
    q_e = 3 / 4
    q_h = 1 / 4

    # ==========================================================================
    # pre-calculations:
    F = (const.k * T / const.q) / Ry_ex  # eq. 29
    a3 = alfa_ex ** 3

    # Normalizing of the densities
    n_e *= a3
    n_h *= a3
    N_D *= a3
    N_A *= a3

    # for eq. 33 (normalized)
    n_sum_xc = n_e + n_h
    n_p_xc = alfa_e * n_e + alfa_h * n_h

    # for eq. 37 (normalized)
    n_sum_i = N_D + N_A  # eq.39 bzw. eq. 29
    n_p_i = alfa_e * N_D + alfa_h * N_A  # eq.39 bzw. eq. 29

    Ui = n_sum_i ** 2 / F ** 2  # eq. 38
    n_ionic = n_sum_i  # McIntosh2010

    # exchange quasi-partical shift Eq33:
    delta_xc_h = -(
        (4 * const.π) ** 3 * n_sum_xc ** 2 * ((48 * n_h / (const.π * g_h)) ** (1 / 3) + c_h * np.log(1 + d_h * n_p_xc ** p_h)) + (
            8 * const.π * alfa_h / g_h) * n_h * F ** 2 + np.sqrt(8 * const.π * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * const.π) ** 3 * n_sum_xc ** 2 + F ** 3 + b_h * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)
    delta_xc_e = -(
        (4 * const.π) ** 3 * n_sum_xc ** 2 * ((48 * n_e / (const.π * g_e)) ** (1 / 3) + c_e * np.log(1 + d_e * n_p_xc ** p_e)) + (
            8 * const.π * alfa_e / g_e) * n_e * F ** 2 + np.sqrt(8 * const.π * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * const.π) ** 3 * n_sum_xc ** 2 + F ** 3 + b_e * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)

    # ionic quasi-partical shift Eq37:
    delta_i_h = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / const.π) * (1 + h_h * np.log(1 + np.sqrt(n_sum_i) / F)) + j_h * Ui * n_p_i ** 0.75 * (
            1 + k_h * n_p_i ** q_h))
    delta_i_e = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / const.π) * (1 + h_e * np.log(1 + np.sqrt(n_sum_i) / F)) + j_e * Ui * n_p_i ** 0.75 * (
            1 + k_e * n_p_i ** q_e))

    # rescale BGN
    dE_gap_h = -Ry_ex * (delta_xc_h + delta_i_h)
    dE_gap_e = -Ry_ex * (delta_xc_e + delta_i_e)
    return dE_gap_e, dE_gap_h


# ******** SubSection: Bulk Recombination *********
def U_radiative(n, p):
    """
    Returns the recombination rate caused by radiative recombination.
    :param n: electron concentration (cm^-3)
    :param p: holes concentration (cm^-3)
    :return: Radiative recombination rate (cm-3 s-1)
    """
    B_rad = 4.73e-15
    return n * p * B_rad


def U_radiative_alt(n0, p0, Δn, T=298.15):
    """
    Returns the recombination rate caused by radiative recombination.
    :param n0: Electron concentration (cm^-3)
    :param p0: Holes concentration (cm^-3)
    :param Δn: Excess of carriers (cm^-3)
    :param T: Temperature (K)
    :return: Radiative recombination rate (cm-3 s-1)
    """
    n_p = n0 + p0 + 2 * Δn
    n = n0 + Δn
    p = p0 + Δn
    B_low = 4.73e-15
    b_min = 0.2 + (0 - 0.2) / (1 + (T / 320) ** 2.5)
    b1 = 1.5E+18 + (10000000 - 1.5E+18) / (1 + (T / 550) ** 3)
    b3 = 4E+18 + (1000000000 - 4E+18) / (1 + (T / 365) ** 3.54)
    B_rel = b_min + (1 - b_min) / (1 + (0.5 * n_p / b1) ** 0.54 + (0.5 * n_p / b3) ** 1.25)
    B_rad = B_low * B_rel
    return n * p * B_rad


def U_SRH(n, p, Et, τ_n, τ_p, ni_eff=8.5e9, T=298.15):
    """
    Return the shockley read hall recombination cm-3
    given Et (eV) trap level from intrinsic
    :param n: Electron concentration (cm^-3)
    :param p: Holes concentration (cm^-3)
    :param Et: Energy of the trap from intrinsic (eV)
    :param τ_n: Electron lifetime (s)
    :param τ_p: Electron lifetime (s)
    :param ni_eff: Intrinsic carrier concentration (cm^-3)
    :param T: Temperature (K)
    :return: SRH recombination rate (cm-3 s-1)
    """
    """"""
    n1 = ni_eff * np.exp(const.q * Et / const.k / T)
    p1 = ni_eff * np.exp(-const.q * Et / const.k / T)
    return (n * p - ni_eff ** 2) / (τ_p * (n + n1) + τ_n * (p + p1))


def U_auger_richter(n0, p0, Δn, ni_eff):
    """
    Return the Auger recombination rate. 18 and 19 https://doi.org/10.1016/j.egypro.2012.07.034
    :param n0: Electron concentration (cm^-3)
    :param p0: Holes concentration (cm^-3)
    :param Δn: Excess of carriers (cm^-3)
    :param ni_eff: Intrinsic carrier concentration (cm^-3)
    :return: Auger recombination rate (cm^-3 s^-1)
    """
    B_n0 = 2.5E-31
    C_n0 = 13
    D_n0 = 3.3E+17
    exp_n0 = 0.66
    B_p0 = 8.5E-32
    C_p0 = 7.5
    D_p0 = 7E+17
    exp_p0 = 0.63
    C_dn = 3E-29
    D_dn = 0.92
    g_eeh = (1 + C_n0 * (1 - np.tanh((n0 / D_n0) ** exp_n0)))
    g_ehh = (1 + C_p0 * (1 - np.tanh((p0 / D_p0) ** exp_p0)))
    np_ni2 = (n0 + Δn) * (p0 + Δn) - ni_eff ** 2
    return np_ni2 * (B_n0 * n0 * g_eeh + B_p0 * p0 * g_ehh + C_dn * Δn ** D_dn)


def U_low_doping(n0, p0, Δn):
    """
    Recombination due to Auger and radiative. Equation 21 in DOI: 10.1103/PhysRevB.86.165202
    :param n0: Electron concentration (cm^-3)
    :param p0: Holes concentration (cm^-3)
    :param Δn: Excess of carriers (cm^-3)
    :return: Recombination rate at low doping (cm^-3 s^-1)
    """
    B_low = 4.73e-15
    n = n0 + Δn
    p = p0 + Δn
    return Δn / (n * p * (8.7e-29 * n0 ** 0.91 + 6.0e-30 * p0 ** 0.94 + 3.0e-29 * Δn ** 0.92 + B_low))


def lifetime(U, Δn):
    """
    Return the lifetime (s) where U is the recombination and Δn is the excess minority carrier density. This is
    the definition of lifetime.
    :param U: Recombination rate(cm^-3 s^-1)
    :param Δn: Excess of carriers (cm^-3)
    :return: Lifetime (s)
    """
    return Δn / U


def lifetime_auger(Δn, Ca=1.66e-30):
    """
    Returns the Auger lifetime (s) at high level injection
    given the injection level (cm-3)
    :param Δn: Excess of carriers (cm^-3)
    :param Ca: ??
    :return: Lifetime (s)
    """
    return 1 / (Ca * Δn ** 2)


def lifetime_SRH(N, Nt, Et, σ_n, σ_p, Δn, T=298.15):
    """
    Returns the lifetime when SRH is dominant over other types of recombination.
    :param N:  Doping concentration (cm^-3)
    :param Nt: Concentration of traps ?? (cm^-3)
    :param Et: Energy of the trap from the intrinsic level (eV)
    :param σ_n: Area of influence of the trap over electrons ??
    :param σ_p: Area of influence of the trap over holes ??
    :param Δn: Excess of carrier (cm^-3)
    :param T: Temperature (K)
    :return:
    """
    Nv = 31000000000000000000 * (T / 300) ** 1.85
    Nc = 28600000000000000000 * (T / 300) ** 1.58
    Eg = 1.1246
    vth = 11000000 * (T / 300) ** 0.5
    p0 = N
    n0 = (ni_Si(T) ** 2) / N
    τ_n0 = 1 / (Nt * σ_n * vth)
    τ_p0 = 1 / (Nt * σ_p * vth)
    n1 = Nc * np.exp(-Et / common.Vt())
    p1 = Nv * np.exp((-Et - Eg) / common.Vt())
    k_ratio = σ_n / σ_p
    return (τ_p0 * (n0 + n1 + Δn) + τ_n0 * (p0 + p1 + Δn)) / (n0 + p0 + Δn)


def U_surface(n, p, Sn, Sp, n1=8.3e9, p1=8.3e9, ni=8.3e9):
    """
    Return the carrier recombination (1/s) at a surface.
    :param n: Electrons concentration (cm^-3)
    :param p: Holes concentration (cm^-3)
    :param Sn: Surface recombination for electrons
    :param Sp: Surface recombination for holes
    :param n1: Electrons instrinsic concentration (cm^-3)
    :param p1: Holes instrinsic concentration (cm^-3)
    :param ni: Intrinsic carrier concentration (cm-3)
    :return:  Carrier recombination at the surface (1/s)
    """
    return Sn * Sp * (n * p - ni ** 2) / (Sn * (n + n1) + Sp * (p + p1))


# ******** Section: Solar Cells *********
def IQE_emitter(ab, We, Le, De, Se):
    """
    Return the internal quantum efficiency of a solar cell emitter.
    :param ab: Absorption coefficient (cm^-1)
    :param We: Thickness of the emitter (cm)
    :param Le: Diffusion length of minority carriers in the emitter (cm)
    :param De: Diffusivity of carriers in the emitter (cm^2 s^-1)
    :param Se: Recombination at the front surface (cm^2 s^-1)
    :return: Emitter quantum efficiency (0 to 1)
    """
    GF = ((Se * Le / De) + ab * Le - np.exp(-ab * We) * ((Se * Le / De) * np.cosh(We / Le) + np.sinh(We / Le))) / (
        (Se * Le / De) * np.sinh(We / Le) + np.cosh(We / Le)) - Le * ab * np.exp(-ab * We)
    QEE = (Le * ab / (ab * ab * Le * Le - 1)) * GF
    return QEE


def IQE_base(ab, We_Wd, Wb, Lb, Db, Sb):
    """
    Return the internal quantum efficiency of a solar cell base.
    :param ab: Absorption coefficient (cm^-1)
    :param We_Wd: Junction depth (cm)
    :param Wb: Base thickness (cm)
    :param Sb: Surface recombination velocity (cm/s)
    :param Lb: Diffusion length of minority carrier in the base (cm)
    :param Db: Diffusivity of minority carriers in the base (cm²/Vs)
    :return: Base quantum efficiency (0 to 1)
    """
    GF = (ab * Lb - (
        (Sb * Lb / Db) * (np.cosh(Wb / Lb) - np.exp(-ab * Wb)) + np.sinh(Wb / Lb) + Lb * ab * np.exp(-ab * Wb)) / (
              (Sb * Lb / Db) * np.sinh(Wb / Lb) + np.cosh(Wb / Lb)))
    QEB = (np.exp(-ab * We_Wd)) * (Lb * ab / (ab ** 2 * Lb ** 2 - 1)) * GF
    return QEB


def IQE_depletion(ab, We, Wd):
    """
    Return the internal quantum efficiency of a solar cell depletion region.
    :param ab: Absorption coefficient (cm^-1)
    :param We: Emitter thickness (cm)
    :param Wd: Depletion region thickness (cm)
    :return: Depletion region quantum efficiency (0 to 1)
    """
    QED = np.exp(-ab * We) * (1 - np.exp(-ab * Wd))
    return QED


def IQE(ab, Wd, Se, Le, De, We, Sb, Wb, Lb, Db):
    """
    Returns the quantum efficiencies by region.
    :param ab: Absorption coefficient (cm^-1)
    :param Wd: Depletion region width (cm)
    :param Se: Recombination at the front surface (cm^2 s^-1)
    :param Le: Diffusion length of minority carriers in the emitter (cm)
    :param De: Diffusivity of carriers in the emitter (cm^2 s^-1)
    :param We: Emitter thickness (cm)
    :param Sb: Recombination at the back surface (cm^2 s^-1)
    :param Wb: Base thickness (cm)
    :param Lb: Diffusion length of minority carriers in the base (cm)
    :param Db: Diffusivity of carriers in the base (cm^2 s^-1)
    :return: Quantum efficiencies (0-1)
    """
    QEE = IQE_emitter(ab, We, Le, De, Se)
    QEB = IQE_base(ab, We + Wd, Wb, Lb, Db, Sb)
    QED = IQE_depletion(ab, We, Wd)
    IQEt = QEE + QEB + QED
    return QEE, QEB, QED, IQEt


def QE2SR(wavelength, QE):
    """
    Converts a QE in units to spectral response given the wavelength (nm)
    :param wavelength: Wavelength (nm)
    :param QE: Quantum efficiency (0 to 1)
    :return: Spectral response (A/W)
    """
    return QE * wavelength / 1239.8


def SR2QE(wavelength, spectral_response):
    """
    Convert SR (A/W) to QE (unit 0 to 1) assumes that the wavelength is in  nm
    :param wavelength: Wavelength (nm)
    :param spectral_response: Spectral response (A/W)
    :return: Quantum efficiency (0 to 1)
    """
    return spectral_response * wavelength / 1239.8


def impliedV(Δn, N, T=298.15):
    """
    Return Voltage (V) where Δn is the excess carrier concentration (cm-3), N is the doping (cm-3) and
    T is the temperature (K). Implied voltage is often used to convert the carrier concentration in a lifetime
    tester to voltage.
    :param Δn: Excess carrier concentration (cm-3)
    :param N: Doping level (cm^-3)
    :param T: Temperature (K)
    :return: Voltage generated (V)
    """
    return common.Vt(T) * np.log((Δn + N) * Δn / ni_Si(T) ** 2)


def implied_carrier(V, N, ni=8.6e9, T=298.15):
    """
    Return excess carrier concentration (cm-3) given voltage and doping.
    :param V: Voltage (V)
    :param N: Doping level (cm^-3)
    :param ni: Intrinsic carrier concentration (cm^-3)
    :param T: Temperature (K)
    :return: Excess carrier concentration (cm-3)
    """
    return (-N + np.sqrt(N ** 2 + 4 * ni ** 2 * np.exp(V / common.Vt(T)))) / 2


def J0_layer(W, N, D, L, S, ni=8.6e9):
    """
    Return the saturation current density (A/cm2) for the narrow case.
    :param W: layer thickness (cm)
    :param N: doping (cm^-3)
    :param D: Diffusivity (cm^2 s^-1)
    :param L: diffusion length (cm)
    :param S: surface recombination velocity (cm s^-1)
    :param ni: intrinsic carrier concentration (cm^-3)
    :return: Saturation current density (A cm^-2)
    """
    F = (S * np.cosh(W / L) + D / L * np.sinh(W * L)) / (D / L * np.cosh(W * L) + S * np.sinh(W / L))
    return const.q * ni ** 2 * F * D / (L * N)


# def J0(ni, We, Ne, De, Le, Se, Nb, Wb, Db, Lb, Sb):
#    '''determines J0, the dark saturation current, under the narrow base diode
# condition where L > W.'''
#    Fe = (Se*np.cosh(We/Le)+De/Le*np.sinh(We*Le))/(De/Le*np.cosh(We*Le)+Se*np.sinh(We/Le))
#    Fb = (Sb*np.cosh(Wb/Lb)+Db/Lb*np.sinh(Wb*Lb))/(Db/Lb*np.cosh(Wb*Lb)+Sb*np.sinh(Wb/Lb))
#    J0 = q*ni**2*(Fe*De/(Le*Ne)+ Fb*Db/(Lb*Nb))
#    return J0


def current2gen(I):
    """
    Return generation (eh pairs/s) given current (amps)
    :param I: Current (A)
    :return: Generation (eh pairs/s)
    """
    return I / const.q


def quantum_efficiencies(ab, Wd, Se, Le, De, We, Sb, Wb, Lb, Db):
    """
    Calculates the quantum efficiency of the solar cell by regions and in total.
    :param ab: Array with absorption coefficient depending on wavelength (cm^-1)
    :param Wd: Depletion region thickness (cm)
    :param Se: Recombination rate at front surface (cm s^-1)
    :param Le: Emitter diffusion length (cm)
    :param De: Emitter diffusivity (cm^2 s-1)
    :param We: Emitter thickness (cm)
    :param Sb: Recombination rate at back surface (cm s^-1)
    :param Wb: Base thickness (cm)
    :param Lb: Base diffusion length (cm)
    :param Db: Base diffusivity (cm^2 s-1)
    :return: quantum efficiency of a solar cell (0 to 1)
    """
    QE = []
    EQE = []
    BQE = []
    DQE = []
    for element in ab:
        QEE, QEB, QED, IQEt = IQE(float(element), Wd, Se, Le, De, We, Sb, Wb, Lb, Db)
        EQE.append(QEE)
        BQE.append(QEB)
        DQE.append(QED)
        QE.append(IQEt)
    return QE, EQE, BQE, DQE


def current_in_band(wavelength, spectrum):
    """
    Returns the Current density in band.
    :param wavelength: Wavelength (nm)
    :param spectrum: Spectrum (W cm^-2 nm^-1)
    :return: Current density in band (A cm^-2 nm^-1)
    """
    previous_wavelength = wavelength[0]
    current_in_band = []
    for current_wavelength, current_power in zip(wavelength, spectrum):
        step = current_wavelength-previous_wavelength
        if current_wavelength > 1698.0 :
            print(common.nm2eV(current_wavelength))
            print(0.1 * current_power * step / common.nm2eV(current_wavelength))
        current_in_band.append(0.1*current_power*step/common.nm2eV(current_wavelength))
        previous_wavelength = current_wavelength
    return current_in_band


def photocurrent(wavelength, current_in_band, Eg):
    """
    Maximum photocurrent that can be generated assuming quantum efficiency of 1.
    :param wavelength: Wavelength (nm)
    :param current_in_band: Current density in band (A cm^-2 nm^-1)
    :param Eg: Bandgap (eV)
    :return:
    """
    accumulated_current = 0
    for current_wavelength, current_power in zip(wavelength, current_in_band):
        if common.nm2eV(current_wavelength) > Eg:
            accumulated_current += current_power
    return accumulated_current


def built_in_voltage(Nd, Na, ni, T):
    """
    Function that returns the built-in voltage of a PN junction for certain
    values.
    :param Nd: Majority carrier concentration on the N side of the junction (cm^-3)
    :param Na: Majority carrier concentration on the P side of the junction (cm^-3)
    :param ni: Intrinsic carrier concentration (cm^-3)
    :param T: Temperature (K)
    :return:  built-in voltage (V)
    """
    return common.Vt(T)*np.log(Na*Nd/(ni**2))


def depletion_region_width(Nd, Na, epsilon, V_0):
    """
    Function that returns the the depletion region width in cm
    (W = x_n + x_p) of a PN junction for certain values of:
    :param Nd: Majority carrier concentration on the N side of the junction (cm^-3)
    :param Na: Majority carrier concentration on the P side of the junction (cm^-3)
    :param epsilon: Relative permittivity of material (No dimensional)
    :param V_0: Built-in voltage (V)
    :return: depletion region width (cm)
    """
    return np.sqrt(2*const.ε*epsilon*V_0*(1/Nd+1/Na)/(100*const.q))


# processing
def phos_active(T):
    """
    Return the active limit of phosphorous in silicon given temperature (K)
    :param T: Temperature (K)
    :return: Active limit of phosphorous in silicon (??)
    """
    return 1.3e22 * np.exp(-0.37 * const.eV/(const.k * T))


def phos_solubility(T):
    """
    Return the solubility limit of phosphorous in silicon given the temperature (K).
    :param T: Temperature (K)
    :return: Solubility limit (??)
    """
    return 2.45e23 * np.exp(-0.62 * const.eV/(const.k * T))
