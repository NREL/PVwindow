#!/usr/bin/env python
# -*- coding: utf-8 -*-

import photovoltaic.semiconductors as semiconductors
import photovoltaic.common as common
import numpy as np


def I_diode(V, I0, T=298.15):
    """Return the current (A) in an ideal diode where I0 is the saturation current (A),
    V is the voltage across the junction (volts), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return I0 * np.exp(V / common.Vt(T) - 1)


def I_cell(V, IL, I0, T=298.15):
    """Return current (amps) of a solar cell
    given voltage, light generated current, I0
    also works for J0
    """
    return IL - I0 * np.exp(V / common.Vt(T))


def I_cell_Rshunt(V, IL, I0, Rshunt, T=298.15):
    """Return current (A) of a solar cell from   """
    return IL - I0 * np.exp(V / common.Vt(T)) - V / Rshunt


def V_Rseries(voltage, I, Rs):
    """Returns the voltage of a solar cells under the effect of series resistance"""
    return voltage - I * Rs


def Voc(IL, I0, n=1, T=298.15):
    """Return the open circuit voltage, Voc, (volts) from IL(A) and I0(A).
    IL and Io must be in the same units, Eg, (A), (mA) etc
    Using (mA/cm**2) uses J0 and JL instead.
    """
    return n * common.Vt(T) * np.log(IL / I0 + 1)


def V_cell(I, IL, I0, T=298.15):
    """Return the voltage (V) in an ideal solar cell where I0 is the saturation current (A),
    I is the current (A), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return common.Vt(T) * np.log((IL - I) / I0 + 1)


def short_circuit_current(current_in_band, QE):
    """
    Function that returns short circuit current of a semiconductor of a PN junction for certain
    values of:
    :param current_in_band: Current density in band (mA cm^-2 nm^-1)
    :param QE: Quantum efficiency (0 to 1)
    :return: short circuit current density (mA cm^-2)
    """

    J_SC = 0
    for c, qe in zip(current_in_band, QE):
        J_SC += c*qe
    return J_SC

def cell_params(V, I):
    """Return key parameters of a solar cell IV curve where V is a voltage array and
    I is a current array, both with type numpy.array.
    Voc (V), Isc (A), FF, Vmp(V), Imp(A) given voltage vector in (volts)
    current vector in (amps) or (A/cm²)
    If I is in (A/cm²) then Isc will be Jsc and Imp will be Jmp. No attempt is made to fit the fill factor.
    """
    Voc = np.interp(0, -I, V)
    Isc = np.interp(0, V, I)
    idx = np.argmax(V * I)
    Vmp = V[idx]
    Imp = I[idx]
    FF = Vmp * Imp / (Voc * Isc)
    Pmp = Vmp * Imp
    return Voc, Isc, FF, Vmp, Imp, Pmp


def finger_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp):
    """
    Return the fractional resistivity power loss in a finger (0 to 1)
    :param L: finger length (cm)
    :param Jmp: currrent density at the max power point in A/cm2
    :param Sf: finger spacing (cm)
    :return:
    """
    return (L ** 2 * Jmp * Sf * resistivity) / (3 * wf * df * Vmp) * 100.0


def finger_shading(wf, Sf):
    """
    Return the fractional power loss due to finger shading (0 to 1) where wf is the width of the finger and \
    Sf is the finger spacing.
    :param wf:
    :param Sf:
    :return:
    """
    return (wf / Sf) * 100.0


def finger_sheet(Sf, Jmp, Rsheet, Vmp):
    """
    Calculate the power loses due to emitter resistance
    :param Sf: Finger spacing (cm)
    :param Jmp: Max power current density (A cm^-2)
    :param Rsheet: Emitter sheet resistance (Ohm/square)
    :param Vmp: Max power voltage (V)
    :return: power loss (%)
    """
    return (Sf ** 2 * Jmp * Rsheet) / (12 * Vmp) * 100.0


def finger_total_loss(L, Jmp, Sf, resistivity, Rsheet, wf, df, Vmp):
    """Return the fractional power loss in a finger
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    Presistivity = finger_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp)
    Pshading = finger_shading(wf, Sf)
    Psheet = finger_sheet(Sf, Jmp, Rsheet, Vmp)
    return Presistivity + Pshading + Psheet, Presistivity, Pshading, Psheet


def FF(Vmp, Imp, Voc, Isc):
    """Return FFv the fill factor of a solar cell.
    given Voc - open circuit voltage (volts)"""
    return (Vmp * Imp) / (Voc * Isc)


def FF_ideal(Voc, ideality=1, T=298.15):
    """Return the FF (units)
    given Voc - open circuit voltage (volts), ideality factor, defaults to 1 (units) """
    voc = normalised_Voc(Voc, ideality, T)
    FF0 = (voc - np.log(voc + 0.72)) / (voc + 1)
    return FF0


def normalised_Voc(Voc, ideality, T=298.15):
    """Return the normalised voc of a solar cell where Voc is the open-circuit voltage, 'ideality' is the ideality factor
     and T is the temperature (K)"""
    return Voc / (ideality * common.Vt(T))


def FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
        Isc - short circuit current (amps)
        Rseries - series resistance (ohms)
        ideality factor (units)
        T - temperature (K)
    """
    # voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rs = Rseries / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - 1.1 * rs) + (rs ** 2 / 5.4)
    return FF


def FF_Rsh(Voc, Isc, Rshunt, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
    """
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - ((voc + 0.7) * FF0) / (voc * rsh))
    return FF


def FF_RsRsh(Voc, Isc, Rseries, Rshunt, ideality=1, T=298.15):
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FFRs = FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15)
    FFRsRsh = FFRs * (1 - ((voc + 0.7) * FFRs) / (voc * rsh))
    return FFRsRsh


def efficiency(Voc, Isc, FF, A=1):
    """Return the efficiency of a solar cell (units not percentage)given Voc (volts), Isc in (amps) and  FF (units).
    also works for Jsc since area of 1 is assumed
    """
    return 1000 * Voc * Isc * FF / A
# Losses

def general_resistivity(Nx,Dx,T):
    """
    Calculates resistivity of a semiconductor based on:
    :param Nx: Doping (cm^-3)
    :param Dx: Diffusivity (cm^2 s^-1)
    :param T: Temperature (K)
    :return: resistivity (Ohm*cm)
    """
    return 1/(semiconductors.q*Nx*Dx/semiconductors.Vt(T))


def base_losses(base_resistivity, base_thickness, Vmp, Jmp):
    """
    Calculate the efficiency drop due to base losses:
    :param base_resistivity: Resistivity of the base (Ohm/cm)
    :param base_thickness: Thickness of the base (cm)
    :param Vmp: Max power voltage (V)
    :param Jmp: Max power Current density (mA/cm^2)
    :return: Power loss due to base losses (%)
    """
    return base_resistivity*base_thickness*Jmp*0.001/Vmp*100


def base_resistance(base_resistivity, base_thickness):
    """
    Calculate the base resistance from the following parameters:
    :param base_resistivity: Resistivity of the base (Ohm*cm)
    :param base_thickness: Thickness of the base (cm)
    :return: resistance of the base (Ohm*cm^2)
    """
    return base_resistivity*base_thickness


def sheet_resistance(emitter_resistivity, thickness_e):
    """
    Calculate sheet resistance of the emitter
    :param emitter_resistivity: Resistivity of the emitter (Ohm*cm)
    :param thickness_e: Thickness of the emitter (um)
    :return: Emitter sheet resistance (Ohm/square)
    """
    return emitter_resistivity/thickness_e


def emitter_resistance(R_sheet, Sf):
    """
    Calculate emitter resistance
    :param R_sheet: Emitter sheet resistance (Ohm/square)
    :param Sf: Finger spacing (cm)
    :return: emitter resistance (Ohm*cm^2)
    """
    return R_sheet*(Sf**2)/12


def finger_resistance(resistivity, Sf, wf, Lf, df):
    """
    Calculates finger resistance
    :param resistivity: Resistivity of the metal (Ohm*cm)
    :param Sf: finger_spacing (cm)
    :param wf: finger_width (cm)
    :param Lf: finger_length (cm)
    :param df: finger_height (cm)
    :return: finger_resistance (Ohm cm^2)
    """
    return resistivity*Lf*Lf*Sf/(3*wf*df)


def bus_losses(W, Z, resistivity, Jmp, Vmp, m, wb, db):
    """
    Calculates the fractional power loss due to bus resistance
    :param W: Solar cell width (cm)
    :param Z: Solar cell height (cm)
    :param resistivity: metal resisitivity (Ohm*cm)
    :param Jmp: Max power Current density (A/cm^2)
    :param Vmp: Max power voltage (V)
    :param m: Number of buses (cm)
    :param wb: width of the bus (cm)
    :param db: height of the bus (cm)
    :return: fractional power loss due to bus resistance (%)
    """
    return W*W*Z*resistivity*Jmp/(3*m*wb*db*Vmp)*100


def bus_resistance(W, resistivity, m, wb, db):
    """
    Calculates the fractional power loss due to bus resistance
    :param W: Solar cell width (cm)
    :param resistivity: metal resisitivity (Ohm*cm)
    :param m: Number of buses (cm)
    :param wb: width of the bus (cm)
    :param db: height of the bus (cm)
    :return: fractional power loss due to bus resistance (%)
    """
    return W*W*W*resistivity/(3*m*wb*db)


def transparency(W, Z, m, n, wf, wb):
    """
    Calculates the percentage of the cell area area covered by metal
    :param W: Solar cell width (cm)
    :param Z: Solar cell height (cm)
    :param m: Number of buses
    :param n: Number of fingers
    :param wf: Width of the finger (cm)
    :param wb: Width of the bus (cm)
    :return: Area of the cell covered by metal (%)
    """
    return (1-(m*wb*W+n*wf*Z-wf*wb*m*n)/(W*Z))*100




