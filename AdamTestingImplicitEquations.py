# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:07:58 2021

@author: aduell
"""

import numpy as np
import tmm
import pandas as pd
#import tmm_vw as tmm
import matplotlib.pyplot as plt
from wpv import Layer, Stack
import scipy.interpolate, scipy.integrate, pandas, sys
from numericalunits import W, K, nm, m, cm, s, eV, meV, V, mA, c0, hPlanck, kB, e, ohm, A, J, C
assert sys.version_info >= (3,6), 'Requires Python 3.6+'
import scipy
from scipy.optimize import fsolve

import sympy
import sympy.solvers.solvers


#x = sympy.Symbol('x')
#w = sympy.solvers.solve(x**2 - 1, x)
#print (w)


Rs = .01 * ohm
Rsh = 1e4 * ohm
Tcell = 300 * K
voltage = 6 * V


#u = sympy.solvers.solve(sympy.exp(I)*2*I - I, I)
#print(u)

#(((e * np.exp((voltage + I * Rs) / (kB * Tcell)) - (voltage + I * Rs)/Rsh)) -I, I)
#print (w)
#I = sympy.Symbol('I')
#test = sympy. solvers.solve(sympy.exp(I)*I-I,I)
#Current = sympy.solvers.solve((sympy.exp(I) - (I) - I, I))
#Current = sympy.solvers.solve((sympy.exp(I) - ( I * Rs)/Rsh) - I, I)
#Current = sympy.solvers.solve((e * sympy.exp(e * (voltage + I * Rs) / (kB * Tcell)) - (voltage + I * Rs)/Rsh) - I, I)
#print(Current)

#def current_density(voltage, eta,Absorbed):
#    I =sympy.Symbol('I')
#    Current = sympy.solvers.solve(e * (Generated - RR0/eta_ext * np.exp(e * (voltage + I * Rs) / (kB * Tcell)) - (voltage + I * Rs)/Rsh) -I, I)
#    return Current

#current_density(6, .9, 50)

#sympy.exp(0)-0-0
#print(test)

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define the expression whose roots we want to find


#unc = lambda I : (e * (np.exp(e * (voltage + I*A * Rs) / (kB * Tcell))) - I*A)
func = lambda I : (e * (np.exp(e * (voltage + I * Rs) / (kB * Tcell)) - (voltage + I * Rs)/Rsh)) - I

# Plot it

print('e = ',e/ C)
print('A = ', A)
initial_guess = -460 * A
solution = fsolve(func, initial_guess)

#rootsolution = scipy.optimize.root(func, initial_guess)

print ("The solution is I = ",  (solution / A))
print ("at which the value of the expression is ", func(solution / A))

#rootsolution = scipy.optimize.root(func, initial_guess)

#print ('using root solution is', rootsolution)


#def numerical_solve_variable(Function, first_guess = 0):
#    initial_guess = first_guess
#    proposed_solution = scipy.optimize.fsolve(Function, initial_guess)
#    solution2 = scipy.optimize.fsolve(Function, proposed_solution)
#    if solution2 == proposed_solution:
#        return solution2
#    else :
#        return numerical_solve_variable(Function, first_guess = solution2)
#WERT = numerical_solve_variable(func)
#print(WERT)


