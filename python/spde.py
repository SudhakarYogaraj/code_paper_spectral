#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy
import numpy

# Boundaries of the domain.
L = -sympy.pi
R = sympy.pi

# Dimension of the kernel of the differential operator.
ns = 2

# Number of fast modes to keep
nf = 4

# Degree of the nonlinearity
degree = 4

# Space variable
x = sympy.Symbol('x')

# Differential operator associated with the equation.
def operator(function):
    return sympy.diff(function, x, 2) + function
    # return -sympy.diff(u,x,4) - sympy.diff(u,x,2)


# Normalized eigenfunctions of the differential operator. (i = 1, 2 ...)
def eigenfunctions(i):
    if (i + 1) % 2 == 0:
        return sympy.cos(sympy.Integer((i+1)/2) * x) * sympy.sqrt(1/sympy.pi)
    else:
        return sympy.sin(sympy.Integer((i+2)/2) * x) * sympy.sqrt(1/sympy.pi)


# Functions used in the expansion. (i = 1, 2 ...)
def expa(i):
    if (i + 1) % 2 == 0:
        return sympy.cos(sympy.Integer((i+1)/2) * x)
    else:
        return sympy.sin(sympy.Integer((i+2)/2) * x)


# Nonlinearity.
def nonlinearity(f):
    # return sympy.Rational(0.5) * sympy.diff(f*f, x, 1)
    return sympy.Rational(0.5) * f*f * sympy.diff(f*f, x, 1)

# END OF USER INPUT

# Coefficients of the expansion
coefficients = [sympy.Symbol('c%d' % i) for i in range(ns + nf)]

# Space variable
x = sympy.Symbol('x')
u = 0

exp = []
eig = []

for i in range(ns + nf):
    exp.append(expa(i))
    eig.append(eigenfunctions(i))
    u += coefficients[i] * exp[i]

# Calculation of the nonlinearity
nonlinear_term = sympy.collect(sympy.expand(sympy.expand_trig(nonlinearity(u))), [sympy.sin(x), sympy.cos(x)])

print("Nonlinearity:")
print(nonlinear_term)

# Projection of the nonlinearity
# TODO: justify division by sqrt(pi)
projections = []
for i in range(ns + nf):
    product = sympy.expand(nonlinear_term * sympy.expand_trig(eig[i]))
    result_integration = sympy.integrate(product, (x,L,R)) / sympy.sqrt(sympy.pi)
    projections.append(sympy.simplify(result_integration))
    print(projections[i])
