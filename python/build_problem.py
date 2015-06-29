#!/usr/bin/env python
# encoding: utf-8

"""
This function builds the Problem that will be used by the C++ program.
It makes use of the python Sympy library.
"""

import os
import math
import sympy

os.system("mkdir -p tmp")

# Slow variables
ns = 2

# Fast variables
nf = 2

# Creation of symbolic variables
x = [sympy.Symbol('x%d' % i) for i in range(ns)]
y = [sympy.Symbol('y%d' % i) for i in range(nf)]
g = [None]*ns
f = [None]*ns
h = [None]*nf
vy = [None]*nf
hy = sympy.MatrixSymbol('hy', nf, nf)
gx = sympy.MatrixSymbol('gx', ns, ns)
fx = sympy.MatrixSymbol('fx', ns, ns)

v = 0.5 * ((y[0] - 1) ** 4 + (y[1] - 1) ** 2 + 0.2*(y[0] - 1) * (y[1] - 1))

# Coefficient of the BM
s = math.sqrt(2)

# Solution of the cell problem
g[0] = sympy.cos(x[0]) * sympy.sin(y[0])
g[1] = sympy.cos(x[1]) * sympy.sin(y[0] + y[1])

# Non-leading order drift of fast process
h[0] = sympy.cos(x[0]) * sympy.cos(y[0]) * sympy.cos(y[1])
h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1])

print v
print g
print h
