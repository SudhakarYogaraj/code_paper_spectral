#!/usr/bin/env python
# encoding: utf-8

"""
This function builds the Problem that will be used by the C++ program.
It makes use of the python Sympy library.
"""

import os
import math
import sympy
import sympy.printing

# Create directory for temporary files
os.system("mkdir -p tmp")

# Slow variables
ns = 2

# Fast variables
nf = 2

# Creation of symbolic variables
x = [sympy.Symbol('x%d' % i) for i in range(ns)]
y = [sympy.Symbol('y%d' % i) for i in range(nf)]
g = [0.]*ns
f = [0.]*ns
h = [0.]*nf
vy = [0.]*nf
hy = [[0.]*nf]*nf
gx = [[0.]*nf]*nf
fx = [[0.]*nf]*nf

v = 0.5 * ((y[0] - 1) ** 4 + (y[1] - 1) ** 2 + 0.2*(y[0] - 1) * (y[1] - 1))

# Coefficient of the BM
s = sympy.sqrt(2)

# Solution of the cell problem
g[0] = sympy.cos(x[0]) * sympy.sin(y[0])
g[1] = sympy.cos(x[1]) * sympy.sin(y[0] + y[1])

# Non-leading order drift of fast process
h[0] = sympy.cos(x[0]) * sympy.cos(y[0]) * sympy.cos(y[1])
h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1])

# DEPENDENT VARIABLES

# Term that appears in the generator
S = s * s

# Derivative of the potential
for i in range(nf):
    vy[i] = sympy.diff(v, y[i])

# Non-normalized invariant measure
rho = sympy.exp(-v)

# Right-hand side of Poisson equation
for i in range(ns):
    for j in range(nf):
        f[i] += 0.5 * sympy.diff(S * rho * sympy.diff(g[i], y[j]), y[j]) / rho
        f[i] = sympy.simplify(f[i])

# Linear term
lin = 0
for i in range(nf):
    lin += sympy.simplify(0.25 * S * sympy.diff(v, y[i], 2))
    lin -= sympy.simplify(0.125 * S * sympy.diff(v, y[i]))

# Derivative of h
for i in range(nf):
    for j in range(nf):
        hy[i][j] = sympy.diff(h[i], y[j])

# differential of g
for i in range(ns):
    for j in range(ns):
        gx[i][j] = sympy.diff(g[i], x[j])

# differential of f
for i in range(ns):
    for j in range(ns):
        fx[i][j] = sympy.diff(f[i], x[j])

# Stardivergence
stardivh = 0
for i in range(nf):
    stardivh += sympy.simplify(vy[i] * h[i] - hy[i][i])


# VARIABLES FOR THE HMM

vyy = [[0.] * nf] * nf
fy = [[0.] * ns] * nf
drif = [0.] * (2*nf)
diff = [0.] * (2*nf)

# Extension of the y vector
y = [sympy.Symbol('y%d' % i) for i in range(2*nf)]

# y-derivative of f
for i in range(ns):
    for j in range(nf):
        fy[i][j] = sympy.diff(f[i], y[j])

# Second derivative of potential
for i in range(nf):
    for i in range(nf):
        vyy[i][j] = sympy.diff(vy[i], y[j])

# Construction of the drift
for i in range(nf):
    drif[i] = -vy[i]
    drif[nf+i] = h[i]
    for j in range(nf):
        drif[nf+i] -= vyy[i][j] * y[j]

for i in range(nf):
    diff[i] = s

# Output file
output = open('tmp/output.cpp', 'w')


def print_double(symbol, fun_name):

    # Generate strings for function declaration
    dec = "double {}(vector<double> x, vector<double> y)".format(fun_name)
    ccode = sympy.ccode(symbol, assign_to="result")

    # Write to output file
    output.write(dec + "{\n" + "    double " + ccode + "\n")
    output.write("    return result; \n}\n\n")


def allocate_function_pointer(fun_base, n, m=0):
    for i in range(n):

        # Case where output is vector
        if (m == 0):
            output.write("    {}[{}] = ".format(fun_base, i))
            output.write("{}{};\n".format(fun_base,  i))

        elif(m > 0):
            for j in range(m):
                output.write("    {}[{}][{}] = ".format(fun_base, i, j))
                output.write("{}{}{};\n".format(fun_base, i, j))

    output.write("\n")

print_double(stardivh, "Problem::stardivh")
print_double(v, "Problem::potential")
print_double(lin, "Problem::linearTerm")
print_double(rho, "Problem::rho")

for i in range(ns):
    print_double(g[i], "phi{}".format(i))
    print_double(f[i], "a{}".format(i))
    for j in range(ns):
        print_double(gx[i][j], "dxphi{}{}".format(i, j))
        print_double(fx[i][j], "dxa{}{}".format(i, j))
    for j in range(nf):
        print_double(fy[i][j], "dya{}{}".format(i, j))

for i in range(nf):
    print_double(vy[i], "dyv{}".format(i))
    print_double(h[i], "h{}".format(i))

for i in range(2*nf):
    print_double(drif[i], "drif{}".format(i))
    print_double(diff[j], "diff{}".format(i))

# Function that allocates function pointers
output.write("void Problem::init_functions() {\n\n")
allocate_function_pointer("a", ns)
allocate_function_pointer("phi", ns)
allocate_function_pointer("dxa", ns, ns)
allocate_function_pointer("dxphi", ns, ns)
allocate_function_pointer("drif", 2*nf)
allocate_function_pointer("diff", 2*nf)
allocate_function_pointer("dya", ns, nf)
output.write("}")


# for i in range(ns):
#     output.write("    fsplit[{}] = a{};\n".format(i, i))

# output.write("\n")
# for i in range(ns):
#     for j in range(ns):
#         output.write("    fxsplit[{}][{}] = dxa{}{};\n".format(i, j, i, j))

# output.write("}")

# Function declarations
# gx_init = "vector< vector<double> > Problem::phi_x\
# (vector<double> x, vector<double> y) {\n\
# "

# vy_init = "vector<double> Problem::grad\
# (vector<double> x, vector<double> y) {\n\
# vector< vector<double> > result(ns,vector<double>(ns,0.));"

# g_init = "vector<double> Problem::phi\
# (vector<double> x, vector<double> y) {\n\
# "

# f_init = "vector<double> Problem::a\
# (vector<double> x, vector<double> y) {\n\
# "

# fx_init = "vector< vector<double> > Problem::dax\
# (vector<double> x, vector<double> y) {\n\
# "

# h_init = "vector<double> Problem::fast_drift_h\
# (vector<double> x, vector<double> y) {\n\
# "

# rho_init = "double Problem::rho\
# (vector<double> x, vector<double> y) {\n\
# "

# fy_init = "vector< vector<double> > Problem::day\
# (vector<double> x, vector<double> y) {\n\
# "

# drif_init = "vector<double> Problem::drif\
# (vector<double> x, vector<double> y) {\n\
# "

# diff_init = "vector<double> Problem::diff\
# (vector<double> x, vector<double> y) {\n\
# "
