#!/usr/bin/env python
# encoding: utf-8

"""
This function builds the Problem that will be used by the C++ program.
It makes use of the python Sympy library.
"""

import os
import sympy
import sympy.printing

# Create directory for temporary files
os.system("mkdir -p tmp")

# USER INPUT: Number of slow and fast variables
ns = 1
nf = 2
# END OF USER INPUT

# Coefficient of the BM
s = sympy.sqrt(2)
S = s * s

# Creation of symbolic Variables
x = [sympy.Symbol('x%d' % i) for i in range(ns)]
y = [sympy.Symbol('y%d' % i) for i in range(nf)]

# Solution of the Poisson equation
g = [0.]*ns
gx = [[0.]*ns for i in range(ns)]

# Right-hand side of Poisson equation
f = [0.]*ns
fx = [[0.]*ns for i in range(ns)]
fy = [[0.]*nf for i in range(ns)]

# Derivatives of potential
vy = [0.]*nf
vyy = [[0.]*nf for i in range(nf)]

# Auxiliary drift term in fast process
h = [0.]*nf
hy = [[0.]*nf for i in range(nf)]

# Extended drift and diffusion coefficients
drif = [0.] * (2*nf)
diff = [0.] * (2*nf)

# USER INPUT: solution of the cell problem
g[0] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1]) - y[0] ** 2 - y[1] ** 2 + y[0]*y[1]
# g[0] = sympy.cos(x[0]) * sympy.cos(y[0]**2 + y[1] ** 2) - y[0] ** 2 - y[1] ** 2
# g[0] = -x[0] * (y[0] ** 5 * y[1] ** 6)
# g[1] = sympy.sin(x[0]) * sympy.sin(y[0] + y[1])
# g[1] = sympy.cos(x[0] + x[1]) * sympy.sin(y[0] + y[1])
# END OF USER INPUT

for i in range(ns):
    for j in range(ns):
        gx[i][j] = sympy.diff(g[i], x[j])

# USER INPUT: Non-leading order drift of fast process
h[0] = sympy.cos(x[0]) * sympy.cos(y[0]) * sympy.cos(y[1])
h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1])
# h[0] = sympy.Symbol('0')
# h[1] = sympy.Symbol('0')
# END OF USER INPUT

stardivh = 0
for i in range(nf):
    stardivh += vy[i] * h[i] - sympy.diff(h[i], y[i])
stardivh = sympy.simplify(stardivh)

# USER INPUT: Potential
v = 0.5 * ((y[0] - 1) ** 4 + (y[1] - 1) ** 2 + 0.2*(y[0] - 1) * (y[1] - 1))
# v = 0.5 * (y[0] - 1) ** 2
# END OF USER INPUT
for i in range(nf):
    vy[i] = sympy.diff(v, y[i])
    for j in range(nf):
        vyy[i][j] = sympy.diff(vy[i], y[j])

# Invariant density
rho = sympy.exp(-v)

# Right-hand side of Poisson equation
for i in range(ns):
    for j in range(nf):
        f[i] -= 0.5 * sympy.diff(S * rho * sympy.diff(g[i], y[j]), y[j]) / rho
        f[i] = sympy.simplify(f[i])

for i in range(ns):
    for j in range(ns):
        fx[i][j] = sympy.diff(f[i], x[j])
    for j in range(nf):
        fy[i][j] = sympy.diff(f[i], y[j])

# Linear term of the Brownian motion
lin = 0
for i in range(nf):
    lin += sympy.simplify(0.25 * S * sympy.diff(v, y[i], 2))
    lin -= sympy.simplify(0.125 * S * sympy.diff(v, y[i])**2)


# Adjoint divergence of h
stardivh = 0
for i in range(nf):
    stardivh += vy[i] * h[i] - sympy.diff(h[i], y[i])
stardivh = sympy.simplify(stardivh)


# Extension of the y vector
y = [sympy.Symbol('y%d' % i) for i in range(2*nf)]

# Construction of the drift
for i in range(nf):
    drif[i] = -vy[i]
    drif[nf+i] = h[i]
    for j in range(nf):
        drif[nf+i] -= vyy[i][j] * y[j]

for i in range(nf):
    diff[i] = s

# Output file
output = open('tmp/output.gen', 'w')


# Print all the necessary functions
def print_double(symbol, fun_name):

    # Generate strings for function declaration
    dec = "double {}(vector<double> x, vector<double> y)".format(fun_name)
    ccode = sympy.ccode(symbol, assign_to="result")

    # Write to output file
    output.write(dec + "{\n" + "    double " + ccode + "\n")
    output.write("    return result; \n}\n\n")


print_double(stardivh, "stardiv_h_n")
print_double(v, "potential_n")
print_double(lin, "linearTerm_n")
print_double(rho, "zrho_n")

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


def allocate_function_pointer(fun_base, n=0, m=0):
    if(n == 0):
        # Case of scalar output
        output.write("    {} = {}_n;\n".format(fun_base, fun_base))

    elif(n > 0):
        for i in range(n):

            # Case where output is vector
            if (m == 0):
                output.write("    {}[{}] = ".format(fun_base, i))
                output.write("{}{};\n".format(fun_base,  i))

            # Matrix case
            elif(m > 0):
                for j in range(m):
                    output.write("    {}[{}][{}] = ".format(fun_base, i, j))
                    output.write("{}{}{};\n".format(fun_base, i, j))

    output.write("\n")

output.write("void Problem::init_functions() {\n\n")

# Declaration of the number of variables
output.write("    ns = {};\n".format(ns))
output.write("    nf = {};\n\n".format(nf))

# Declaration of function pointers
output.write("dyv = vector<double (*) (vector<double> x,\
             vector<double> y)> (nf);\n")
output.write("h = vector<double (*) (vector<double> x,\
             vector<double> y)> (nf);\n")
output.write("a = vector<double (*) (vector<double> x,\
             vector<double> y)> (ns);\n")
output.write("dxa = vector< vector<double (*) (vector<double> x,\
             vector<double> y)> >(ns, vector<double (*) (vector<double> x,\
             vector<double> y)> (ns));\n")
output.write("dya = vector< vector<double (*) (vector<double> x,\
             vector<double> y)> >(ns, vector<double (*) (vector<double> x,\
             vector<double> y)> (nf));\n")
output.write("phi = vector<double (*) (vector<double> x,\
             vector<double> y)> (ns);\n")
output.write("dxphi = vector< vector<double (*) (vector<double> x,\
             vector<double> y)> >(ns, vector<double (*) (vector<double> x,\
             vector<double> y)> (ns));\n")
output.write("drif = vector<double (*) (vector<double> x,\
             vector<double> y)> (2*nf);\n")
output.write("diff = vector<double (*) (vector<double> x,\
             vector<double> y)> (2*nf);\n")

# Allocation of function pointers
allocate_function_pointer("stardiv_h")
allocate_function_pointer("zrho")
allocate_function_pointer("linearTerm")
allocate_function_pointer("potential")
allocate_function_pointer("dyv", nf)
allocate_function_pointer("h", nf)
allocate_function_pointer("a", ns)
allocate_function_pointer("dxa", ns, ns)
allocate_function_pointer("dya", ns, nf)
allocate_function_pointer("phi", ns)
allocate_function_pointer("dxphi", ns, ns)
allocate_function_pointer("drif", 2*nf)
allocate_function_pointer("diff", 2*nf)
output.write("}")