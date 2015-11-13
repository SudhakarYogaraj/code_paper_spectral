# Description of the variables
#
# ns : number of slow processes
# nf : number of fast processes
#
# g : solution of the cell problem (dimension = ns)
#
# h : Non-leading order part of the drift of the fast processes.
# (The one that multiplies 1/eps.)
#
# v : The potential of the problem
#
# s : Diffusion coefficient
#

## 3D problem

# user input : dimensions
ns = 2
nf = 3
# end

# user input : diffusion
s = sympy.sqrt(2)
# end

# user input : solution
g[0] = sympy.cos(x[0] + y[0] + y[1])
g[1] = sympy.sin(x[1]) * sympy.sin(y[0] + y[1] + 0.5*y[2]);
# end

# user input : second order drift
h[0] = sympy.cos(x[0]) * sympy.cos(y[0] + y[2])  * sympy.cos(y[1]);
h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1]);
h[2] = sympy.Symbol('0');
# end

# user input : potential
v = (y[0] ** 4/4. - y[0] ** 2/2. + y[1] ** 4 + y[2] ** 2)
# end
