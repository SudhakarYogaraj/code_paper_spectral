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

## 1D problem

# user input : dimensions
ns = 1
nf = 1
# end

# user input : diffusion
s = sympy.Symbol('1')
# end

# user input : solution
g[0] = x[0] * sympy.sin(y[0])
# end

# user input : second order drift
# h[0] = sympy.cos(x[0]) * sympy.cos(0.5*y[0]);
h[0] = sympy.Symbol('0')
# end

# user input : potential
v = y[0] ** 4 / 4. - y[0] ** 2 / 2.
# end
