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

## 2D problem

# user input : dimensions
ns = 2
nf = 2
# end

# user input : diffusion
s = sympy.sqrt(2)
# end

# user input : solution
g[0] = x[0] * sympy.cos(y[0])
g[1] = x[1] * sympy.cos(y[1])
# end

# user input : second order drift
h[0] = sympy.Symbol('0')
h[1] = sympy.Symbol('0')
# end

# user input : potential
v = y[0] ** 2 + y[1] ** 2
# end
