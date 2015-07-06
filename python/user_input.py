import sympy

# user input : dimensions
ns = 1
nf = 1
# end


# user input : solution
g[0] = sympy.cos(x[0] + y[0])
# end


# user input : second order drift
h[0] = sympy.cos(x[0]) * sympy.cos(y[0]) * sympy.exp(y[0])
# end


# user input : potential
v = 0.5 * (y[0] - 1) ** 2
# end
