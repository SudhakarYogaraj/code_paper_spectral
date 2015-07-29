# user input : dimensions
ns = 1
nf = 1
# end

# user input : solution
g[0] = sympy.cos(x[0] + y[0])
# g[0] = sympy.cos(x[0]) * sympy.sin(y[0]);
# g[0] = sympy.cos(x[0]) * sympy.sin(y[0] + y[1]);
# end

# user input : second order drift
h[0] = sympy.cos(x[0]) * sympy.cos(y[0])
# h[0] = sympy.cls(x[0]) * sympy.cos(y[0]) * sympy.cos(y[1]);
# h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1]);
# h[0] = sympy.Symbol('0');
# h[1] = sympy.Symbol('0');
# end

# user input : potential
v =  (y[0] - 10) ** 4 / 4. - (y[0] - 10) ** 2 / 2.
# v = 0.5 * ((y[0] - 1) ** 4 + (y[1]-1) ** 2 + 0.2*(y[0] - 1)*(y[1] -1))
# end
