import sympy

# USER INPUT: Number of slow and fast variables
ns = 1
nf = 1
# END OF USER INPUT


# USER INPUT: solution of the cell problem
# g[0] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1]) - y[0] ** 2 - y[1] ** 2 + y[0]*y[1]
# g[0] = sympy.cos(x[0]) * sympy.cos(y[0]**2 + y[1] ** 2) - y[0] ** 2 - y[1] ** 2
g[0] = sympy.cos(x[0] + y[0])
# g[0] = -x[0] * (y[0] ** 5 * y[1] ** 6)
# g[1] = sympy.sin(x[0]) * sympy.sin(y[0] + y[1])
# g[1] = sympy.cos(x[0] + x[1]) * sympy.sin(y[0] + y[1])
# END OF USER INPUT


# USER INPUT: Non-leading order drift of fast process
h[0] = sympy.cos(x[0]) * sympy.cos(y[0]) * sympy.exp(y[0])
# h[1] = sympy.cos(x[0]) * sympy.cos(y[0] + y[1])
# h[0] = sympy.Symbol('0')
# h[1] = sympy.Symbol('0')
# END OF USER INPUT


# USER INPUT: Potential
# v = 0.5 * ((y[0] - 1) ** 4 + (y[1] - 1) ** 2 + 0.2*(y[0] - 1) * (y[1] - 1))
v = 0.5 * (y[0] - 1) ** 2
# END OF USER INPUT
