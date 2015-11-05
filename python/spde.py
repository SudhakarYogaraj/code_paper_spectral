import numpy
import sympy

# Boundaries of the domain.
L = -numpy.pi
R = numpy.pi

# Dimension of the kernel of the differential operator.
d = 2

# Truncation index.
N = 6

# Degree of the nonlinearity
degree = 2

x = sympy.Symbol('x')

# Differential operator associated with the equation.
def operator(function):
    return sympy.diff(function, x, 2) + function
    # return -sympy.diff(u,x,4) - sympy.diff(u,x,2)

# Normalized eigenfunctions of the differential operator.
def eigenfunctions(i):
    if i%2 == 0:
        return sympy.cos(i/2 * x) * numpy.sqrt(1/numpy.pi)
    else:
        return sympy.sin((i+1)/2 * x) * numpy.sqrt(1/numpy.pi)

# Functions used in the expansion.
def expa(i):
    if i%2 == 0:
        return sympy.cos(i/2 * x) * numpy.sqrt(1/numpy.pi)
    else:
        return sympy.sin((i+1)/2 * x) * numpy.sqrt(1/numpy.pi)

# Nonlinearity.
def nonlin(f):
    return 0.5 * sympy.diff(f*f, x, 1)

## End of user input ##

# Projection of the nonlinearity on the eigenspaces

    fs = sym(zeros(1,N));
    es = sym(zeros(1,N));
    u  = sym(0);
    for i = 1:N;
        fs(i) = eigfun(i,x);
        es(i) = expa(i,x);
        u = u + xs(i)*es(i);
    end

    nonLin = simplify(F(u));
    toInt  = simplify(expand(nonLin*transpose(fs)));
    tic; proj = simplify(int(toInt, x, L, R)/sqrt(pi)); toc
