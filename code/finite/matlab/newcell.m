function newcell()

format longEng

% Symbolic variables
syms y x pi;

% Right hand side
f = cos(1.2) * (2*y*cos(y) + sin(y));

% Potential
V = y*y + log(pi)/2.

% Associated invariant measure
rho = exp(-V)

% Standard gaussian used for basis
sigma = 0.5; %1/sqrt(2.);
gaussian = 1/sqrt(2*sym(pi)*sigma^2) * exp(-y^2/(2*sigma^2))

% --
s = sqrt(2);
S = 2;

% Maximal degree
N = 8;

% Usual basis of polynomials
mon_y = sym(zeros(N + 1,1));
for i = 0:N; mon_y(i+1) = y^i; end       

% Gram matrix
G = int(mon_y*transpose(mon_y)*gaussian,y,-inf,inf);

% Cholesky decomposition
U = chol(G,'upper','nocheck'); 

% Change of basis matrix
B = transpose(inv(U));

% First eigenfunctions of the generator of the OU process with 
% parameters lambda(i) and q(i))
basis = simplify(B*mon_y)

% Check of the orthogonality
% simplify(int(basis(:,i)*transpose(basis(:,i))*rho,y,-inf,inf));

% Linear term
lin = 1/4*S*diff(V,y,2) - 1/8*S*(diff(V,y)^2)
linn = matlabFunction(lin);

% Generator
L = @(f) -0.5*S*diff(f,y,2) - lin*f;

coefficients = zeros(N+1, 1);
for i = 0:N
    toint = simplify(f*sqrt(rho)*basis(i+1)*sqrt(gaussian));
    coefficients(i + 1) = integral(matlabFunction(toint), -inf, inf)
end

mat = zeros(N + 1, N + 1);
mat2 = zeros(N + 1, N + 1);
for i = 0:N
    for j = 0:N
        toint = (1/(2*sigma^2) - y^2/(4*sigma^4) - lin) * basis(i+1) * basis(j+1) * gaussian;
        matfunc = matlabFunction(simplify(toint));
        mat(i+1, j+1) = integral( matfunc, -inf, inf)
        
        %mat(i+1 , j+1) = double( int ( toint, -inf, inf) )
        toint = simplify( L(basis(i+1)*sqrt(gaussian)) * basis(j+1) * sqrt(gaussian) );
        matfunc = matlabFunction(toint);
        mat2(i+1 , j+1) = integral ( matfunc, -inf, inf)
    end
    mat(i+1,i+1) = mat(i+1,i+1) + i/(sigma^2);

end

w = zeros(1, N+1);
for i = 0:N
    toint = matlabFunction( simplify(basis(i + 1) * sqrt(rho) * sqrt(gaussian)) );
    w(i+1) = integral( toint, -inf, inf);
end

A = 0*mat;
b = 0*coefficients;
A(1,1) = 1.;
for i = 2:N+1
    for j = 2:N+1
        A(i,j) = mat(i,j) - w(i)/w(1) * mat(i,1) - w(j)/w(1) * mat(j,1)    + w(i)/w(1) * w(j)/w(1) * mat(1,1);
        A(1,i) = 0.;
        A(i,1) = 0.;
        b(i) = coefficients(i) - w(i)/w(1) * coefficients(1);
    end
end

sol = A\b;
result = sol;
for i = 2:N+1
    result(1) = result(1) - w(i)/w(1) * sol(i)
end

end


