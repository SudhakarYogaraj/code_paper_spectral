function cellprobFinal()

clc; clear all; close all;

%% Parameters of the simulation
alpha = 1;

lambda = [3;1]; q = [1;1]; 
sigma = (sqrt(q.^2./(2*lambda)));

% Symbolic variables
syms x w real;

% It is nice to declare pi as symmetric to have nice output
syms pi real positive;


%% Hermite basis of polynomials

% Number of variables
nvars = 2;

ys = sym(zeros(1, nvars));
for k=1:nvars;  
    ys(k) = sym(sprintf('y%d', k),'real'); 
end

% Maximal degree
N = 10; 

% Basis
basis = sym(zeros(N + 1, nvars));

% Usual basis of polynomials
mon_w = sym(zeros(N + 1,1));
for i = 0:N; mon_w(i+1) = w^i; end       

for i = 1:nvars
    if (abs(q(i)) > 1e-10)

        % Invariant density of the OU process with parameters lambda(i)
        % and q(i), associated with the i-th fast variable.
        rho = sqrt(lambda(i)/((sym(pi))*q(i)^2))*exp(-lambda(i)*w^2/q(i)^2);

        % Gram matrix
        G = int(mon_w*transpose(mon_w)*rho,w,-inf,inf);

        % Cholesky decomposition
        U = chol(G,'upper','nocheck'); 

        % Change of basis matrix
        B = transpose(inv(U));

        % First eigenfunctions of the generator of the OU process with 
        % parameters lambda(i) and q(i))
        basis(:,i) = simplify(B*mon_w);

        % Check of the orthogonality
        % simplify(int(basis(:,i)*transpose(basis(:,i))*rho,w,-inf,inf))
    else

        % In this case, the eigenfunctions of the operator are all the
        % polynomials with only one term.
        basis(:,i) = mon_w;
    end

    % Substitution of w by the variable
    basis(:,i) = subs(basis(:,i), w, ys(i));
    basis(1,i) = basis(1,i) + ys(i)/1e40;
end

% Basis in linear vector
I = 1; bvec = sym(zeros(1,factorial(N+nvars)/factorial(N)/factorial(nvars)));
for i = 1:N+1
    for j = 1:N+2-i
        bvec(I) = basis(i,1)*basis(j,2);
        I = I + 1;
    end
end

%% Construction of the FEM matrix

oper = @(f) ys(1)*lambda(1)*diff(f,ys(1)) ...
    + ys(2)*lambda(2)*diff(f,ys(2)) ...
    + 0.5*q(1)*diff(f,ys(1),2) ...
    + 0.5*q(2)*diff(f,ys(2),2);

f = sin(ys(1))*ys(2);

A = zeros(length(bvec));
for i = 1:length(bvec)
    for j = 1:length(bvec)
        I = oper(bvec(i))*bvec(i);
        
    end
end

%% functions of the problem

f = @(x,y) (x-2)*y(:,1).*y(:,2);
fd = @(x,y,z) y(:,1).*y(:,2);

%% Calculation of the Hermite coefficients of the function

index = (0:N)';

% Vector of the coefficients
c  = (zeros(length(index)));
cd = (zeros(length(index)));

%% Resolution of the cell problem

% Resolution based on the fact that the Hermite polynomials are the
% eigenfunction of the backward Kolmogorov operator.
sol = c./index/lambda;

%% Euler-Maruyama for the effective equation

% Final time
T = 1;

% time step 
dt = .01;
% Vector of times
t = 0:dt:T;

% Vectors to store the solution
x = zeros(1,length(t));
X = x;

% Initial condition
x(1) = .5;
X(1) = .5;

% Euler's method
for i = 2:length(t)    
    i
    p = 3;
    % Computation of the coefficients
    for j = 1:length(index);
        basis_j = matlabFunction(basis(j));
        I_c = @(s) basis_j(s).*f(x(i-1),s);
        randvec = [sigma(1)*randn(1000*2^(2*p),1), sigma(2)*randn(1000*2^(2*p),1)]
        c(j) = mean(I_c(sigma(1)*randn(1,1000*2^(2*p))));
        I_c = @(s) basis_j(s).*fd(x(i-1),s);
        cd(j) = mean(I_c(sigma*randn(1,1000*2^(2*p))));
    end
    
    % Brownian increment
    r = randn;
    
    phi = c./index/lambda; phi(1) = 0;
    phi_d = cd./index/lambda; phi_d(1) = 0;
    
    F = c'*phi_d;
    A0 = 2*c'*phi;
    A = sqrt(A0);
    
    % Spectral method
    x(i) = x(i-1) + F*dt + r*sqrt(dt)*A;   
    
    % Simplified equation
    X(i) = X(i-1) + X(i-1)/alpha*dt + r*sqrt(dt)*sqrt(2/alpha)*X(i-1);   
end

plot(t,x,t,X);

end

