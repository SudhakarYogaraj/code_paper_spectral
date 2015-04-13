function cellprobFinal()

set(0,'DefaultAxesFontSize',15);
set(0,'defaulttextinterpreter','latex');

clc; clear all; close all;

%% Parameters of the simulation
alpha = 1;

lambda = [3;1]; q = [1;1]; 
sigma = (sqrt(q.^2./(2*lambda)));

% Symbolic variables
syms x y real;

% It is nice to declare pi as symmetric to have nice output
syms pi real positive;


%% Hermite basis of polynomials

% Maximal degree
N = 10; 

% Weight function
gamma = 1/sqrt(sym(pi))*exp(-y^2);

% Usual basis of polynomials

mon_y = sym(zeros(N + 1,1));
for i = 0:N; mon_y(i+1) = y^i; end       

% Gram matrix
G = int(mon_y*transpose(mon_y)*gamma,y,-inf,inf);

% Cholesky decomposition
U = chol(G,'upper','nocheck');

% Step to construct the basis
B = transpose(inv(U));

% Basis of Hermite polynomials
basis = simplify(B*mon_y);

% Remove 1 from basis
basis = basis(2:end);

for i = 1:length(basis)
    I = matlabFunction(gamma*basis(i));
    c(i) = integral(I,-0.5,inf);
end
c'

%% Construction of the FEM matrix

% basis = mon_y(2:end);

V = y^2 - 10*y^4*exp(-2*y^2);
rho = exp(-V);
oper = @(f) diff(V)*diff(f) - diff(f,2);

aux = exp(y-y^2); 
Z = integral(matlabFunction(rho),-inf,inf);
I = matlabFunction(aux*rho/Z);
aux = aux - integral(I,-inf,inf);
ex_sol =  aux;
f = oper(ex_sol);
ex_sol_n = matlabFunction(ex_sol);

A = zeros(N);
b = zeros(N,1);
for i = 1:N
    i/N
    for j = 1:N
        I = simplify(oper(basis(i))*basis(j)*rho);
        I = matlabFunction(I);
        A(i,j) = integral(I,-15,15);
    end
    I = matlabFunction(f*basis(i)*rho);
    b(i) = integral(I,-inf,inf);
end

s = (A\b).*basis;
s = cumsum(s);

% Solution needs to be centered
for i = 1:length(s)
    s_n = matlabFunction(s(i)*rho/Z);
    s(i) = s(i) - integral(s_n,-inf,inf);
end

err = s - ex_sol; 
err_L2 = zeros(length(err),1);
for i = 1:length(err);
    I = matlabFunction(err(i).^2*rho);
    err_L2(i) = integral(I,-inf,inf);
end
figure; loglog(1:length(err),err_L2,'r.');
xlabel('$d$'); 
ylabel('$\|\Phi-\Phi_d\|_{\rho}$');



%% Computation of the residual
res = matlabFunction(oper(s) - f);
sol = matlabFunction(s);
t = -1.5:.01:1.5;
figure; plot(t,res(t));
figure; plot(t,sol(t)); hold on; plot(t,ex_sol_n(t),'linewidth',2);
xlabel('$y$'); ylabel('Approximations of $\Phi$');

rho_n = matlabFunction(rho/Z);
figure; plot(t,rho_n(t));
xlabel('$y$'); ylabel('$e^{-V}/Z$');

% A2 = zeros(N);
% b2 = zeros(N,1);
% for i = 1:N
%     i/N
%     for j = 1:N
%         I = simplify(oper(basis(i))*basis(j)*rho);
%         I = matlabFunction(I);
%         A2(i,j) = integral(I,-inf,inf);
%     end
%     I = matlabFunction(f*basis(i)*rho);
%     b2(i) = integral(I,-inf,inf);
% % end
% 
% ba1 = mon_y(2:end);
% ba2 = basis;
% 
% s1 = (A\b)'*ba1;
% s2 = (A2\b2)'*ba2;
% t = -1:.01:1;
% tt = matlabFunction(s2-s1);
% plot(t,tt(t))


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

