clc; clear all; close all;

set(0,'DefaultAxesFontSize',14);

%% Parameters of the simulation
lambda = 2; q = 2; sigma = sqrt(q^2/(2*lambda));

% Symbolic variables
syms x y real;

% It is nice to declare pi as symmetric to have nice output
syms pi real positive;


%% Hermite basis of polynomials

% Maximal degree
N = 20; 

% Weight function
rho = sqrt(lambda/((sym(pi))*q^2))*exp(-lambda*y^2/q^2);

% Usual basis of polynomials
tic
mon_y = sym(zeros(N + 1,1));
for i = 0:N; mon_y(i+1) = y^i; end       

% Gram matrix
G2 = int(mon_y*transpose(mon_y)*rho,y,-inf,inf);
toc

tic
% Vector used to construct polynomials
aux_vec = sym(zeros(2*N+1,1));
for i = 1:size(aux_vec)
    aux_vec(i) = int(y^(i-1)*rho,y,-inf,inf);
end

% Gram matrix
G = sym(zeros(N+1,N+1));
for i = 1:N+1
    for j = 1:N+1
        G(i,j) = aux_vec(i+j-1);
    end
end
toc

% Cholesky decomposition
U = chol(G,'upper','nocheck');

% Step to construct the basis
B = transpose(inv(U));

% Basis of Hermite polynomials
basis = simplify(B*mon_y);

%% functions of the problem

% Small parameter
eps = .05;

% 1/eps drift term of the slow process
fy = 1/sqrt(rho)*exp(-abs(y));
fx = cos(x);
f  = fx*fy;

% 1/eps drift term of the fast process
h  = 10*exp(-y^2)*sin(x)*x^2;

% 1/eps^2 term of the fast process
g = @(x,y) -lambda*y;

% 1/eps diffusion term of the fast process
dy = @(x,y) q/eps;

% Non symbolic drift terms
aux1 = matlabFunction(fx);
aux2 = matlabFunction(fy);
aux3 = matlabFunction(h);
fn = @(x,y) 1/eps*aux1(x)*aux2(y);
gn = @(x,y) 1/eps^2*g(x,y) + 1/eps*aux3(x,y);


%% Calculation of the Hermite coefficients of the function

index = (0:N)';

% Vector of the coefficients
c  = (zeros(length(index),1));
cd = (zeros(length(index),1));
c1 = c;
c2 = c;

% Differentiation of the function
n = 3; fd = diff(fy,y,n);

% Numerical (non-symbolic) functions
f_n   = matlabFunction(fy);
f_nrho   = matlabFunction(simplify(fy*rho));
rho_n = matlabFunction(rho);
fd_n  = matlabFunction(fd);

method = 'NUMERICAL';
switch method
    
    % Symbolic integration for the projection
    case 'ANALYTIC' 
        for i = 1:length(index);
            c(i) = int(fy*rho*basis(i),y,-inf,inf) 
        end
        
        for i = 1:length(index);
            cd(i) = int(fd*rho*basis(i),-inf,inf)
        end
        
    % Numerical integration (faster but less precise)
    case 'NUMERICAL'
        for i = 2:length(index);
            basis_i = matlabFunction(basis(i));
            I_c = @(s) basis_i(s).*f_nrho(s);
            c(i) = integral(I_c,-inf,inf)             
        end
        
        for i = 1:length(index);
            basis_i = matlabFunction(basis(i));
            I_c = matlabFunction(basis(i)*rho*fd);
            cd(i) = integral(I_c,-inf,inf)
        end
    case 'MONTECARLO'
        for i = 2:length(index);
            basis_i = matlabFunction(basis(i));
            I_c = @(s) basis_i(s).*f_n(s);
            tic
            p = 2; c1(i) = mean(I_c(sigma*randn(1,10000)))
            toc
            tic
            p = 4; c2(i) = mean(I_c(sigma*randn(1,20000)))
            toc
            tic
            I_c = @(s) basis_i(s).*rho_n(s).*f_n(s);
            cn(i) = integral(I_c,-inf,inf); cn'       
            toc
        end
        mean((cn'-c1).^2)
        mean((cn'-c2).^2)
end

plots = 1;
if plots
    % Approximate function
    fN = simplify(cumsum(c.*basis));

    % Plot of the difference
    s = -5:.01:5;

    % Construction of the whole function
    gN = matlabFunction(fN(2:end));

    % Visualization
    figure; plot(s,f_n(s),s,gN(s)); 
end

% Sobolev norm of the derivative
au_c = matlabFunction(fd*fd*rho);
norm = sqrt(integral(au_c,-inf,inf));

%% Checking proposition 2.2

% ratio given by the proposition
R = inf*ones(length(index),1);
for i = n+1:length(index)
    R(i) = sigma^n * sqrt(factorial(index(i)-n)/factorial(index(i)));
end

% Check that the formula holds
alter = [zeros(n,1); cd(1:end-n).*R(n+1:end)];
max(abs(c(n+1:end) - alter(n+1:end)))


%% Resolution of the cell problem

% Resolution based on the fact that the Hermite polynomials are the
% eigenfunction of the backward Kolmogorov operator.
sol = c./index/lambda;

% Constraints that the solution of the cell problem integrates to zero with
% respect to the invariant measure of the fast process. 
sol(1) = 0;

% Backward Kolmogorov operator
L = @(func) -lambda*y*diff(func,y) + 0.5*q^2*diff(func,y,2);

% Approximate solution
Phi = simplify(sol.*basis);
Phi_tot = simplify(cumsum(Phi));
Phi_n = matlabFunction(Phi_tot(2:end));
figure; plot(s,Phi_n(s));

% Adding the dependence in x
syms x;
Phi = Phi*fx;
fN  = fN*fx;

% Drift term
aux1 = f*diff(Phi,x); aux1 = simplify(int(aux1*rho,y,-inf,inf));
aux2 = h*diff(Phi,y); aux2 = simplify(int(aux2*rho,y,-inf,inf));
dr = aux1 + aux2;

% Diffusion term
di = sqrt(2*simplify(int(f*Phi*rho,y,-inf,inf)));

% Sum to have the contributions of all modes
dr = cumsum(dr); di = cumsum(di);

% Vizualization of the drift coefficient
x = -6:.1:6;
dr_n = matlabFunction(dr(2:end));
figure; plot(x,dr_n(x));

%% Euler-Maruyama for the effective equation

% Final time
T = 1;

% Cell that will contain all the drift and diffusion functions
Dr = cell(length(dr)-1,1);
Di = cell(length(dr)-1,1);

for i = 2:length(dr)
    Dr{i-1} = matlabFunction(dr(i));   
    Di{i-1} = matlabFunction(di(i));
end

% time step 
dt = eps*eps*.01;

% Vector of times
t = 0:dt:T;

% Vectors to store the solution
x = zeros(length(dr)-1,length(t));
xe = zeros(1,length(t));
ye = zeros(1,length(t));

% Initial condition
x(:,1) = .5;
xe(1)  = .5;
ye(1) = 1;

% Euler's method
for i = 2:length(t)    
    
    % Brownian increment
    r = randn;
    
    % Spectral method
    for j = 1:length(x(:,1)) 
        x(j,i) = x(j,i-1) + Dr{j}(x(j,i-1))*dt ...
                 + r*sqrt(dt)*Di{j}(x(j,i-1));        
    end
    
    % Exact solution
    xx = xe(i-1);
    yy = ye(i-1);
    
    xe(i) = xx + fn(xx,yy)*dt;
    ye(i) = yy + gn(xx,yy)*dt + dy(xx,yy)*sqrt(dt)*r;
end

figure; plot(t,x); hold on; plot(t,xe,'linewidth',2);
