% close all;
% index = 0:20;
% r = sigma.^index.*sqrt(1./factorial(index));
% c = exp(0.5)*r;
% for i = 2:19
% n = i; 
% bound = sigma^n*sqrt(n^n/factorial(n))*index.^(-n/2)*norm;
% bound2 = sigma.^index.*sqrt(1./factorial(index))*norm;
% loglog(index(n+1:end), abs(c(n+1:end)), 'r.' , index(n+1:end), abs(bound(n+1:end)),'b-');
% xlabel('Index of the coefficients', 'FontSize',15,'inTerpreter','Latex');
% ylabel('Value of the coefficiens','Fontsize', 15,'interpreter','latex');
% hold on; end

clc; clear all; close all;

%% Parameters of the simulation
lambda = 2; q = 2; sigma = sqrt(q^2/(2*lambda));

% Symbolic variable
syms x real;

% Final time
T = 1;

% Small parameter
eps = .05;

% Drifts
fx = @(x,y) cos(x)*sin(y)/eps;
fy = @(x,y) -lambda*y/(eps*eps) + exp(x)*sin(y)/eps;

% Diffusion fast process
gy = @(x,y) q/eps;


%% Simulation using forward Euler method

% time step 
dt = eps*eps*.01;

% Vector of times
t = 0:dt:T;

% Vector to store the solution
x = zeros(1,length(t));
y = zeros(1,length(t));

% Initial condition
x(1) = .5;

% Euler's method
for i = 2:length(t)
    
    xx = x(i-1);
    yy = y(i-1);
    
    x(i) = xx + fx(xx,yy)*dt;
    y(i) = yy + fy(xx,yy)*dt + gy(xx,yy)*sqrt(dt)*randn;
end

figure; plot(t,x);

%% Hermite basis of polynomials

% Maximal degree
N = 20; 

% Symbolic variables for the computation
syms w; syms pi real positive;

% Weight function
rho = sqrt(lambda/((sym(pi))*q^2))*exp(-lambda*w^2/q^2);

if (exist('cellprob_basis.mat','file'))
    load('cellprob_basis')
else
    % Usual basis of polynomials
    mon_w = sym(zeros(N + 1,1));
    for i = 0:N; mon_w(i+1) = w^i; end       

    % Gram matrix
    G = int(mon_w*transpose(mon_w)*rho,w,-inf,inf);

    % Cholesky decomposition
    U = chol(G,'upper','nocheck');

    % Step to construct the basis
    B = transpose(inv(U));

    % Basis of Hermite polynomials
    basis = simplify(B*mon_w);
    
    % Save basis in file
    save('cellprob_basis','basis');
end

%% Construction of the function

% Regularity C^deg
n = 3; 

% Position of the discontinuity in d^{deg+1}f/dt^{deg+1}
r = 100;  % Left 
l = -100; % Right

% Function between l and r
f = sin(w);

% Value of the function for x < l and x > r
fr = sym(0); fl = sym(0);

for i = 0:n
    fr = fr + subs(diff(f,w,i),w,r)*(w-r)^i/factorial(i);
    fl = fl + subs(diff(f,w,i),w,l)*(w-l)^i/factorial(i);
end

%% Calculation of the Hermite coefficients of the function

index = (0:N)';
c = (zeros(length(index),1));
cd = (zeros(length(index),1));

% Differentiation of the function
fd = diff(f,w,n);
fdl = diff(fl,w,n);
fdr = diff(fr,w,n);

% Numerical function
fl_n = matlabFunction(fl);
f_n  = matlabFunction(f);
fr_n = matlabFunction(fr);
rho_n = matlabFunction(rho);
fd_n  = matlabFunction(fd);
fdl_n = matlabFunction(fdl);
fdr_n = matlabFunction(fdr);

method = 'NUMERICAL';
switch method
    
    % Symbolic integration for the projection
    case 'ANALYTIC' 
        for i = 1:length(index);
            c(i) = int(simplify(fl*rho*basis(i)),w,-inf,l) ...
                   + int(f*rho*basis(i),w,l,r) ...
                   + int(simplify(fr*rho*basis(i)),w,r,inf)
        end
        
        for i = 1:length(index);
            cd(i) = int(simplify(fdl*rho*basis(i)),-inf,l) ...
                    + int(fd*rho*basis(i),l,r) ...
                    + int(simplify(fdr*rho*basis(i)),r,inf)
        end
        
    % Numerical integration (faster but less precise)
    case 'NUMERICAL'
        for i = 2:length(index);
            basis_i = matlabFunction(basis(i));
            I_c = @(s) basis_i(s).*rho_n(s).*f_n(s);
            I_l = @(s) basis_i(s).*rho_n(s).*fl_n(s);
            I_r = @(s) basis_i(s).*rho_n(s).*fr_n(s);
            c(i) = integral(I_l,-inf,l) ...
                   + integral(I_c,l,r) ...
                   + integral(I_r,r,inf)
        end
        
        for i = 1:length(index);
            basis_i = matlabFunction(basis(i));
            I_c = matlabFunction(basis(i)*rho*fd);
            I_l = matlabFunction(basis(i)*rho*fdl);
            I_r = matlabFunction(basis(i)*rho*fdr);
            cd(i) = integral(I_l,-inf,l) ...
                    + integral(I_c,l,r) ...
                    + integral(I_r,r,inf)
        end
end

plots = 0;
if plots
    % Approximate function
    fN = simplify(cumsum(c.*basis));

    % Plot of the difference
    t = -5:.01:5;

    % Construction of the whole function
    g  = @(s) (s < l).*fl_n(s) + (s >= l).*(s <= r).*f_n(s) + (s > r).*fr_n(s);
    gN = matlabFunction(fN(2:end));

    % Visualization
    figure; plot(t,g(t),t,gN(t)); 
end

% Sobolev norm of the derivative
au_l = matlabFunction(fdl*fdl*rho);
au_c = matlabFunction(fd*fd*rho);
au_r = matlabFunction(fdr*fdr*rho);
norm = sqrt( integral(au_l,-inf,l) ...
             + integral(au_c,l,r) ...
             + integral(au_r,r,inf))

%% Checking proposition 2.3

% ratio given by the proposition
R = inf*ones(length(index),1);
for i = n+1:length(index)
    R(i) = sigma^n * sqrt(factorial(index(i)-n)/factorial(index(i)));
end
% Check that the formula holds
alter = [zeros(n,1); cd(1:end-n).*R(n+1:end)];
max(abs(c(n+1:end) - alter(n+1:end)));

C = sqrt(factorial(n))*sigma^n;
bound = C*index.^(-n/2)*norm;
figure; plot(log(index(n+1:end)),log(abs(c(n+1:end))),'r.',...
    log(index(n+1:end)),log(abs(bound(n+1:end))),'b-');

%% Resolution of the cell problem

% Resolution based on the fact that the Hermite polynomials are the
% eigenfunction of the backward Kolmogorov operator.
sol = c./index/lambda;

% Note that c(1) = 0, and the first eignevalue is hence sol(0) = 0
sol(1) = 0;

% Backward Kolmogorov operator
L = @(func) -lambda*w*diff(func,w) + 0.5*q^2*diff(func,w,2);

% Approximate solution
Phi = simplify(cumsum(sol.*basis));
Phi_n = matlabFunction(Phi(2:end));
figure; plot(t,Phi_n(t));

% Adding the dependence in xi
syms xi;
Phi = Phi*cos(xi);
f   = f*cos(xi);
fN  = fN*cos(xi);

% 1/eps drift term of the fast process
h = exp(xi)*sin(w);

% Drift term
au = f*diff(Phi,xi) + h*diff(Phi,w);
dr = simplify(int(au*rho,w,-inf,inf));

% Vizualization
x = -10:.1:10;
dr_n = matlabFunction(dr(2:end));
figure; plot(x,dr_n(x));

%loglog(1:length(c)-1,abs(c(2:end))); hold on; loglog(1:length(c)-1, bo(2:end));
%loglog(1:length(c)-1,pb(2:end));

% n = 6;
% g = @(x) (x < l).*fl(x) + (x >= l).*(x <= r).*f(x) + (x > r).*fr(x);
% 
% c = zeros(length(basis),1);
% cn = zeros(length(basis),1);
% ca = zeros(length(basis),1);
% bo = zeros(length(basis),1);
% pb = zeros(length(basis),1);
% 
% err = double(int(diff(fl,n)*diff(fl,n)*rho,-inf,l) + ...
%     int(diff(fr,n)*diff(fr,n)*rho,r,inf) + ...
%     int(diff(f,n)*diff(f,n)*rho,l,r));
% 
% x = -3:.01:3;
% for i = 1:length(basis);
%     
%     I = sym(0);
%     I = I + int(simplify(fl*rho*basis(i)),-inf,l);
%     I = I + int(f*rho*basis(i),l,r);
%     I = I + int(simplify(fr*rho*basis(i)),r,inf);
%     c(i) = double(I)
%     
%     I = sym(0);
%     I = I + int(simplify(diff(fl,n)*rho*basis(i)),-inf,l);
%     I = I + int(diff(f,n)*rho*basis(i),l,r);
%     I = I + int(simplify(diff(fr,n)*rho*basis(i)),r,inf);
%     cn(i) = double(I)
%     
%     if (i >= (n + 1))
%         ca(i) = (q^2/(2*lambda))^(n/2)*sqrt(factorial(i-1-n)/factorial(i-1))*cn(i-n);
%         bo(i) = (q^2/(2*lambda))^(n/2)*sqrt(factorial(i-1-n)/factorial(i-1))*err;
%         pb(i) = (q^2/(2*lambda))^(n/2)*sqrt(factorial(n)/(i-1)^n)*err;
%     end
%     
%     
%     
%     % basis_i = matlabFunction(basis(i));
%     % h = @(x) basis_i(x).*rh(x).*g(x);
%     % c(i) = quad(h,-14,14)
% end