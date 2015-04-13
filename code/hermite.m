% In the whole program, the processes are numbered from 1 and not 0.

%% MAIN PART
% Function computing the auxiliary system associated with the multiscale Burger
% equation. The eigenvector involved in the expansion are not normalized. 
% TODO : Adapt formulae of the coefficients of the effective equation to the
% case of non-normalized eigenfunctions.

function output()

%% Initialisation of vectors and parameters
    
    % Clearing previously defined variables, command windows,
    % and closing all open figures.    clear all; close all; clc;

    % Symbolic variable and truncation index
    syms epsilon nu real;
    global x; x = sym('x','real');
    global y; y = sym('y','real');
    global pi; pi = sym('pi','real');
    global alpha; alpha = sym('alpha','real');
    global beta; beta = sym('beta','real');
    
    % Truncation index = total number of processes that are kept
    global N; 
    
    % Number of slow processes
    global d;

    % Left and right limits of the domain
    global L;  
    global R;    
    global degree;
    
    % Initialization of user-defined parameters
    initParams();
    
    % Export data to C++ files
    command_nf = ['sed -i "s/\(this->nf =\) [0-9]*/\1 ', num2str(N - d), '/g" /home/urbain/mres/source/program/Problem_init.c'];
    command_d = ['sed -i "s/\(this->d =\) [0-9]*/\1 ', num2str(d), '/g" /home/urbain/mres/source/program/Problem_init.c'];
    system(command_nf); system(command_d);
    
    % Creating the symbolic variables
    xs = sym(zeros(1,N));
    for k=1:N;  xs(k) = sym(sprintf('x%d', k-1),'real'); end

    % Symbolic variables for the auxiliary fast processes
    zs = sym(zeros(1,N-d));
    for k=1:N-d;  zs(k) = sym(sprintf('z%d', k-1),'real'); end

    % Symbolic variables to match with the notations of the c++ program (see
    % end of this program.
    ys = sym(zeros(1, 2*(N-d)));
    for k=1:2*(N-d);  ys(k) = sym(sprintf('y%d', k-1),'real'); end

%% Projection of the nonlinearity on the eigenspaces
    
    % Eigenfunctions of the differential operator, and functions used for the
    % seris expansion of the solution.
    fs = sym(zeros(1,N));
    es = sym(zeros(1,N));
    u  = sym(0);
    for i = 1:N; 
        fs(i) = eigfun(i);
        es(i) = expa(i);
        u = u + xs(i)*es(i);
    end
    
    nonLin = simplify(F(u));
    tic; proj = simplify(innerProduct(nonLin, transpose(fs))/sqrt(pi)); toc
    
%% Creating the system of SDEs from the original SPDE

% Construction of the drift-term at the right-hand side of the original equation
    RHS = sym(zeros(1,N));
      
    % Contribution of the differential operator
    As = sym(zeros(1,N));

    for i = 1:N
        % Part of the RHS for the i-th unknown
        As(i) = -eigval(i)*xs(i)/epsilon^2;
    end

    % Combination of both
    for i = 1:N; RHS(i) = nu*xs(i) + As(i) + proj(i)/epsilon; end

    % Scaled drift of the equations corresponding with the fast processes.
    drif = simplify(epsilon^2*RHS(d+1:N)); % Called b in the article HMM
    
    % Additional step due to the fact that Matlab simplify does not work
    % as we expect at times
    drif = simplify(expand(drif));
    
    % Construction of the auxiliary system
    drif_aux = sym(zeros(1,2*(N-d)));

    % Decomposition of the drift, using the notations of HMM
    b0 = subs(drif,epsilon,0);
    b1 = simplify((drif - b0)/epsilon);
    
    % Construction of the gradient of b0
    gradb0 = sym(zeros(N-d,N-d));
    for i = 1:N-d
        for j = 1:N-d
            gradb0(i,j) = diff(b0(i),xs(j+d));
        end
    end

    drif_aux(1:N-d) = b0;
    drif_aux(N-d+1:end) = transpose(gradb0*transpose(zs)) + b1;

    % Invariant measure must be evaluated at epsilon = 0
    drif_aux = subs(expand(drif_aux), epsilon, 0);

    % Drift coefficient of the slow processes
    a = simplify(expand(epsilon*RHS(1:d)));

    % Separation of a in 2 terms according to the power of epsilon.
    a_aux = subs(a,epsilon,0);
    a_nu  = simplify((a - a_aux)/epsilon);
    a = a_aux; % for convenience in future notations

    % Computation of d_y a and d_x a   
    da = sym(zeros(d,N));

    for i = 1:d
        for j = 1:N
            da(i,j) = diff(a(i),xs(j));
        end
    end

    dax = da(:,1:d);
    day = da(:,d+1:N);
    dae = diff(a,epsilon);

    % Variables to export: drif_aux, dax, day  
    % Renaming of the variables to match C++ program
    dax = subsVar(dax,xs(d+1:N),ys(1:N-d));  
    dax = subsVar(dax,zs,ys(N-d+1:end));

    day = subsVar(day,xs(d+1:N),ys(1:N-d));  
    day = subsVar(day,zs,ys(N-d+1:end));

    drif_aux = subsVar(drif_aux,xs(d+1:N),ys(1:N-d));  
    drif_aux = subsVar(drif_aux,zs,ys(N-d+1:end));

    a = subsVar(a,xs(d+1:N),ys(1:N-d));  
    a = subsVar(a,zs,ys(N-d+1:end));

    diffu = zeros(1,2*(N-d));
    for i = 1:N-d; diffu(i) = noiseTerm(i+d); end

    f1 = fopen('build/a.gen','w');
    f2 = fopen('build/a_nu.gen','w');
    f3 = fopen('build/dax.gen','w');
    f4 = fopen('build/day.gen','w');
    f5 = fopen('build/drif.gen','w');
    f6 = fopen('build/diff.gen','w');
    f7 = fopen('build/drif_init.cpp','w');

    fprintf(f1,ccode(a));
    fprintf(f2,ccode(a_nu));
    fprintf(f3,ccode(dax));
    fprintf(f4,ccode(day));
    fprintf(f5,ccode(drif_aux));
    fprintf(f6,ccode(sym(diffu)));
    fprintf(f7,ccode(drif));

    fclose(f1); system('echo -e "" >> build/a.gen');
    fclose(f2); system('echo -e "" >> build/a_nu.gen');
    fclose(f3); system('echo -e "" >> build/dax.gen');
    fclose(f4); system('echo -e "" >> build/day.gen');
    fclose(f5); system('echo -e "" >> build/drif.gen');
    fclose(f6); system('echo -e "" >> build/diff.gen');
    fclose(f7); system('echo -e "" >> build/drif_init.cpp');

    
%% Constructing the Hermite polynomials    
    
    vars = xs(d+1:end);
    lambda = zeros(1,length(vars));
    q = zeros(1,length(vars));
    for i = d+1:length(xs)
        lambda(i-d) = eigval(i);
        q(i-d) = noiseTerm(i);
    end
    
    basis = sym(zeros(degree + 1, length(vars)));
    
    syms w;
    mon_w = sym(zeros(degree + 1,1));
    for i = 0:degree; mon_w(i+1) = w^i; end   
    
    syms pi real positive
    
    for i = 1:length(vars)
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
            
            % Adaptation in order to have monomials, which will be more
            % convenient to work with at a later stage. As a consequence,
            % the resulting polynomials will not be normalized anymore (but
            % they will remain orthogonal).
            B = diag(1./diag(B))*B;
            
            % First eigenfunctions of the generator of the OU process with 
            % parameters lambda(i) and q(i))
            basis(:,i) = simplify(B*mon_w);
            
            % Check of the orthogonality
            simplify(int(basis(:,i)*transpose(basis(:,i))*rho,w,-inf,inf));
        else
            
            % In this case, the eigenfunctions of the operator are all the
            % polynomials with only one term.
            basis(:,i) = mon_w;
        end
        
        % Substitution of w by the variable
        basis(:,i) = subs(basis(:,i), w, vars(i));
    end

%% Solving the cell problem
    

    % p = x  +  y + y^3  + x^3 + x*y + x*y^2 + y^4 + x^4*y + x^5 + 1 + z*z + ...
    %    z*x*y + z*x + z*x^2 + x^2*2*y + z*z*3.1;
    
    phi = sym(zeros(d,1));
    for k = 1:d
        p = proj(k);
        mat = poly2mat(p,vars,basis);

        % Check that the polynomial can be recovered from the matrix.
        if(simplify(p - mat2poly(mat,basis)) ~= 0)
            error('poly2mat did not work');
        end

        eigsGen = zeros(size(mat))+inf;
        for i = 1:numel(mat)
            i/numel(mat)
            if (mat(i) ~= 0)
                sub = myInd2Sub(size(mat),i);
                p_i = sym(1);
                for j = 1:length(basis(1,:))
                    p_i = p_i * basis(sub(j),j);
                end
                gp_i = simplify(gen(p_i,vars));
                eigsGen(i) = simplify(gp_i/p_i);
            end
            
            if(abs(eigsGen(i)) < 10*eps)
               error('Centering condition not satisfied');
            end
        end

        phi(k) = -mat2poly(mat./eigsGen,basis);
        simplify(gen(phi(k),vars) + p)
    end

%% Calculation of the coefficients of the effective equation

    % Vector of varibables of the original system
    z = [xs(1:d), ys(1:N-d)];

    % Part of the drift of the fast process that is due to the nonlinearity
    b2 = subs(b1,epsilon,0);
    b2 = subsVar(b2,xs(d+1:end),ys(1:N-d));
    
    % Adaptation of the variable name for phi
    phi = subsVar(phi,xs(d+1:end),ys(1:N-d));
    
    % Initialization of the (Square of) the diffusion coefficient of the
    % effective equation.
    g = 2*simplify(phi*a);
    
    % Initialization of the drift coefficient of the effective equation
    h = sym(zeros(1,d));
    for i = 1:d
        for j = 1:d
            h(i) = h(i) + diff(phi(i),xs(j))*a(j);
        end
    end
    
    for i = 1:d;
        for j = 1:N-d
            h(i) = h(i) + diff(phi(i),ys(j))*b2(j);
        end
    end
    
    % Integration with respect to the invariant measure to obtain the
    % coefficients
    for i = 1:N-d
        l = eigval(i+d);
        q = noiseTerm(i+d);
        if(abs(q)>1e-10)
            rho = sqrt(l/(pi*q^2))*exp(-l*ys(i)^2/q^2);
            g = simplify(int(g*rho,ys(i),-inf,inf));
            h = simplify(int(h*rho,ys(i),-inf,inf));
        else
            g = subs(g,ys(i),0);
            h = subs(h,ys(i),0);
        end
        g
        h
    end
    
    g = 0.5*(g+transpose(g));
    
    g = ccode(chol(g,'lower','noCheck'));
    h = ccode(h);

    f1 = fopen('build/solDiff.gen','w');
    f2 = fopen('build/solDrif.gen','w');
    
    fprintf(f1,g);
    fprintf(f2,h);

    fclose(f1); system('echo -e "" >> build/solDiff.gen');
    fclose(f2); system('echo -e "" >> build/solDrif.gen');
    
exit(0);
end

%% Auxiliary functions

% Function that returns the matrix of coefficients of a polynomial. The
% last argument contains the basis functions associated with each of the
% variable. Those have to be monomials for the program to work. The
% multivariate basis of polynomials is obtain by all the possible products
% of the elements of each column. 
%
% E.g. if p is a polynomial in 3 variables, x, y, z, and we take the
% classical monomial basis functions, then the matrix basis would look
% like:  [1   1   1  ]
%        [x   y   z  ]
%        [x^2 y^2 z^2]
% Furthermore the element (i,j,k) of the resulting matrix would be equal to
% the coefficient x^i y^j z^k of the polynomial.
function mat = poly2mat(p,vars,basis)
    degree = length(basis(:,1)) - 1;
    matOld = p;
    for i = 1:length(vars)
        si   = (degree + 1)*ones(1,i);
        mat = sym(zeros([si 1]));
        for j = 1:numel(matOld)     
            j/numel(matOld)
            sub = myInd2Sub(si(1:end-1),j);
            degree_j = degree - sum(sub - 1);
            vec = poly2vec(matOld(j), vars(i), basis(1:degree_j + 1,i));
            for k = 1:degree_j + 1;
                indNew = mySub2Ind(si, [sub k]);
                mat(indNew) = vec(k);
            end

        end
        matOld = mat;
    end
end

% Reverse function to the previous one. Starting from the matrix of
% coefficients of a polynomial in a base specified by the last argument,
% the symbolic expression is obtained.
function p = mat2poly(mat,basis)
    p = sym(0); siz = size(mat);
    for i = 1:numel(mat)
        i/numel(mat)
        pi = sym(mat(i));
        if (pi == 0); continue; end
        sub = myInd2Sub(siz,i);    
        for j = 1:length(basis(1,:))
            pi = pi * basis(sub(j),j);
        end
        p = p + pi;
    end
    p = simplify(p);
end

% Get subscript from linear index
function sub = myInd2Sub(siz, ind)
    sub = zeros(1,length(siz));
    k = [1 cumprod(siz(1:end-1))];
    for i = length(siz):-1:1,
        vi = rem(ind-1, k(i)) + 1;
        vj = (ind - vi)/k(i) + 1;

        sub(i) = double(vj);
        ind = vi;
    end
end

% Get linear index from subscript.
function ind = mySub2Ind(siz, sub)
    k = [1 cumprod(siz(1:end-1))];
    ind = 1;
    for i = 1:length(sub)
        v = sub(i);
        ind = ind + (double(v)-1)*k(i);
    end
end

% Function to obtain the coefficients of a polynomial in one variable, in
% the basis of monomials specified in basis_mon (by increasing degrees);
function vec = poly2vec(p,var,basis_mon)

    degree = length(basis_mon) - 1;
    vec = sym(zeros(degree + 1,1));
    for i = degree:-1:0
        aux = expand(p/var^i);
        vec(i+1) = subs(aux,var,inf);
        p = simplify(var^i*(aux) - vec(i+1)*basis_mon(i+1));
    end
end

% Generator of the stochastic partial differential equation
function G = gen(f,ys)

    G = 0; global d;
    for i = 1:length(ys)
        G = G - ys(i)*eigval(i+d)*diff(f,ys(i)) + ...
            0.5*noiseTerm(i+d)^2*diff(f,ys(i),2);
    end
    % !!!!!! noiseTerm i+d or just i????
end

% Substitute variable x by y in symbolic function f
function f = subsVar(f,x,y)
    N = length(x);
    for i = 1:N
        f = subs(f,x(i),y(i));
    end
end

% (Opposite of the) eigenvalues of the differential operator.
function li = eigval(i)
    global x;
    f  = eigfun(i);
    li = - oper(f)/f;
end

% Normalized eigenfunctions of the differential operator.
function f = eigfun(i)
    f_aux = expa(i);
    f = f_aux/sqrt(innerProduct(f_aux,f_aux));
end

function ip = innerProduct(u,v)
    global x;
    global y;
    global L;
    global R;
    
    ip = int(int(u*v,x,L,R),y,L,R);
end

%% INPUT VARIABLES OF THE PROGRAM
% Non-linearity in the equation

function nonlin = F(f)
    global x y;
    nonlin = 0.5*diff(f,x)*diff(f,y);
end

% Differential operator associated with the equation. 
function Au = oper(u)
    global x y;
    Au = diff(u,x,2) + diff(u,y,2) + 2*u;
end

% Functions involved in the expansion ( = non-normalized eigenfunctions)
function f = expa(I)
    global x y;
    
    [i,j] = meshgrid(1:I);
    lambda = i.^2 + j.^2;
    
    for loop = 1:I
        [mm, ii] = min(lambda);
        [~, jj] = min(mm);
        lambda(ii(jj),jj) = inf;
    end

    f = sin(ii(jj)*x)*sin(jj*y);
end

% Assumed to be constant
function q = noiseTerm(i)
    q = 1/i*(i<=5); %(i == 2);
end

function initParams()
    global degree L R N d;
    degree = 2;
    L = 0;
    R = pi;
    N = 4;
    d = 1;
end
