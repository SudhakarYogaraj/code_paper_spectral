function build_problem_nd()

% Format of output
format longEng

% Symbolic treatment of pi
pi = sym('pi','real');

%% USER INPUT %%

% Slow variables
ns = 1;

% Fast variables
nf = 2;

% Creation of symbolic variables
x  = sym(zeros(1, ns));
y  = sym(zeros(1, nf));
g  = sym(zeros(1, ns));
f  = sym(zeros(1, ns));
h  = sym(zeros(1, nf));
vy = sym(zeros(1, nf));
hy = sym(zeros(nf, nf));
gx = sym(zeros(ns, ns));
fx = sym(zeros(ns, ns));

for k = 1:nf; x(k) = sym(sprintf('x%d', k-1), 'real'); end
for k = 1:nf; y(k) = sym(sprintf('y%d', k-1), 'real'); end

% Potential
v = y(1)^2 + y(2)^2

% Coefficient of the BM
s = sqrt(2);

% Solution of the cell problem
g(1) = cos(x(1)) * sin(y(1));
g(2) = cos(x(1)) * sin(y(2));

% Non-leading order drift of fast process
h(1) = cos(x(1)) * cos(y(1)) * cos(y(2));
h(2) = cos(x(1)) * cos(y(1) + y(2));

%% DEPENDENT VARIABLES

% Square of covariance matrix
S = s*s

% Derivative of the potential
for i = 1:nf;
    vy(i) = diff(v,y(i));
end

% Associated invariant measure
rho = exp(-v)
rho_n = matlabFunction(rho);
rho = rho/integral2(rho_n, -inf, inf, -inf, inf);

% Generator in weighted space
Lw = @(f) 0.5 *( diff( S * rho * diff(f,y(1)) , y(1)) + diff( S * rho * diff(f,y(2)) , y(2)) ) / rho;

% Linear term
lin = simplify( 1/4*S * (diff(v,y(1),2) + diff(v,y(2),2)) - 1/8*S*( diff(v,y(1))^2 + diff(v,y(2))^2 ) )

% derivative of h
for i = 1:nf
    for j = 1:nf
        hy(i,j) = diff(h(i),y(j));
    end
end

% Differential of g
for i = 1:ns
    for j = 1:ns
        gx(i,j) = diff(g(i),x(j));
    end
end

% Associated rhs
for i = 1:ns
    f(i) = - simplify( Lw(g(i)) )
end

% x-derivative of rhs
for i = 1:ns
    for j = 1:ns
        fx(i,j) = simplify( diff(f(i),x(j)) )
    end
end

%% GENERATION OF FILES FOR C++ PROGRAM %%

% Output to files
f0 = fopen('tmp/gx.gen','w');
f1 = fopen('tmp/v.gen','w');
f2 = fopen('tmp/vy.gen','w');
f3 = fopen('tmp/g.gen','w');
f4 = fopen('tmp/f.gen','w');
f5 = fopen('tmp/fx.gen','w');
f6 = fopen('tmp/h.gen','w');
f7 = fopen('tmp/hy.gen','w');
f8 = fopen('tmp/lin.gen','w');
f9 = fopen('tmp/rho.gen','w');

fprintf(f0, ccode(gx));
fprintf(f1, ccode(v));
fprintf(f2, ccode(vy));
fprintf(f3, ccode(g));
fprintf(f4, ccode(f));
fprintf(f5, ccode(fx));
fprintf(f6, ccode(h));
fprintf(f7, ccode(hy));
fprintf(f8, ccode(lin));
fprintf(f9, ccode(rho));

fclose(f0); system('echo -e "" >> tmp/gx.gen');
fclose(f1); system('echo -e "" >> tmp/v.gen');
fclose(f2); system('echo -e "" >> tmp/vy.gen');
fclose(f3); system('echo -e "" >> tmp/g.gen');
fclose(f4); system('echo -e "" >> tmp/f.gen');
fclose(f5); system('echo -e "" >> tmp/fx.gen');
fclose(f6); system('echo -e "" >> tmp/h.gen');
fclose(f7); system('echo -e "" >> tmp/hy.gen');
fclose(f8); system('echo -e "" >> tmp/lin.gen');
fclose(f9); system('echo -e "" >> tmp/rho.gen');

exit(0);
