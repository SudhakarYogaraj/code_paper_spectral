function build_problem_1d()

format longEng

% Symbolic variables
syms y x pi;
x = sym('x','real');
pi = sym('pi','real');


%% USER INPUT %%

% Fast variables
nf = 2;

% 
y = sym( zeros(1, nf) );
for k = 1:nf; 
    y(k) = sym(sprintf('y%d', k-1), 'real'); 
end

% Potential
% v = y^4/4 - y^2/2;
v = y(1)^2 + y(2)^2

% Coefficient of the BM
s = sqrt(2);

%% END OF USER INPUT %%

% Square of covariant matrix
S = s*s

% Derivative of the potential
vy = diff(v,y);

% Associated invariant measure
rho = exp(-v)
rho_n = matlabFunction(rho);
rho = rho/integral(rho_n, -inf, inf);

% Generator in weighted space
Lw = @(f) 0.5 * diff( S * rho * diff(f,y) , y) / rho;

% Solution of the cell problem
g = cos(x) * sin(y);

% Differential of g
gx = diff(g,x);

% Associated rhs
f = - simplify( Lw(g) )

% x-derivative of rhs
fx = simplify( diff(f,x) )

% non-leading order drift of fast process
h = cos(x) * cos(y);

% derivative of h
hy = diff(h,y);

% Linear term
lin = simplify( 1/4*S*diff(v,y,2) - 1/8*S*(diff(v,y)^2) )

% Standard deviation of approximating gaussian
sigma = 1.2;

% Approximating gaussian
gaussian = 1/sqrt(2*sym(pi)*sigma^2) * exp(-y^2/(2*sigma^2))

%% Generation of files for C++ program %%

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
