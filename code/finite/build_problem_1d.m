function build_problem_1d()

format longEng

% Symbolic variables
syms y x pi;

% Potential
v = y*y + log(pi)/2.

% Derivative of the potential
vy = diff(v,y);

% Associated invariant measure
rho = exp(-v)

% Coefficient of the BM
s = sqrt(2); S = 2;

% Generator in weighted space
Lw = @(f) 0.5 * diff( S * rho * diff(f,y) , y) / rho;

% Solution of the cell problem
g = cos(x) * sin(y);

% Associated rhs
f = - Lw(g)

% x-derivative of rhs
fx = diff(f,x)

% non-leading order drift of fast process
h = cos(x) * cos(y);

% derivative of h
hy = diff(h,y);

% Linear term
lin = 1/4*S*diff(v,y,2) - 1/8*S*(diff(v,y)^2)

% Standard deviation of approximating gaussian
sigma = 1.2;

% Approximating gaussian
gaussian = 1/sqrt(2*sym(pi)*sigma^2) * exp(-y^2/(2*sigma^2))

% Output to files
f1 = fopen('tmp/v.gen','w');
f2 = fopen('tmp/vy.gen','w');
f3 = fopen('tmp/g.gen','w');
f4 = fopen('tmp/f.gen','w');
f5 = fopen('tmp/fx.gen','w');
f6 = fopen('tmp/h.gen','w');
f7 = fopen('tmp/hy.gen','w');
f8 = fopen('tmp/lin.gen','w');
f9 = fopen('tmp/rho.gen','w');

fprintf(f1, ccode(v));
fprintf(f2, ccode(vy));
fprintf(f3, ccode(g));
fprintf(f4, ccode(f));
fprintf(f5, ccode(fx));
fprintf(f6, ccode(h));
fprintf(f7, ccode(hy));
fprintf(f8, ccode(lin));
fprintf(f9, ccode(rho));

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
