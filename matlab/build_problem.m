function build_problem()

% Ensure that folder for temporary files exist
system('mkdir -p tmp');

% Format of output
format longEng

% Symbolic treatment of pi
pi = sym('pi');

%% USER INPUT %%

% Slow variables
ns = 2;

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
% v = y(1)^2/2 + y(2)^2/2 + y(3)^2/2
% v = y(1)^4/4 - y(1)^2/2 + (y(2)-2)^2/2 + y(3)^4
v = 0.5 * ((y(1) - 1)^4 + (y(2)-1)^2 + 0.2*(y(1) - 1)*(y(2) -1)); % + 2*(y(3) + 1)^2
% v = y(1)^2 + y(2)^2
% v = y(1)*y(1) + log(pi)/2.
% v = y(1)^4/4 - y(1)^2/2;

% Coefficient of the BM
s = sqrt(2);

% Solution of the cell problem
g(1) = cos(x(1)) * sin(y(1)*y(2));
g(2) = cos(x(2) + y(2)) * sin(y(1) + y(2));
% g(2) = cos(x(1)) * sin(y(2) + y(1));

% Non-leading order drift of fast process
h(1) = cos(x(1)) * cos(y(1))  * cos(y(2));
h(2) = cos(x(1) + x(2)) * cos(y(1) + y(2));
% h(3) = cos(x(3)) * cos(y(3));

%% DEPENDENT VARIABLES

% Square of covariance matrix
S = s*s

% Derivative of the potential
for i = 1:nf;
    vy(i) = diff(v,y(i));
end

% Associated invariant measure
rho = simplify (exp(-v))
rho_n = matlabFunction(rho);
if (nf == 1)
    rho = rho/integral(rho_n, -inf, inf);
elseif (nf == 2)
    rho = rho/integral2(rho_n, -inf, inf, -inf, inf);
elseif (nf == 3)
    rho = rho/integral3(rho_n, -inf, inf, -inf, inf, -inf, inf);
end

% Generator in weighted space
if (nf == 1)
    Lw = @(f) 0.5 * diff( S * rho * diff(f,y(1)) , y(1)) / rho;
elseif (nf == 2)
    Lw = @(f) 0.5 *( diff( S * rho * diff(f,y(1)) , y(1)) + diff( S * rho * diff(f,y(2)) , y(2)) ) / rho;
elseif (nf == 3)
    Lw = @(f) 0.5 *( diff( S * rho * diff(f,y(1)) , y(1)) + diff( S * rho * diff(f,y(2)) , y(2)) + diff( S * rho * diff(f,y(3)) , y(3)) ) / rho;
end

% Linear term
if (nf == 1)
    lin = simplify( 1/4*S*diff(v,y(1),2) - 1/8*S*(diff(v,y(1))^2) )
elseif (nf == 2)
    lin = simplify( 1/4*S * (diff(v,y(1),2) + diff(v,y(2),2)) - 1/8*S*( diff(v,y(1))^2 + diff(v,y(2))^2 ) )
elseif (nf == 3)
    lin = simplify( 1/4*S * (diff(v,y(1),2) + diff(v,y(2),2) + diff(v,y(3),2)) - 1/8*S*( diff(v,y(1))^2 + diff(v,y(2))^2 + diff(v,y(3))^2 ) )
end

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

stardivh = sym('0');
for i = 1:nf
    stardivh = stardivh + simplify(vy(i)*h(i) - hy(i,i));
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

%% VARIABLES FOR THE HMM

vyy  = sym(zeros(nf, nf));
fy   = sym(zeros(ns, nf));
drif = sym(zeros(1, 2*nf));
diffu = sym(zeros(1, 2*nf));

% Extension of the y vector
for k = 1:(2 * nf); y(k) = sym(sprintf('y%d', k-1), 'real'); end

% Creation of fy
for i = 1:ns
    for j = 1:nf
        fy(i,j) = diff(f(i),y(j));
    end
end

% Second derivative of potential
for i = 1:nf
    for j = 1:nf
        vyy(i,j) = diff(vy(i), y(j));
    end
end

drif(1:nf) = -vy;
drif(nf+1:end) = -y(nf+1:end)*vyy + h;

for i = 1:nf
    diffu(i) = s;
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
f10 = fopen('tmp/fy.gen','w');
f11 = fopen('tmp/drif.gen','w');
f12 = fopen('tmp/diff.gen','w');

fprintf(f0, ccode(gx));
fprintf(f1, ccode(v));
fprintf(f2, ccode(vy));
fprintf(f3, ccode(g));
fprintf(f4, ccode(f));
fprintf(f5, ccode(fx));
fprintf(f6, ccode(h));
fprintf(f7, ccode(stardivh));
fprintf(f8, ccode(lin));
fprintf(f9, ccode(rho));
fprintf(f10, ccode(fy));
fprintf(f11, ccode(drif));
fprintf(f12, ccode(diffu));

command_nf = ['sed -i "s/\(this->nf =\) [0-9]*/\1 ', num2str(nf), '/g" problem.init'];
command_ns = ['sed -i "s/\(this->ns =\) [0-9]*/\1 ', num2str(ns), '/g" problem.init'];
system(command_nf);
system(command_ns);

% Feature
split = fopen('tmp/split.gen','w');

% Auxiliary functions
for i = 1:ns
    fprintf(split, 'double a%d (vector<double> x, vector<double> y) {\n', i);
    fprintf(split, ccode(f(i)));
    fprintf(split, '\n}\n\n');
end

for i = 1:ns
    for j = 1:ns
        fprintf(split, 'double dxa%d%d (vector<double> x, vector<double> y) {\n', i, j);
        fprintf(split, ccode(fx(i,j)));
        fprintf(split, '\n}\n\n');
    end
end

% Function to init functions
fprintf(split, 'void Problem::init_functions() {\n\n');

for i = 1:ns
    fprintf(split, '    fsplit[%d] = a%d;\n', i-1, i);
end

fprintf(split, '\n');

for i = 1:ns
    for j = 1:ns
        fprintf(split, '    fxsplit[%d][%d] = dxa%d%d;\n', i-1, j-1, i, j);
    end
end

fprintf(split, '}');


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
fclose(f10); system('echo -e "" >> tmp/fy.gen');
fclose(f11); system('echo -e "" >> tmp/drif.gen');
fclose(f12); system('echo -e "" >> tmp/diff.gen');
fclose(split);

exit(0);
