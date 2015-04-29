%% INPUT VARIABLES OF THE PROGRAM
% Non-linearity in the equation

function nonlin = F(f)
    nonlin = 0.5*diff(f^2);
end

% Differential operator associated with the equation. 
function Au = oper(u)
    %Au = -diff(u,4) -diff(u,2);
    Au = diff(u,2) + u;
end

% Normalized eigenfunctions of the differential operator.
function f = eigfun(i,x)
    global pi;
    if mod(i,2) == 0;
        f = cos(i/2*x)*sqrt(1/pi);    
    else
        f = sin((i+1)/2*x)*sqrt(1/pi);
    end
   
    % f = sin(i*x)*sqrt(2/pi);
end

% Functions involved in the expansion ( = non-normalized eigenfunctions)
function f = expa(i,x)
    if mod(i,2) == 0;
        f = cos(i/2*x);
    else
        f = sin((i+1)/2*x);
    end
    
    % f = sin(i*x);
end

% Assumed to be constant
function q = noiseTerm(i)
    q = 1/i;
end

function initParams()
    global degree; degree = 2;
    global L; L = -pi;
    global R; R = pi;
    global N; N = 5;
    global d; d = 2;
end
