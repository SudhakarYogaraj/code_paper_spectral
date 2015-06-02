function compute_hermite()

    % Maximal degree
    N = 40;

    % Symbolic variables for the computation
    syms x; syms pi;

    % Weight function
    rho = 1/sqrt(2*sym(pi))*exp(-x^2/2);

    % Usual basis of polynomials
    monomials = sym(zeros(N + 1,1));

    for i = 0:N
        monomials(i+1) = x^i; 
    end

    % Gram matrix
    G = int(monomials*transpose(monomials)*rho,x,-inf,inf);

    % Cholesky decomposition
    U = chol(G,'upper','nocheck');

    % Step to construct the basis
    B = transpose(inv(U));

    % Basis of Hermite polynomials
    basis = simplify(B*monomials);

    % Save basis in file
    save('cellprob_basis','basis');

    ccode(B)

end
