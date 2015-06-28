/* TODO: Problem with exact sigma (urbain, Wed 20 May 2015 21:00:20 BST) */

#include "structures.hpp"
#include "toolbox.hpp"
#include "Problem.hpp"
#include "Gaussian_integrator.hpp"
#include "Solver_spectral.hpp"
#include "templates.hpp"
#include "io.hpp"
#include "global.hpp"
#include <iomanip>

using namespace std;

/*! Function that calculates the 0-order term in the Schrodinger equation
 *
 * Assuming that L is the generator of an OU process with
 *
 * TODO: add description (urbain, Thu 25 Jun 2015 01:27:09 CEST)
 */
double Solver_spectral::gaussian_linear_term(vector<double> z) {

    // Laplacian term.
    double laplacian = 0.;

    // Square of the gradient.
    double grad2 = 0.;

    // Square of the coefficient of the noise.
    double S = 2.;

    // Computation of the term
    for (int k = 0; k < this->nf; ++k) {
        laplacian += 1/this->eig_val_cov[k];
        grad2 += z[k]*z[k] /this->eig_val_cov[k];
    }

    // Standard formula for the potential
    return (0.25 * S * laplacian - 0.125 * S * grad2);
}

/*! Change of variable y = Cz + m
 *
 * This function implements the change of variables y = Cz + m.. If z is
 * distributed as G(0,I), y will be distributed as G(CC^T,m).
 */
vector<double> Solver_spectral::map_to_real(vector<double> z) {

    // Initialization
    vector<double> result(z.size(), 0.);

    for (int i = 0; i < z.size(); ++i) {

        // Left-multiply z by covariance matrix
        for (int j = 0; j < z.size(); ++j) {
            result[i] += this->sqrt_cov[i][j] * z[j];
        }

        // Add bias
        result[i] += this->bias[i];
    }

    return result;
}

/*! Update the statistics of the Gaussian related to Hermite functions
 *
 * This function updates the statistics of the Gaussian from which the Hermite
 * functions are calculated. In simple words, the Hermite functions are
 * concentrated where the density of the Gaussian is non-zero.
 */
void Solver_spectral::update_stats(Problem &problem, vector<double> var_scaling) {

    // Update of bias and covariance
    this->bias = problem.bias;
    this->eig_vec_cov = problem.eig_vec_cov;
    this->eig_val_cov = problem.eig_val_cov;

    // Apply user-defined extra-scaling
    for (int i = 0; i < nf; ++i) {
        this->eig_val_cov[i] *= var_scaling[i]*var_scaling[i];
    }

    // Square root of covariance matrix
    this->sqrt_cov = this->eig_vec_cov;
    for (int i = 0; i < nf; i++) {
        for (int j = 0; j < nf; j++) {
            this->sqrt_cov[i][j] *= sqrt(this->eig_val_cov[j]);
        }
    }

    // Determinant of the covariance matrix
    this->det_cov = 1.;
    for (int i = 0; i < nf; ++i) {
        this->det_cov *= this->eig_val_cov[i];
    }
}

vector< vector<double> >  Solver_spectral::discretize_a(Problem &problem, vector<double> x, Gaussian_integrator& gauss) {

    // Weights and points of the integrator.
    vector< vector<double> > quad_points = gauss.nodes;
    vector<double>  quad_weights = gauss.weights;

    // Parameters of the problem
    int ns = problem.ns;
    int ni = gauss.weights.size();

    // Discretization of a at gridpoints
    vector< vector<double> > a_discretized(ns, vector<double> (ni, 1.));

    // Computation of the function to integrate at the gridpoints.
    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ni; ++j) {

            // Integration point and rescaled version for integrattion
            vector<double> z = quad_points[j];
            vector<double> y = map_to_real(z);

            // Scaling to pass to Schrodinger equation;
            a_discretized[i][j] *= sqrt(problem.rho(x,y));

            // Scaling needed for the integration;
            a_discretized[i][j] /= sqrt(gaussian(z));

            // Discretization of a
            a_discretized[i][j] *= problem.a(x,y)[i];

            // Weight of the integration
            a_discretized[i][j] *= quad_weights[j];
        }
    }

    return a_discretized;
}

vector<double> Solver_spectral::discretize(Problem &problem, vector<double> x, Gaussian_integrator& gauss, double(*f)(vector<double>, vector<double> )) {

    // Weights and points of the integrator.
    vector< vector<double> > quad_points = gauss.nodes;
    vector<double>  quad_weights = gauss.weights;

    // Number of integration points
    int ni = gauss.weights.size();

    // Discretization of a at gridpoints
    vector<double> f_discretized(ni, 1.);

    // Computation of the function to integrate at the gridpoints.
    for (int j = 0; j < ni; ++j) {

        // Integration point and rescaled version for integrattion
        vector<double> z = quad_points[j];
        vector<double> y = map_to_real(z);

        // Scaling to pass to Schrodinger equation;
        f_discretized[j] *= sqrt(problem.rho(x,y));

        // Scaling needed for the integration;
        f_discretized[j] /= sqrt(gaussian(z));

        // Discretization of a
        f_discretized[j] *= f(x,y);

        // Weight of the integration
        f_discretized[j] *= quad_weights[j];
    }

    return f_discretized;
}

vector<double> Solver_spectral::project(int nf, int degree, Gaussian_integrator& gauss, vector<double> f_discretized) {

    // Number of polynomials
    int nb = bin(degree + nf, nf);

    // Number of integration points
    int ni = gauss.weights.size();

    // Nodes of the gaussian quadrature
    vector< vector<double> > quad_points = gauss.nodes;

    // Vector of coefficients of the projection
    vector<double> coefficients(nb, 0.);

    // Loop over all the monomials
    for (int i = 0; i < nb; ++i) {

        // Multi-index associated with i
        vector<int> m = ind2mult(i, nf);

        // Loop over the points of the quadrature
        for (int j = 0; j < ni; ++j) {

            // Evaluate monomial in point
            double mon_val = monomial(m, quad_points[j]);

            // Update coefficient
            coefficients[i] += f_discretized[j] * mon_val;
        }

        // Scaling due to change of variable
        coefficients[i] *= sqrt(sqrt(det_cov));
    }
    return coefficients;
}

SDE_coeffs Solver_spectral::estimator(Problem &problem, vector<double> x, double t) {

    // Vectors to store the coefficients of the sde
    SDE_coeffs sde_coeffs;

    // Integrator
    Gaussian_integrator gauss = Gaussian_integrator(nNodes,nf);

    // Parameters of the problem
    int nf = problem.nf;
    int ns = problem.ns;
    int nb = bin(degree + nf, nf);

    // User defined rescaling of the eigenvalues
    vector<double> var_scaling = {0.8, 0.7};

    // Update statistics of Gaussian
    this->update_stats(problem, var_scaling);

    // Discretization of a at gridpoints
    vector< vector<double> >a_discretized =  discretize_a(problem, x, gauss);

    // Weights and points of the integrator.
    vector< vector<double> > quad_points = gauss.nodes;

    // Coefficients of the projection of a.
    vector< vector<double> > coefficients(ns);

    // Expansion of right-hand side of the Poisson equation
    for (int i = 0; i < ns; ++i) {
        coefficients[i] = project(nf, degree, gauss, a_discretized[i]);
    }

    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));
    vector<double> coefficients_h(nb, 0.);
    if(VERBOSE) cout << "* Calculating the right-hand side of the linear system." << endl << endl;
    for (int i = 0; i < nb; ++i) {
        if(VERBOSE) progress_bar(((double) i )/((double) nb));
        vector<int> multIndex = ind2mult(i, nf);

        vector<double> v0(ns,0.);
        vector< vector<double> > m0(ns, vector<double> (ns,0.));

        auto lambda_dx = [&] (vector<double> z) -> vector< vector<double> > {
            vector<double> y = this->map_to_real(z);
            return problem.dax(x,y) * monomial(multIndex, z) * sqrt( problem.rho(x,y) / gaussian(z) );
        };
        auto lambda_h = [&] (vector<double> z) -> double {
            vector<double> y = this->map_to_real(z);
            return problem.stardiv_h(x,y) * monomial(multIndex, z) * sqrt( problem.rho(x,y) / gaussian(z) );
        };
        // relation stardiv hmm formula.

        // ! Hermite polynomials don't have scaling (H(x,s) = H(x/s, 1))
        vector< vector<double> > result_dx = gauss.quadnd(lambda_dx, m0) * sqrt(sqrt(det_cov));
        double result_h = gauss.quadnd(lambda_h) * sqrt(sqrt(det_cov));;

        for (int j = 0; j < ns; ++j) {
            for (int k = 0; k < ns; ++k) {
                coefficients_dx[j][k][i] = result_dx[j][k];
            }
        }
        coefficients_h[i] = result_h;
    }
    if(VERBOSE) end_progress_bar();

    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j) {
            coefficients_dx[i][j] = basis2herm(coefficients_dx[i][j],nf,degree);
        }
        coefficients[i] = basis2herm(coefficients[i],nf,degree);
    }
    coefficients_h = basis2herm(coefficients_h,nf,degree);

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nb,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));

    int n_tmp  = bin(2*degree + nf, nf);
    int n_done = 0;
    vector<double> tmp_vec(n_tmp, 0.);
    vector<int> position_visited(n_tmp, 0);
    if(VERBOSE) cout << "* Calculating the necessary products < L mi , mj >." << endl << endl;
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i, nf);
        for (int j = 0; j < nb; ++j) {
            vector<int> m2 = ind2mult(j, nf);
            int index = mult2ind(m1 + m2);
            if (position_visited[index] == 0) {
                n_done ++;
                if(VERBOSE) progress_bar(((double) n_done)/((double) n_tmp));
                position_visited[index] = 1;
                auto lambda = [&] (vector<double> z) -> double {
                    vector<double> y = map_to_real(z);
                    double tmp_term = (gaussian_linear_term(z) - problem.linearTerm(x,y));
                    /* if (fabs(tmp_term) > 1e-7) { */
                    /*     cout << tmp_term; */
                    /* } */
                    return tmp_term * monomial(m1 + m2, z);
                };
                tmp_vec[index] += gauss.quadnd(lambda);
            }
        }
    }
    if(VERBOSE) end_progress_bar();

    vector< vector<double> > prod_mat(nb, vector<double>(nb, 0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i, nf);
        for (int j = 0; j < nb; ++j) {
            vector<int> m2 = ind2mult(j, nf);
            int index = mult2ind(m1 + m2);
            prod_mat[i][j] = tmp_vec[index];
        }
    }

    if(VERBOSE) cout << "* Creating the matrix of the linear system." << endl << endl;
    vector< vector<double> > tmp_mat(nb, vector<double>(nb, 0.));
    vector< vector<double> > mat(nb, vector<double>(nb, 0.));
    for (int i = 0; i < nb; ++i) {
        if(VERBOSE) progress_bar(( (double) (i*nb) )/( (double) (nb*(2*nb + 1)) ));
        for (int j = 0; j < nb; ++j) {
            for (int k = 0; k < nb; ++k) {
                tmp_mat[i][j] += hermiteCoeffs_nd[j][k]*prod_mat[i][k];
            }
        }
    }
    for (int i = 0; i < nb; ++i) {
        if(VERBOSE) progress_bar(( (double) (nb*nb + i*(i+1)) )/( (double) (nb*(2*nb+1)) ));
        vector<int> m1 = ind2mult(i, nf);
        /* FIXME: sigmas? (urbain, Thu 28 May 2015 18:15:19 BST) */
        /* FIXME: Normalization (urbain, Thu 28 May 2015 20:55:46 BST) */

        for (int j = 0; j < nf; ++j) {
            mat[i][i] += m1[j]/this->eig_val_cov[j];
        }
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k < nb; ++k) {
                mat[i][j] += hermiteCoeffs_nd[i][k]*tmp_mat[k][j];
            }
            mat[j][i] = mat[i][j];
        }
    }
    if(VERBOSE) end_progress_bar();

    if(VERBOSE) cout << "* Solving linear system." << endl;
    for (int i = 0; i < ns; ++i) {
        solution[i] = solve(mat, coefficients[i]);
        for (int j = 0; j < ns; ++j) {
            solution_dx[i][j] = solve(mat, coefficients_dx[i][j]);
        }
    }
    if(VERBOSE) cout << "done!" << endl;

    // Calculation of the coefficients of the simplified equation
    vector<double> F1(ns, 0.);
    vector<double> F2(ns, 0.);
    vector< vector <double> > A0(ns, vector<double>(ns,0.));

    // Calculation of the coefficients of the effective equation
    for (int i = 0; i < nb; ++i) {

        // First part of the drift coefficient
        for (int j = 0; j < ns; ++j) {
            for (int k = 0; k < ns; ++k)
                F1[j] += solution_dx[j][k][i]*coefficients[k][i];
        }

        // Second part of the drift coefficient
        for (int j = 0; j < ns; ++j) {
            F2[j] += solution[j][i]*coefficients_h[i];
        }

        // Diffusion coefficient
        for (int j = 0; j < ns; ++j) {
            for (int k = 0; k < ns; ++k) {
                A0[j][k] += 2*solution[j][i]*coefficients[k][i];
            }
        }
    }

    sde_coeffs.diff =  cholesky(symmetric(A0));
    sde_coeffs.drif = F1 + F2;

    return sde_coeffs;
}

/*! Constructor of the spectral solver
 *
 * The contructor takes 3 arguments:
 * - The degree of the polynomial approximation to use,
 * - The number of nodes to use in the Gauss-Hermite quadrature,
 * - The number of variables in the fast processes.
 */
Solver_spectral::Solver_spectral(int degree, int nNodes, int n_vars) {

    // Number of fast processes
    this->nf = n_vars;

    // Dimension of the approximation space
    int nb = bin(degree + n_vars, n_vars);

    // Matrices that will contain the coefficients of the 1- and multi- dimensional
    // Hermite polynomials
    vector< vector<double> > mat1d (degree + 1, vector<double>(degree + 1,0.));
    vector< vector<double> > matnd(nb, vector<double>(nb,0.));

    // Fill the matrix for the unidimensional coefficients.
    hermite_coefficients(degree, mat1d);

    // Calculate the multi-dimensional coefficients by tensor products.
    for (int i = 0; i < nb; ++i) {

        // Multi-index corresponding to linear index i
        vector<int> m1 = ind2mult(i,n_vars);

        for (int j = 0; j < nb; ++j) {

            // Multi-index corresponding to linear index j
            vector<int> m2 = ind2mult(j,n_vars);

            // Calculation of the coefficients of x^m2 in h^m1
            matnd[i][j] = 1.;

            for (int k = 0; k < n_vars; ++k) {
                matnd[i][j] *= mat1d[m1[k]][m2[k]];
            }
        }
    }

    // Initialize variables of the solver
    this->degree = degree;
    this->nNodes = nNodes;
    this->hermiteCoeffs_1d = mat1d;
    this->hermiteCoeffs_nd = matnd;
}

/*! Compute linear index from multi-index
 *
 * This functions computes the linear index associated with a given multi-index.
 * Note that the ordering chosen is such that the resulting multi-index does not
 * depend on the degree of polynomial approximation.
 */
int Solver_spectral::mult2ind(vector<int> m) {

    // Number of variables
    int n = m.size();

    // Univariate case
    if(n == 1) {
        return m[0];
    }

    // Degree of m
    int degree = 0;

    // Calculation of the degree
    for (int i = 0; i < m.size(); ++i) {
        degree += m[i];
    }

    // Case degree == 0
    if (degree == 0) {
        return 0;
    }

    // Dimension of the space of degree (degree - 1)
    int result = bin((degree-1) + n, n);

    // Remove last element of m
    m.pop_back();

    return result + mult2ind(m);
}

/*! Compute multi-index from linear index
 *
 * This function computes the multi-index corresponding to the linear index
 * passed in argument. The arguments are
 * - ind: the linear index,
 * - n: size of the multiindex.
 */
vector<int> Solver_spectral::ind2mult(int ind, int n) {

    // Univariate case
    if (n == 1) {
        return vector<int>(1,ind);
    }

    // Degree of the corresponding multi-index
    int degree = 0;

    // Calculation of the degree
    while ( bin(degree + n, n) <= ind) {
        degree ++;
    }

    // Case when degree == 0
    if (degree == 0) {
        return vector<int>(n,0);
    }

    // Dimension of the space of strictly lower degree
    int dim = bin((degree-1) + n, n);

    // Initialization of the vector that will contain the result
    vector<int> result = ind2mult(ind - dim, n-1);

    // Add an element at the end
    result.push_back(degree);

    // Calculate value of last element
    for (int i = 0; i < n-1; ++i) {
        result[n-1] -= result[i];
    }

    return result;
}

/*! Function to evaluate monomial in a point
 *
 * This function returns the evaluation of a monomial, parametrized by its
 * multi-index, in a point.
 */
double Solver_spectral::monomial(vector<int> mult, vector<double> x) {

    double result = 1.;
    for (unsigned int i = 0; i < mult.size(); ++i) {
        result *= ipow(x[i], mult[i]);
    }
    return result;
}

vector<double> Solver_spectral::basis2herm (vector<double> bcoeffs, int n, int ns) {
    vector<double> result(bcoeffs.size(), 0.);
    for (unsigned int i = 0; i < bcoeffs.size(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            result[i] += hermiteCoeffs_nd[i][j] * bcoeffs[j];
        }
    }
    return result;
}

/*! Function to compute hermite coefficients
 *
 * This fills the matrix passed in argument with the hermite coefficients of
 * degree 0 to the degree passed to the method. The coefficients on the first
 * line represent the coefficients of the.
 */
void Solver_spectral::hermite_coefficients (int degree, vector< vector<double> >& matrix) {

    if (degree >= 0) {
        matrix[0][0] = 1.0;
    }

    if (degree >= 1) {
        matrix[1][1] = 1.0;
    }

    if (degree >= 2) {
        matrix[2][0] = sqrt(2.0)*(-1.0/2.0);
        matrix[2][2] = sqrt(2.0)*(1.0/2.0);
    }

    if (degree >= 3) {
        matrix[3][1] = sqrt(6.0)*(-1.0/2.0);
        matrix[3][3] = sqrt(6.0)*(1.0/6.0);
    }

    if (degree >= 4) {
        matrix[4][0] = sqrt(6.0)*(1.0/4.0);
        matrix[4][2] = sqrt(6.0)*(-1.0/2.0);
        matrix[4][4] = sqrt(6.0)*(1.0/1.2E1);
    }

    if (degree >= 5) {
        matrix[5][1] = sqrt(3.0E1)*(1.0/4.0);
        matrix[5][3] = sqrt(3.0E1)*(-1.0/6.0);
        matrix[5][5] = sqrt(3.0E1)*(1.0/6.0E1);
    }

    if (degree >= 6) {
        matrix[6][0] = sqrt(5.0)*(-1.0/4.0);
        matrix[6][2] = sqrt(5.0)*(3.0/4.0);
        matrix[6][4] = sqrt(5.0)*(-1.0/4.0);
        matrix[6][6] = sqrt(5.0)*(1.0/6.0E1);
    }

    if (degree >= 7) {
        matrix[7][1] = sqrt(3.5E1)*(-1.0/4.0);
        matrix[7][3] = sqrt(3.5E1)*(1.0/4.0);
        matrix[7][5] = sqrt(3.5E1)*(-1.0/2.0E1);
        matrix[7][7] = sqrt(3.5E1)*(1.0/4.2E2);
    }

    if (degree >= 8) {
        matrix[8][0] = sqrt(7.0E1)*(1.0/1.6E1);
        matrix[8][2] = sqrt(7.0E1)*(-1.0/4.0);
        matrix[8][4] = sqrt(7.0E1)*(1.0/8.0);
        matrix[8][6] = sqrt(7.0E1)*(-1.0/6.0E1);
        matrix[8][8] = sqrt(7.0E1)*5.952380952380952E-4;
    }

    if (degree >= 9) {
        matrix[9][1] = sqrt(7.0E1)*(3.0/1.6E1);
        matrix[9][3] = sqrt(7.0E1)*(-1.0/4.0);
        matrix[9][5] = sqrt(7.0E1)*(3.0/4.0E1);
        matrix[9][7] = sqrt(7.0E1)*(-1.0/1.4E2);
        matrix[9][9] = sqrt(7.0E1)*1.984126984126984E-4;
    }

    if (degree >= 10) {
        matrix[10][0] = sqrt(7.0)*(-3.0/1.6E1);
        matrix[10][2] = sqrt(7.0)*(1.5E1/1.6E1);
        matrix[10][4] = sqrt(7.0)*(-5.0/8.0);
        matrix[10][6] = sqrt(7.0)*(1.0/8.0);
        matrix[10][8] = sqrt(7.0)*(-1.0/1.12E2);
        matrix[10][10] = sqrt(7.0)*1.984126984126984E-4;
    }

    if (degree >= 11) {
        matrix[11][1] = sqrt(7.7E1)*(-3.0/1.6E1);
        matrix[11][3] = sqrt(7.7E1)*(5.0/1.6E1);
        matrix[11][5] = sqrt(7.7E1)*(-1.0/8.0);
        matrix[11][7] = sqrt(7.7E1)*(1.0/5.6E1);
        matrix[11][9] = sqrt(7.7E1)*(-9.920634920634921E-4);
        matrix[11][11] = sqrt(7.7E1)*1.803751803751804E-5;
    }

    if (degree >= 12) {
        matrix[12][0] = sqrt(2.31E2)*(1.0/3.2E1);
        matrix[12][2] = sqrt(2.31E2)*(-3.0/1.6E1);
        matrix[12][4] = sqrt(2.31E2)*(5.0/3.2E1);
        matrix[12][6] = sqrt(2.31E2)*(-1.0/2.4E1);
        matrix[12][8] = sqrt(2.31E2)*(1.0/2.24E2);
        matrix[12][10] = sqrt(2.31E2)*(-1.984126984126984E-4);
        matrix[12][12] = sqrt(2.31E2)*3.006253006253006E-6;
    }

    if (degree >= 13) {
        matrix[13][1] = 5.479963503528103E1*(1.0/3.2E1);
        matrix[13][3] = 5.479963503528103E1*(-1.0/1.6E1);
        matrix[13][5] = 5.479963503528103E1*(1.0/3.2E1);
        matrix[13][7] = 5.479963503528103E1*(-1.0/1.68E2);
        matrix[13][9] = 5.479963503528103E1*4.96031746031746E-4;
        matrix[13][11] = 5.479963503528103E1*(-1.803751803751804E-5);
        matrix[13][13] = 5.479963503528103E1*2.312502312502313E-7;
    }

    if (degree >= 14) {
        matrix[14][0] = sqrt(8.58E2)*(-1.0/6.4E1);
        matrix[14][2] = sqrt(8.58E2)*(7.0/6.4E1);
        matrix[14][4] = sqrt(8.58E2)*(-7.0/6.4E1);
        matrix[14][6] = sqrt(8.58E2)*(7.0/1.92E2);
        matrix[14][8] = sqrt(8.58E2)*(-1.0/1.92E2);
        matrix[14][10] = sqrt(8.58E2)*3.472222222222222E-4;
        matrix[14][12] = sqrt(8.58E2)*(-1.052188552188552E-5);
        matrix[14][14] = sqrt(8.58E2)*1.156251156251156E-7;
    }

    if (degree >= 15) {
        matrix[15][1] = 3.781534080237807E1*(-3.0/6.4E1);
        matrix[15][3] = 3.781534080237807E1*(7.0/6.4E1);
        matrix[15][5] = 3.781534080237807E1*(-2.1E1/3.2E2);
        matrix[15][7] = 3.781534080237807E1*(1.0/6.4E1);
        matrix[15][9] = 3.781534080237807E1*(-1.0/5.76E2);
        matrix[15][11] = 3.781534080237807E1*9.46969696969697E-5;
        matrix[15][13] = 3.781534080237807E1*(-2.428127428127428E-6);
        matrix[15][15] = 3.781534080237807E1*2.312502312502313E-8;
    }

    if (degree >= 16) {
        matrix[16][0] = 3.781534080237807E1*(3.0/2.56E2);
        matrix[16][2] = 3.781534080237807E1*(-3.0/3.2E1);
        matrix[16][4] = 3.781534080237807E1*(7.0/6.4E1);
        matrix[16][6] = 3.781534080237807E1*(-7.0/1.6E2);
        matrix[16][8] = 3.781534080237807E1*(1.0/1.28E2);
        matrix[16][10] = 3.781534080237807E1*(-6.944444444444444E-4);
        matrix[16][12] = 3.781534080237807E1*3.156565656565657E-5;
        matrix[16][14] = 3.781534080237807E1*(-6.937506937506938E-7);
        matrix[16][16] = 3.781534080237807E1*5.781255781255781E-9;
    }

    if (degree >= 17) {
        matrix[17][1] = 1.559166443969341E2*(3.0/2.56E2);
        matrix[17][3] = 1.559166443969341E2*(-1.0/3.2E1);
        matrix[17][5] = 1.559166443969341E2*(7.0/3.2E2);
        matrix[17][7] = 1.559166443969341E2*(-1.0/1.6E2);
        matrix[17][9] = 1.559166443969341E2*8.680555555555556E-4;
        matrix[17][11] = 1.559166443969341E2*(-6.313131313131313E-5);
        matrix[17][13] = 1.559166443969341E2*2.428127428127428E-6;
        matrix[17][15] = 1.559166443969341E2*(-4.625004625004625E-8);
        matrix[17][17] = 1.559166443969341E2*3.400738694856342E-10;
    }

    if (degree >= 18) {
        matrix[18][0] = 1.102497165529236E2*(-1.0/2.56E2);
        matrix[18][2] = 1.102497165529236E2*(9.0/2.56E2);
        matrix[18][4] = 1.102497165529236E2*(-3.0/6.4E1);
        matrix[18][6] = 1.102497165529236E2*(7.0/3.2E2);
        matrix[18][8] = 1.102497165529236E2*(-3.0/6.4E2);
        matrix[18][10] = 1.102497165529236E2*5.208333333333333E-4;
        matrix[18][12] = 1.102497165529236E2*(-3.156565656565657E-5);
        matrix[18][14] = 1.102497165529236E2*1.040626040626041E-6;
        matrix[18][16] = 1.102497165529236E2*(-1.734376734376734E-8);
        matrix[18][18] = 1.102497165529236E2*1.133579564952114E-10;
    }

    if (degree >= 19) {
        matrix[19][1] = 4.805673730081975E2*(-1.0/2.56E2);
        matrix[19][3] = 4.805673730081975E2*(3.0/2.56E2);
        matrix[19][5] = 4.805673730081975E2*(-3.0/3.2E2);
        matrix[19][7] = 4.805673730081975E2*(1.0/3.2E2);
        matrix[19][9] = 4.805673730081975E2*(-5.208333333333333E-4);
        matrix[19][11] = 4.805673730081975E2*4.734848484848485E-5;
        matrix[19][13] = 4.805673730081975E2*(-2.428127428127428E-6);
        matrix[19][15] = 4.805673730081975E2*6.937506937506938E-8;
        matrix[19][17] = 4.805673730081975E2*(-1.020221608456903E-9);
        matrix[19][19] = 4.805673730081975E2*5.966208236590074E-12;
    }

    if (degree >= 20) {
        matrix[20][0] = 2.149162627629654E2*(1.0/5.12E2);
        matrix[20][2] = 2.149162627629654E2*(-5.0/2.56E2);
        matrix[20][4] = 2.149162627629654E2*(1.5E1/5.12E2);
        matrix[20][6] = 2.149162627629654E2*(-1.0/6.4E1);
        matrix[20][8] = 2.149162627629654E2*(1.0/2.56E2);
        matrix[20][10] = 2.149162627629654E2*(-5.208333333333333E-4);
        matrix[20][12] = 2.149162627629654E2*3.945707070707071E-5;
        matrix[20][14] = 2.149162627629654E2*(-1.734376734376734E-6);
        matrix[20][16] = 2.149162627629654E2*4.335941835941836E-8;
        matrix[20][18] = 2.149162627629654E2*(-5.66789782476057E-10);
        matrix[20][20] = 2.149162627629654E2*2.983104118295037E-12;
    }

    if (degree >= 21) {
        matrix[21][1] = 9.848700421883082E2*(1.0/5.12E2);
        matrix[21][3] = 9.848700421883082E2*(-5.0/7.68E2);
        matrix[21][5] = 9.848700421883082E2*(3.0/5.12E2);
        matrix[21][7] = 9.848700421883082E2*(-1.0/4.48E2);
        matrix[21][9] = 9.848700421883082E2*4.340277777777778E-4;
        matrix[21][11] = 9.848700421883082E2*(-4.734848484848485E-5);
        matrix[21][13] = 9.848700421883082E2*3.035159285159285E-6;
        matrix[21][15] = 9.848700421883082E2*(-1.156251156251156E-7);
        matrix[21][17] = 9.848700421883082E2*2.550554021142256E-9;
        matrix[21][19] = 9.848700421883082E2*(-2.983104118295037E-11);
        matrix[21][21] = 9.848700421883082E2*1.420525770616684E-13;
    }

    if (degree >= 22) {
        matrix[22][0] = 4.199499970234552E2*(-9.765625E-4);
        matrix[22][2] = 4.199499970234552E2*1.07421875E-2;
        matrix[22][4] = 4.199499970234552E2*(-1.790364583333333E-2);
        matrix[22][6] = 4.199499970234552E2*1.07421875E-2;
        matrix[22][8] = 4.199499970234552E2*(-3.069196428571429E-3);
        matrix[22][10] = 4.199499970234552E2*4.774305555555556E-4;
        matrix[22][12] = 4.199499970234552E2*(-4.340277777777778E-5);
        matrix[22][14] = 4.199499970234552E2*2.38476800976801E-6;
        matrix[22][16] = 4.199499970234552E2*(-7.949226699226699E-8);
        matrix[22][18] = 4.199499970234552E2*1.558671901809157E-9;
        matrix[22][20] = 4.199499970234552E2*(-1.64070726506227E-11);
        matrix[22][22] = 4.199499970234552E2*7.102628853083421E-14;
    }

    if (degree >= 23) {
        matrix[23][1] = 2.014009433940169E3*(-9.765625E-4);
        matrix[23][3] = 2.014009433940169E3*3.580729166666667E-3;
        matrix[23][5] = 2.014009433940169E3*(-3.580729166666667E-3);
        matrix[23][7] = 2.014009433940169E3*1.534598214285714E-3;
        matrix[23][9] = 2.014009433940169E3*(-3.410218253968254E-4);
        matrix[23][11] = 2.014009433940169E3*4.340277777777778E-5;
        matrix[23][13] = 2.014009433940169E3*(-3.338675213675214E-6);
        matrix[23][15] = 2.014009433940169E3*1.58984533984534E-7;
        matrix[23][17] = 2.014009433940169E3*(-4.67601570542747E-9);
        matrix[23][19] = 2.014009433940169E3*8.203536325311351E-11;
        matrix[23][21] = 2.014009433940169E3*(-7.812891738391763E-13);
        matrix[23][23] = 2.014009433940169E3*3.088099501340618E-15;
    }

    if (degree >= 24) {
        matrix[24][0] = 8.222159083841665E2*4.8828125E-4;
        matrix[24][2] = 8.222159083841665E2*(-3.0/5.12E2);
        matrix[24][4] = 8.222159083841665E2*1.07421875E-2;
        matrix[24][6] = 8.222159083841665E2*(-7.161458333333333E-3);
        matrix[24][8] = 8.222159083841665E2*2.301897321428571E-3;
        matrix[24][10] = 8.222159083841665E2*(-4.092261904761905E-4);
        matrix[24][12] = 8.222159083841665E2*4.340277777777778E-5;
        matrix[24][14] = 8.222159083841665E2*(-2.861721611721612E-6);
        matrix[24][16] = 8.222159083841665E2*1.192384004884005E-7;
        matrix[24][18] = 8.222159083841665E2*(-3.117343803618313E-9);
        matrix[24][20] = 8.222159083841665E2*4.922121795186811E-11;
        matrix[24][22] = 8.222159083841665E2*(-4.261577311850053E-13);
        matrix[24][24] = 8.222159083841665E2*1.544049750670309E-15;
    }

    if (degree >= 25) {
        matrix[25][1] = 8.222159083841665E2*2.44140625E-3;
        matrix[25][3] = 8.222159083841665E2*(-5.0/5.12E2);
        matrix[25][5] = 8.222159083841665E2*1.07421875E-2;
        matrix[25][7] = 8.222159083841665E2*(-5.115327380952381E-3);
        matrix[25][9] = 8.222159083841665E2*1.278831845238095E-3;
        matrix[25][11] = 8.222159083841665E2*(-1.860119047619048E-4);
        matrix[25][13] = 8.222159083841665E2*1.669337606837607E-5;
        matrix[25][15] = 8.222159083841665E2*(-9.539072039072039E-7);
        matrix[25][17] = 8.222159083841665E2*3.507011779070603E-8;
        matrix[25][19] = 8.222159083841665E2*(-8.203536325311351E-10);
        matrix[25][21] = 8.222159083841665E2*1.171933760758764E-11;
        matrix[25][23] = 8.222159083841665E2*(-9.264298504021853E-14);
        matrix[25][25] = 8.222159083841665E2*3.088099501340618E-16;
    }

    if (degree >= 26) {
        matrix[26][0] = 3.224996124028679E2*(-1.220703125E-3);
        matrix[26][2] = 3.224996124028679E2*1.5869140625E-2;
        matrix[26][4] = 3.224996124028679E2*(-3.173828125E-2);
        matrix[26][6] = 3.224996124028679E2*2.327473958333333E-2;
        matrix[26][8] = 3.224996124028679E2*(-8.312406994047619E-3);
        matrix[26][10] = 3.224996124028679E2*1.662481398809524E-3;
        matrix[26][12] = 3.224996124028679E2*(-2.015128968253968E-4);
        matrix[26][14] = 3.224996124028679E2*1.550099206349206E-5;
        matrix[26][16] = 3.224996124028679E2*(-7.750496031746032E-7);
        matrix[26][18] = 3.224996124028679E2*2.53284184043988E-8;
        matrix[26][20] = 3.224996124028679E2*(-5.332298611452378E-10);
        matrix[26][22] = 3.224996124028679E2*6.925063131756335E-12;
        matrix[26][24] = 3.224996124028679E2*(-5.018161689678504E-14);
        matrix[26][26] = 3.224996124028679E2*1.544049750670309E-16;
    }

    if (degree >= 27) {
        matrix[27][1] = 5.585857141030372E2*(-3.662109375E-3);
        matrix[27][3] = 5.585857141030372E2*1.5869140625E-2;
        matrix[27][5] = 5.585857141030372E2*(-1.904296875E-2);
        matrix[27][7] = 5.585857141030372E2*9.974888392857143E-3;
        matrix[27][9] = 5.585857141030372E2*(-2.770802331349206E-3);
        matrix[27][11] = 5.585857141030372E2*4.534040178571429E-4;
        matrix[27][13] = 5.585857141030372E2*(-4.650297619047619E-5);
        matrix[27][15] = 5.585857141030372E2*3.100198412698413E-6;
        matrix[27][17] = 5.585857141030372E2*(-1.367734593837535E-7);
        matrix[27][19] = 5.585857141030372E2*3.999223958589284E-9;
        matrix[27][21] = 5.585857141030372E2*(-7.617569444931969E-11);
        matrix[27][23] = 5.585857141030372E2*9.032691041421307E-13;
        matrix[27][25] = 5.585857141030372E2*(-6.021794027614205E-15);
        matrix[27][27] = 5.585857141030372E2*1.715610834078121E-17;
    }

    if (degree >= 28) {
        matrix[28][0] = 2.111255550614373E2*1.8310546875E-3;
        matrix[28][2] = 2.111255550614373E2*(-2.5634765625E-2);
        matrix[28][4] = 2.111255550614373E2*5.55419921875E-2;
        matrix[28][6] = 2.111255550614373E2*(-4.443359375E-2);
        matrix[28][8] = 2.111255550614373E2*1.74560546875E-2;
        matrix[28][10] = 2.111255550614373E2*(-3.879123263888889E-3);
        matrix[28][12] = 2.111255550614373E2*5.289713541666667E-4;
        matrix[28][14] = 2.111255550614373E2*(-4.650297619047619E-5);
        matrix[28][16] = 2.111255550614373E2*2.712673611111111E-6;
        matrix[28][18] = 2.111255550614373E2*(-1.063793572984749E-7);
        matrix[28][20] = 2.111255550614373E2*2.799456771012499E-9;
        matrix[28][22] = 2.111255550614373E2*(-4.847544192229435E-11);
        matrix[28][24] = 2.111255550614373E2*5.269069774162429E-13;
        matrix[28][26] = 2.111255550614373E2*(-3.242504476407649E-15);
        matrix[28][28] = 2.111255550614373E2*8.578054170390605E-18;
    }

    if (degree >= 29) {
        matrix[29][1] = 1.13694590900359E3*1.8310546875E-3;
        matrix[29][3] = 1.13694590900359E3*(-8.544921875E-3);
        matrix[29][5] = 1.13694590900359E3*1.11083984375E-2;
        matrix[29][7] = 1.13694590900359E3*(-6.34765625E-3);
        matrix[29][9] = 1.13694590900359E3*1.939561631944444E-3;
        matrix[29][11] = 1.13694590900359E3*(-3.526475694444444E-4);
        matrix[29][13] = 1.13694590900359E3*4.069010416666667E-5;
        matrix[29][15] = 1.13694590900359E3*(-3.100198412698413E-6);
        matrix[29][17] = 1.13694590900359E3*1.595690359477124E-7;
        matrix[29][19] = 1.13694590900359E3*(-5.598913542024997E-9);
        matrix[29][21] = 1.13694590900359E3*1.333074652863095E-10;
        matrix[29][23] = 1.13694590900359E3*(-2.107627909664972E-12);
        matrix[29][25] = 1.13694590900359E3*2.107627909664972E-14;
        matrix[29][27] = 1.13694590900359E3*(-1.200927583854685E-16);
        matrix[29][29] = 1.13694590900359E3*2.957949713927795E-19;
    }

    if (degree >= 30) {
        matrix[30][0] = 1.037884868374137E3*(-3.662109375E-4);
        matrix[30][2] = 1.037884868374137E3*5.4931640625E-3;
        matrix[30][4] = 1.037884868374137E3*(-1.28173828125E-2);
        matrix[30][6] = 1.037884868374137E3*1.11083984375E-2;
        matrix[30][8] = 1.037884868374137E3*(-4.7607421875E-3);
        matrix[30][10] = 1.037884868374137E3*1.163736979166667E-3;
        matrix[30][12] = 1.037884868374137E3*(-1.763237847222222E-4);
        matrix[30][14] = 1.037884868374137E3*1.743861607142857E-5;
        matrix[30][16] = 1.037884868374137E3*(-1.162574404761905E-6);
        matrix[30][18] = 1.037884868374137E3*5.318967864923747E-8;
        matrix[30][20] = 1.037884868374137E3*(-1.679674062607499E-9);
        matrix[30][22] = 1.037884868374137E3*3.635658144172076E-11;
        matrix[30][24] = 1.037884868374137E3*(-5.269069774162429E-13);
        matrix[30][26] = 1.037884868374137E3*4.863756714611473E-15;
        matrix[30][28] = 1.037884868374137E3*(-2.573416251117181E-17);
        matrix[30][30] = 1.037884868374137E3*5.91589942785559E-20;
    }


    if (degree >= 31) {
        matrix[31][1] = 5.778698382854049E3*(-3.662109375E-4);
        matrix[31][3] = 5.778698382854049E3*1.8310546875E-3;
        matrix[31][5] = 5.778698382854049E3*(-2.5634765625E-3);
        matrix[31][7] = 5.778698382854049E3*1.5869140625E-3;
        matrix[31][9] = 5.778698382854049E3*(-5.289713541666667E-4);
        matrix[31][11] = 5.778698382854049E3*1.057942708333333E-4;
        matrix[31][13] = 5.778698382854049E3*(-1.356336805555556E-5);
        matrix[31][15] = 5.778698382854049E3*1.162574404761905E-6;
        matrix[31][17] = 5.778698382854049E3*(-6.838672969187675E-8);
        matrix[31][19] = 5.778698382854049E3*2.799456771012499E-9;
        matrix[31][21] = 5.778698382854049E3*(-7.998447917178567E-11);
        matrix[31][23] = 5.778698382854049E3*1.580720932248729E-12;
        matrix[31][25] = 5.778698382854049E3*(-2.107627909664972E-14);
        matrix[31][27] = 5.778698382854049E3*1.801391375782027E-16;
        matrix[31][29] = 5.778698382854049E3*(-8.873849141783384E-19);
        matrix[31][31] = 5.778698382854049E3*1.908354654146964E-21;
    }

    if (degree >= 32) {
        matrix[32][0] = 8.172313625895668E3*4.57763671875E-5;
        matrix[32][2] = 8.172313625895668E3*(-7.32421875E-4);
        matrix[32][4] = 8.172313625895668E3*1.8310546875E-3;
        matrix[32][6] = 8.172313625895668E3*(-1.708984375E-3);
        matrix[32][8] = 8.172313625895668E3*7.9345703125E-4;
        matrix[32][10] = 8.172313625895668E3*(-2.115885416666667E-4);
        matrix[32][12] = 8.172313625895668E3*3.526475694444444E-5;
        matrix[32][14] = 8.172313625895668E3*(-3.875248015873016E-6);
        matrix[32][16] = 8.172313625895668E3*2.906436011904762E-7;
        matrix[32][18] = 8.172313625895668E3*(-1.519705104263928E-8);
        matrix[32][20] = 8.172313625895668E3*5.598913542024997E-10;
        matrix[32][22] = 8.172313625895668E3*(-1.45426325766883E-11);
        matrix[32][24] = 8.172313625895668E3*2.634534887081215E-13;
        matrix[32][26] = 8.172313625895668E3*(-3.242504476407649E-15);
        matrix[32][28] = 8.172313625895668E3*2.573416251117181E-17;
        matrix[32][30] = 8.172313625895668E3*(-1.183179885571118E-19);
        matrix[32][32] = 8.172313625895668E3*2.385443317683705E-22;
    }

    if (degree >= 33) {
        matrix[33][1] = 4.694636759111401E4*4.57763671875E-5;
        matrix[33][3] = 4.694636759111401E4*(-2.44140625E-4);
        matrix[33][5] = 4.694636759111401E4*3.662109375E-4;
        matrix[33][7] = 4.694636759111401E4*(-2.44140625E-4);
        matrix[33][9] = 4.694636759111401E4*8.816189236111111E-5;
        matrix[33][11] = 4.694636759111401E4*(-1.923532196969697E-5);
        matrix[33][13] = 4.694636759111401E4*2.712673611111111E-6;
        matrix[33][15] = 4.694636759111401E4*(-2.583498677248677E-7);
        matrix[33][17] = 4.694636759111401E4*1.709668242296919E-8;
        matrix[33][19] = 4.694636759111401E4*(-7.998447917178567E-10);
        matrix[33][21] = 4.694636759111401E4*2.666149305726189E-11;
        matrix[33][23] = 4.694636759111401E4*(-6.322883728994915E-13);
        matrix[33][25] = 4.694636759111401E4*1.053813954832486E-14;
        matrix[33][27] = 4.694636759111401E4*(-1.200927583854685E-16);
        matrix[33][29] = 4.694636759111401E4*8.873849141783384E-19;
        matrix[33][31] = 4.694636759111401E4*(-3.816709308293929E-21);
        matrix[33][33] = 4.694636759111401E4*7.228616114193047E-24;
    }

    if (degree >= 34) {
        matrix[34][0] = 8.051235619456184E3*(-4.57763671875E-5);
        matrix[34][2] = 8.051235619456184E3*7.781982421875E-4;
        matrix[34][4] = 8.051235619456184E3*(-2.0751953125E-3);
        matrix[34][6] = 8.051235619456184E3*2.0751953125E-3;
        matrix[34][8] = 8.051235619456184E3*(-1.03759765625E-3);
        matrix[34][10] = 8.051235619456184E3*2.997504340277778E-4;
        matrix[34][12] = 8.051235619456184E3*(-5.450007891414141E-5);
        matrix[34][14] = 8.051235619456184E3*6.587921626984127E-6;
        matrix[34][16] = 8.051235619456184E3*(-5.489934689153439E-7);
        matrix[34][18] = 8.051235619456184E3*3.229373346560847E-8;
        matrix[34][20] = 8.051235619456184E3*(-1.359736145920356E-9);
        matrix[34][22] = 8.051235619456184E3*4.12041256339502E-11;
        matrix[34][24] = 8.051235619456184E3*(-8.957418616076129E-13);
        matrix[34][26] = 8.051235619456184E3*1.378064402473251E-14;
        matrix[34][28] = 8.051235619456184E3*(-1.458269208966403E-16);
        matrix[34][30] = 8.051235619456184E3*1.00570290273545E-18;
        matrix[34][32] = 8.051235619456184E3*(-4.055253640062299E-21);
        matrix[34][34] = 8.051235619456184E3*7.228616114193047E-24;
    }

    if (degree >= 35) {
        matrix[35][1] = 9.526350455447249E3*(-2.288818359375E-4);
        matrix[35][3] = 9.526350455447249E3*1.2969970703125E-3;
        matrix[35][5] = 9.526350455447249E3*(-2.0751953125E-3);
        matrix[35][7] = 9.526350455447249E3*1.482282366071429E-3;
        matrix[35][9] = 9.526350455447249E3*(-5.764431423611111E-4);
        matrix[35][11] = 9.526350455447249E3*1.362501972853535E-4;
        matrix[35][13] = 9.526350455447249E3*(-2.096156881313131E-5);
        matrix[35][15] = 9.526350455447249E3*2.195973875661376E-6;
        matrix[35][17] = 9.526350455447249E3*(-1.614686673280423E-7);
        matrix[35][19] = 9.526350455447249E3*8.498350912002228E-9;
        matrix[35][21] = 9.526350455447249E3*(-3.237467014096087E-10);
        matrix[35][23] = 9.526350455447249E3*8.957418616076129E-12;
        matrix[35][25] = 9.526350455447249E3*(-1.791483723215226E-13);
        matrix[35][27] = 9.526350455447249E3*2.551971115691205E-15;
        matrix[35][29] = 9.526350455447249E3*(-2.514257256838626E-17);
        matrix[35][31] = 9.526350455447249E3*1.62210145602492E-19;
        matrix[35][33] = 9.526350455447249E3*(-6.14432369706409E-22);
        matrix[35][35] = 9.526350455447249E3*1.032659444884721E-24;
    }

    if (degree >= 36) {
        matrix[36][0] = 9.526350455447249E3*3.814697265625E-5;
        matrix[36][2] = 9.526350455447249E3*(-6.866455078125E-4);
        matrix[36][4] = 9.526350455447249E3*1.94549560546875E-3;
        matrix[36][6] = 9.526350455447249E3*(-2.0751953125E-3);
        matrix[36][8] = 9.526350455447249E3*1.111711774553571E-3;
        matrix[36][10] = 9.526350455447249E3*(-3.458658854166667E-4);
        matrix[36][12] = 9.526350455447249E3*6.812509864267677E-5;
        matrix[36][14] = 9.526350455447249E3*(-8.983529491341991E-6);
        matrix[36][16] = 9.526350455447249E3*8.234902033730159E-7;
        matrix[36][18] = 9.526350455447249E3*(-5.382288910934744E-8);
        matrix[36][20] = 9.526350455447249E3*2.549505273600668E-9;
        matrix[36][22] = 9.526350455447249E3*(-8.829455492989328E-11);
        matrix[36][24] = 9.526350455447249E3*2.239354654019032E-12;
        matrix[36][26] = 9.526350455447249E3*(-4.134193207419752E-14);
        matrix[36][28] = 9.526350455447249E3*5.468509533624011E-16;
        matrix[36][30] = 9.526350455447249E3*(-5.028514513677251E-18);
        matrix[36][32] = 9.526350455447249E3*3.041440230046724E-20;
        matrix[36][34] = 9.526350455447249E3*(-1.084292417128957E-22);
        matrix[36][36] = 9.526350455447249E3*1.721099074807868E-25;
    }

    if (degree >= 37) {
        matrix[37][1] = 5.79465276008839E4*3.814697265625E-5;
        matrix[37][3] = 5.79465276008839E4*(-2.288818359375E-4);
        matrix[37][5] = 5.79465276008839E4*3.8909912109375E-4;
        matrix[37][7] = 5.79465276008839E4*(-2.964564732142857E-4);
        matrix[37][9] = 5.79465276008839E4*1.235235305059524E-4;
        matrix[37][11] = 5.79465276008839E4*(-3.144235321969697E-5);
        matrix[37][13] = 5.79465276008839E4*5.240392203282828E-6;
        matrix[37][15] = 5.79465276008839E4*(-5.989019660894661E-7);
        matrix[37][17] = 5.79465276008839E4*4.84406001984127E-8;
        matrix[37][19] = 5.79465276008839E4*(-2.832783637334076E-9);
        matrix[37][21] = 5.79465276008839E4*1.214050130286033E-10;
        matrix[37][23] = 5.79465276008839E4*(-3.838893692604055E-12);
        matrix[37][25] = 5.79465276008839E4*8.957418616076129E-14;
        matrix[37][27] = 5.79465276008839E4*(-1.531182669414723E-15);
        matrix[37][29] = 5.79465276008839E4*1.885692942628969E-17;
        matrix[37][31] = 5.79465276008839E4*(-1.62210145602492E-19);
        matrix[37][33] = 5.79465276008839E4*9.216485545596135E-22;
        matrix[37][35] = 5.79465276008839E4*(-3.097978334654163E-24);
        matrix[37][37] = 5.79465276008839E4*4.651619121102347E-27;
    }

    if (degree >= 38) {
        matrix[38][0] = 1.880033611401669E4*(-1.9073486328125E-5);
        matrix[38][2] = 1.880033611401669E4*3.62396240234375E-4;
        matrix[38][4] = 1.880033611401669E4*(-1.087188720703125E-3);
        matrix[38][6] = 1.880033611401669E4*1.232147216796875E-3;
        matrix[38][8] = 1.880033611401669E4*(-7.040841238839286E-4);
        matrix[38][10] = 1.880033611401669E4*2.346947079613095E-4;
        matrix[38][12] = 1.880033611401669E4*(-4.978372593118687E-5);
        matrix[38][14] = 1.880033611401669E4*7.11196084731241E-6;
        matrix[38][16] = 1.880033611401669E4*(-7.11196084731241E-7);
        matrix[38][18] = 1.880033611401669E4*5.113174465388007E-8;
        matrix[38][20] = 1.880033611401669E4*(-2.691144455467372E-9);
        matrix[38][22] = 1.880033611401669E4*1.048497839792483E-10;
        matrix[38][24] = 1.880033611401669E4*(-3.039124173311544E-12);
        matrix[38][26] = 1.880033611401669E4*6.545805911747941E-14;
        matrix[38][28] = 1.880033611401669E4*(-1.039016811388562E-15);
        matrix[38][30] = 1.880033611401669E4*1.194272196998347E-17;
        matrix[38][32] = 1.880033611401669E4*(-9.631227395147961E-20);
        matrix[38][34] = 1.880033611401669E4*5.150388981362546E-22;
        matrix[38][36] = 1.880033611401669E4*(-1.635044121067475E-24);
        matrix[38][38] = 1.880033611401669E4*2.325809560551173E-27;
    }

    if (degree >= 39) {
        matrix[39][1] = 3.913602046708377E4*(-5.7220458984375E-5);
        matrix[39][3] = 3.913602046708377E4*3.62396240234375E-4;
        matrix[39][5] = 3.913602046708377E4*(-6.52313232421875E-4);
        matrix[39][7] = 3.913602046708377E4*5.280630929129464E-4;
        matrix[39][9] = 3.913602046708377E4*(-2.346947079613095E-4);
        matrix[39][11] = 3.913602046708377E4*6.400764762581169E-5;
        matrix[39][13] = 3.913602046708377E4*(-1.14885521379662E-5);
        matrix[39][15] = 3.913602046708377E4*1.422392169462482E-6;
        matrix[39][17] = 3.913602046708377E4*(-1.255051914231602E-7);
        matrix[39][19] = 3.913602046708377E4*8.073433366402116E-9;
        matrix[39][21] = 3.913602046708377E4*(-3.844492079239103E-10);
        matrix[39][23] = 3.913602046708377E4*1.367605877990195E-11;
        matrix[39][25] = 3.913602046708377E4*(-3.646949007973853E-13);
        matrix[39][27] = 3.913602046708377E4*7.273117679719934E-15;
        matrix[39][29] = 3.913602046708377E4*(-1.074844977298512E-16);
        matrix[39][31] = 3.913602046708377E4*1.155747287417755E-18;
        matrix[39][33] = 3.913602046708377E4*(-8.755661268316328E-21);
        matrix[39][35] = 3.913602046708377E4*4.414619126882182E-23;
        matrix[39][37] = 3.913602046708377E4*(-1.325711449514169E-25);
        matrix[39][39] = 3.913602046708377E4*1.789084277347056E-28;
    }

    if (degree >= 40) {
        matrix[40][0] = 6.187948161547574E4*5.7220458984375E-6;
        matrix[40][2] = 6.187948161547574E4*(-1.1444091796875E-4);
        matrix[40][4] = 6.187948161547574E4*3.62396240234375E-4;
        matrix[40][6] = 6.187948161547574E4*(-4.3487548828125E-4);
        matrix[40][8] = 6.187948161547574E4*2.640315464564732E-4;
        matrix[40][10] = 6.187948161547574E4*(-9.387788318452381E-5);
        matrix[40][12] = 6.187948161547574E4*2.133588254193723E-5;
        matrix[40][14] = 6.187948161547574E4*(-3.282443467990343E-6);
        matrix[40][16] = 6.187948161547574E4*3.555980423656205E-7;
        matrix[40][18] = 6.187948161547574E4*(-2.789004253848004E-8);
        matrix[40][20] = 6.187948161547574E4*1.614686673280423E-9;
        matrix[40][22] = 6.187948161547574E4*(-6.989985598616551E-11);
        matrix[40][24] = 6.187948161547574E4*2.279343129983658E-12;
        matrix[40][26] = 6.187948161547574E4*(-5.610690781498235E-14);
        matrix[40][28] = 6.187948161547574E4*1.039016811388562E-15;
        matrix[40][30] = 6.187948161547574E4*(-1.433126636398017E-17);
        matrix[40][32] = 6.187948161547574E4*1.444684109272194E-19;
        matrix[40][34] = 6.187948161547574E4*(-1.030077796272509E-21);
        matrix[40][36] = 6.187948161547574E4*4.905132363202425E-24;
        matrix[40][38] = 6.187948161547574E4*(-1.395485736330704E-26);
        matrix[40][40] = 6.187948161547574E4*1.789084277347056E-29;
    }

    if (degree > 40) {
        cout << "Degree too high" << endl; exit(0);
    }
}
