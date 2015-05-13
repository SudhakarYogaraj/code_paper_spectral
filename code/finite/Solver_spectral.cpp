/* TODO: Linear system solution (urbain, Sat 09 May 2015 14:32:40 BST) */
/* TODO: Integration over infinite domain: which variance for gaussian(urbain, Sat 09 May 2015 14:33:06 BST) */
/* TODO: Convergence of hermite functions? (urbain, Sat 09 May 2015 14:33:28 BST) */

#include "structures.hpp"
#include "aux.hpp"
#include "Problem.hpp"
#include "Gaussian_integrator.hpp"
#include "Solver_spectral.hpp"
#include "templates.hpp"
#include <iomanip>

using namespace std;

double Solver_spectral::basis(vector<int> mult, vector<double> x, vector<double> sigmas) {

    double result;
    int monomial = 0;

    if (monomial) {
        result = 1.;
        for (unsigned int i = 0; i < mult.size(); ++i) {
            result *= ipow(x[i]/sigmas[i], mult[i]);
        }
    }

    if (!monomial) {
        result = 1.;
        for (unsigned int i = 0; i < mult.size(); ++i) {
            double result_i = 0.;
            double xi = x[i]/sigmas[i];
            for (unsigned int j = 0; j < hermiteCoeffs_1d.size(); ++j) {
                result_i += hermiteCoeffs_1d[mult[i]][j] * ipow(xi, j);
            }
            result *= result_i;
        }
    }
    return result;
}

vector<double> aux( vector< vector<double> > mat, vector<double> coefficients, vector<double> weights) {

    /* FIXME: bad conditioning (urbain, Tue 12 May 2015 18:53:50 BST) */

    int nb = mat.size();

    /* for (int i = 0; i < nb; ++i) { */
    /*     cout << weights[i] << endl; */
    /* } */

    vector< vector<double> > system_mat(nb, vector<double>(nb, 0.));
    for (int i = 1; i < nb; ++i) {
        for (int j = 1; j < nb; ++j) {
            system_mat[i][j] = mat[i][j] - weights[i]/weights[0]*mat[i][0] - weights[j]/weights[0]*mat[j][0] + weights[i]*weights[j]/(weights[0]*weights[0])*mat[0][0];
        }
    }
    system_mat[0][0] = 1.;


    /*     // print */
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            cout << setw(10) <<  mat[i][j] << " " ;
        }
        cout << endl;
    }
    cout << "*** " << endl;

    /* // print */
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            cout << setw(10) <<  system_mat[i][j] << " " ;
        }
        cout << endl;
    }


    vector<double> system_rhs(nb, 0.);
    for (int i = 0; i < nb; ++i) {
        system_rhs[i] = coefficients[i] - weights[i]/weights[0] * coefficients[0];
    }

    /* for (int i = 0; i < nb; ++i) { */
    /*     cout << coefficients[i] << endl; */
    /* } */
    /* cout << " --- " << endl; */

    /* for (int i = 0; i < nb; ++i) { */
    /*     cout << system_rhs[i] << endl; */
    /* } */

    vector<double> result = solve(system_mat, system_rhs);
    for (int i = 1; i < nb; ++i) {
        result[0] += weights[i]/weights[0] * result[i];
    }
    return result;
}

void Solver_spectral::estimator(Problem &problem, vector<double> x,  vector<SDE_coeffs>& c, double t) {

    int nf     = problem.nf;
    int ns     = problem.d;
    int nb     = bin(degree + nf, nf);

    c = vector<SDE_coeffs> (nb);
    Gaussian_integrator gauss = Gaussian_integrator(nNodes,nf);

    // Eigenvalues
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        sigmas[i] = problem.betas[i]*sqrt(0.5/problem.lambdas[i]);
    }

    // Variance parameter for hermite functions
    /* vector<double> sigmas_hf(1, 1./sqrt(2.)); */
    vector<double> sigmas_hf(1, 1./sqrt(2.));

    // Variance parameter for integration
    vector<double> sigmas_quad(1, 1./sqrt(2.));

    // Expansion of right-hand side of the Poisson equation
    vector< vector<double> > coefficients(ns, vector<double>(nb, 0.));
    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));
    vector< vector<double> > coefficients_h(nf, vector<double>(nb,0.));
    for (int i = 0; i < nb; ++i) {

        // Graphical progression bar.
        /* cout << "["; */
        /* int bw = 103; */
        /* double progress = ( (double) i) / ( (double) nb); */
        /* int pos = bw * progress; */
        /* for (int j = 0; j < bw; ++j) { */
        /*     if (j < pos) cout << "="; */
        /*     else if (j == pos) cout << ">"; */
        /*     else cout << " "; */
        /* } */
        /* cout << "] " << int(progress * 100.0) << " %\r"; */

        vector<int> multIndex = ind2mult(i, degree, nf);

        vector<double> v0(ns,0.);
        vector<double> h0(nf,0.);
        vector< vector<double> > m0(ns, vector<double> (ns,0.));

        auto lambda = [&] (vector<double> y) -> vector<double> {
            return problem.a(x,y)*basis(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };
        auto lambda_dx = [&] (vector<double> y) -> vector< vector<double> > {
            return problem.dax(x,y)*basis(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };
        auto lambda_h = [&] (vector<double> y) -> vector<double> {
            return problem.fast_drift_h(x,y)*basis(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };

        vector<double> result = gauss.flatquadnd(lambda, sigmas_quad, v0);
        vector< vector<double> > result_dx = gauss.flatquadnd(lambda_dx, sigmas_quad, m0);
        vector<double> result_h = gauss.flatquadnd(lambda_h, sigmas_quad, h0);

        for (int j = 0; j < ns; ++j) {
            coefficients[j][i] = result[j];
            for (int k = 0; k < ns; ++k) {
                coefficients_dx[j][k][i] = result_dx[j][k];
            }
        }
        for (int j = 0; j < nf; ++j) {
            coefficients_h[j][i] = result_h[j];
        }
        cout.flush();
    }

    /* for (int i = 0; i < ns; ++i) { */
    /*     for (int j = 0; j < ns; ++j) { */
    /*         coefficients_dx[i][j] = mon2herm(coefficients_dx[i][j],nf,degree); */
    /*     } */
    /*     coefficients[i] = mon2herm(coefficients[i],nf,degree); */
    /* } */
    /* for (int i = 0; i < nf; ++i) { */
    /*     coefficients_h[i] = mon2herm(coefficients_h[i],nf,degree); */
    /* } */

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nb,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nb,0.)));


    // weights(i,j) = int (phi_i, e^(-V) )
    vector<double> weights(nb, 0.);
    for (int i = 0; i < nb; ++i) {
        vector<int> m = ind2mult(i, degree, nf);
        auto lambda = [&] (vector<double> y) -> double {
            return basis(m, y, sigmas_hf) * sqrt(problem.rho(x,y) * gaussian(y,sigmas_hf)); };
        weights[i] = gauss.flatquadnd(lambda, sigmas_hf);
    }

    // mat(i,j) = int ( L phi_i, phi_j)
    vector< vector<double> > mat(nb, vector<double>(nb, 0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i, degree, nf);
        for (int j = 0; j < nf; ++j) {
            mat[i][i] += m1[j]/(sigmas_hf[j]*sigmas_hf[j]);
        }
        for (int j = 0; j < nb; ++j) {
            vector<int> m2 = ind2mult(j, degree, nf);
            auto lambda = [&] (vector<double> y) -> double {
                double tmp = 1/(2*pow(sigmas_hf[0], 2)) - pow(y[0], 2)/(4*pow(sigmas_hf[0], 4));
                return (problem.linearTerm(x,y) - tmp) * basis(m1, y, sigmas_hf) * basis(m2, y, sigmas_hf) * gaussian(y,sigmas_hf); };
            mat[i][j] += gauss.flatquadnd(lambda, sigmas_hf);
        }
    }

    // print
    /* for (int i = 0; i < nb; ++i) { */
    /*     for (int j = 0; j < nb; ++j) { */
    /*         cout << setw(10) <<  mat[i][j] << " " ; */
    /*     } */
    /*     cout << endl; */
    /* } */

    for (int i = 0; i < ns; ++i) {
        /* vector<double> tmp = aux(mat, coefficients[i], weights); */
        /* mat[0][0] = 1.0; solution[i] = solve(mat, coefficients[i]); */
        /* cout << " ***" << endl; */
        /* for (int j = 0; j < nb; ++j) { */
        /*     cout << solution[0][j] - tmp[j] << endl; */
        /* } */
        /* cout << " ***" << endl; */
        /* for (int j = 0; j < nb; ++j) { */
        /*     cout << solution[0][j] << endl; */
        /* } */
        /* exit(0); */
        solution[i] = aux(mat, coefficients[i], weights);
        for (int j = 0; j < ns; ++j) {
            solution_dx[i][j] = aux(mat, coefficients_dx[i][j], weights);
        }
    }

    // y-Derivatives of the solution
    for (int i = 0; i < nb; ++i) {
        vector<int> mult = ind2mult(i,degree,nf);
        if (accumulate (mult.begin(), mult.end(),0) == degree) {
            continue;
        }
        for (int j = 0; j < nf; ++j) {
            vector<int> newMult = mult;
            newMult[j] ++;
            int newInd = mult2ind(newMult, degree);
            for (int k = 0; k < ns; ++k) {
                solution_dy[k][j][i] = solution[k][newInd]*sqrt(newMult[j])/sigmas[j];
            }
        }
    }

    // Calculation of the coefficients of the simplified equation
    vector<double> F1(ns, 0.);
    vector<double> F2(ns, 0.);
    vector< vector <double> > A0(ns, vector<double>(ns,0.));
    for (int j = 0; j < nb; ++j) {
        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l)
                F1[k] += solution_dx[k][l][j]*coefficients[l][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < nf; ++l)
                F2[k] += solution_dy[k][l][j]*coefficients_h[l][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l) {
                A0[k][l] += 2*solution[k][j]*coefficients[l][j];
            }
        }
        c[j].diff =  cholesky(symmetric(A0));
        c[j].drif = F1 + F2;
    }
}

Solver_spectral::Solver_spectral(int degree, int nNodes, int n_vars)
{
    this->degree = degree;
    this->nNodes = nNodes;

    vector< vector<double> > mat1d(degree + 1, vector<double>(degree + 1,0.));

    if (degree >= 0) {
        mat1d[0][0] = 1.0;
    }

    if (degree >= 1) {
        mat1d[1][1] = 1.0;
    }

    if (degree >= 2) {
        mat1d[2][0] = sqrt(2.0)/2.0*(-1.0);
        mat1d[2][2] = sqrt(2.0)/2.0*( 1.0);
    }

    if (degree >= 3) {
        mat1d[3][1] = sqrt(6.0)/6.0*(-3.0);
        mat1d[3][3] = sqrt(6.0)/6.0*( 1.0);
    }

    if (degree >= 4) {
        mat1d[4][0] = sqrt(6.0)/12.0*( 3.0);
        mat1d[4][2] = sqrt(6.0)/12.0*(-6.0);
        mat1d[4][4] = sqrt(6.0)/12.0*( 1.0);
    }

    if (degree >= 5) {
        mat1d[5][1] = sqrt(30.0)/60.0*( 15.0);
        mat1d[5][3] = sqrt(30.0)/60.0*(-10.0);
        mat1d[5][5] = sqrt(30.0)/60.0*(  1.0);
    }

    if (degree >= 6) {
        mat1d[6][0] = sqrt(5.0)/60.0*(-15.0);
        mat1d[6][2] = sqrt(5.0)/60.0*( 45.0);
        mat1d[6][4] = sqrt(5.0)/60.0*(-15.0);
        mat1d[6][6] = sqrt(5.0)/60.0*(  1.0);
    }

    if (degree >= 7) {
        mat1d[7][1] = sqrt(35.0)/420.0*(-105.0);
        mat1d[7][3] = sqrt(35.0)/420.0*( 105.0);
        mat1d[7][5] = sqrt(35.0)/420.0*( -21.0);
        mat1d[7][7] = sqrt(35.0)/420.0*(   1.0);
    }

    if (degree >= 8) {
        mat1d[8][0] = sqrt(70.0)/1680.0*( 105.0);
        mat1d[8][2] = sqrt(70.0)/1680.0*(-420.0);
        mat1d[8][4] = sqrt(70.0)/1680.0*( 210.0);
        mat1d[8][6] = sqrt(70.0)/1680.0*( -28.0);
        mat1d[8][8] = sqrt(70.0)/1680.0*(   1.0);
    }

    if (degree >= 9) {
        mat1d[9][1] = sqrt(70.0)/5040.0*(  945.0);
        mat1d[9][3] = sqrt(70.0)/5040.0*(-1260.0);
        mat1d[9][5] = sqrt(70.0)/5040.0*(  378.0);
        mat1d[9][7] = sqrt(70.0)/5040.0*(  -36.0);
        mat1d[9][9] = sqrt(70.0)/5040.0*(    1.0);
    }

    if (degree >= 10) {
        mat1d[10][0]  = sqrt(7.0)/5040.0*( -945.0);
        mat1d[10][2]  = sqrt(7.0)/5040.0*( 4725.0);
        mat1d[10][4]  = sqrt(7.0)/5040.0*(-3150.0);
        mat1d[10][6]  = sqrt(7.0)/5040.0*(  630.0);
        mat1d[10][8]  = sqrt(7.0)/5040.0*(  -45.0);
        mat1d[10][10] = sqrt(7.0)/5040.0*(    1.0);
    }

    if (degree >= 11) {
        mat1d[11][11] = (sqrt(77.))/55440.;
        mat1d[11][9]  = - (sqrt(77.))/1008.;
        mat1d[11][7]  =  (sqrt(77.))/56.;
        mat1d[11][5]  = - (sqrt(77.))/8.;
        mat1d[11][3]  =  (5.*sqrt(77.))/16.;
        mat1d[11][1]  = - (3.*sqrt(77.))/16.;
    }

    if (degree >= 12) {
        mat1d[12][12] = (sqrt(231.))/332640.;
        mat1d[12][10] = - (sqrt(231.))/5040.;
        mat1d[12][8]  =  (sqrt(231.))/224.;
        mat1d[12][6]  = - (sqrt(231.))/24.;
        mat1d[12][4]  =  (5.*sqrt(231.))/32.;
        mat1d[12][2]  = - (3.*sqrt(231.))/16.;
        mat1d[12][0]  =  sqrt(231.)/32.;
    }

    if (degree >= 13) {
        mat1d[13][13] = (sqrt(3003.))/4324320.;
        mat1d[13][11] = - (sqrt(3003.))/55440.;
        mat1d[13][9]  =  (sqrt(3003.))/2016.;
        mat1d[13][7]  = - (sqrt(3003.))/168.;
        mat1d[13][5]  =  (sqrt(3003.))/32.;
        mat1d[13][3]  = - (sqrt(3003.))/16.;
        mat1d[13][1]  =  (sqrt(3003.))/32.;
    }

    if (degree >= 14) {
        mat1d[14][14] = (sqrt(858.))/8648640.;
        mat1d[14][12] = - (sqrt(858.))/95040.;
        mat1d[14][10] =  (sqrt(858.))/2880.;
        mat1d[14][8]  = - (sqrt(858.))/192.;
        mat1d[14][6]  =  (7.*sqrt(858.))/192.;
        mat1d[14][4]  = - (7.*sqrt(858.))/64.;
        mat1d[14][2]  =  (7.*sqrt(858.))/64.;
        mat1d[14][0]  = - sqrt(858.)/64.;
    }

    if (degree >= 15) {
        mat1d[15][15] = (sqrt(1430.))/43243200.;
        mat1d[15][13] = - (sqrt(1430.))/411840.;
        mat1d[15][11] =  (sqrt(1430.))/10560.;
        mat1d[15][9]  = - (sqrt(1430.))/576.;
        mat1d[15][7]  =  (sqrt(1430.))/64.;
        mat1d[15][5]  = - (21.*sqrt(1430.))/320.;
        mat1d[15][3]  =  (7.*sqrt(1430.))/64.;
        mat1d[15][1]  = - (3.*sqrt(1430.))/64.;
    }

    if (degree >= 16) {
        mat1d[16][16] = (sqrt(1430.))/172972800.;
        mat1d[16][14] = - (sqrt(1430.))/1441440.;
        mat1d[16][12] =  (sqrt(1430.))/31680.;
        mat1d[16][10] = - (sqrt(1430.))/1440.;
        mat1d[16][8]  =  (sqrt(1430.))/128.;
        mat1d[16][6]  = - (7.*sqrt(1430.))/160.;
        mat1d[16][4]  =  (7.*sqrt(1430.))/64.;
        mat1d[16][2]  = - (3.*sqrt(1430.))/32.;
        mat1d[16][0]  =  (3.*sqrt(1430.))/256.;
    }

    if (degree >= 17) {
        mat1d[17][17] = (sqrt(24310.))/2940537600.;
        mat1d[17][15] = - (sqrt(24310.))/21621600.;
        mat1d[17][13] =  (sqrt(24310.))/411840.;
        mat1d[17][11] = - (sqrt(24310.))/15840.;
        mat1d[17][9]  =  (sqrt(24310.))/1152.;
        mat1d[17][7]  = - (sqrt(24310.))/160.;
        mat1d[17][5]  =  (7.*sqrt(24310.))/320.;
        mat1d[17][3]  = - (sqrt(24310.))/32.;
        mat1d[17][1]  =  (3.*sqrt(24310.))/256.;
    }

    if (degree >= 18) {
        mat1d[18][18] = (sqrt(12155.))/8821612800.;
        mat1d[18][16] = - (sqrt(12155.))/57657600.;
        mat1d[18][14] =  (sqrt(12155.))/960960.;
        mat1d[18][12] = - (sqrt(12155.))/31680.;
        mat1d[18][10] =  (sqrt(12155.))/1920.;
        mat1d[18][8]  = - (3.*sqrt(12155.))/640.;
        mat1d[18][6]  =  (7.*sqrt(12155.))/320.;
        mat1d[18][4]  = - (3.*sqrt(12155.))/64.;
        mat1d[18][2]  =  (9.*sqrt(12155.))/256.;
        mat1d[18][0]  = - sqrt(12155.)/256.;
    }

    if (degree >= 19) {
        mat1d[19][19] = (sqrt(230945.))/167610643200.;
        mat1d[19][17] = - (sqrt(230945.))/980179200.;
        mat1d[19][15] =  (sqrt(230945.))/14414400.;
        mat1d[19][13] = - (sqrt(230945.))/411840.;
        mat1d[19][11] =  (sqrt(230945.))/21120.;
        mat1d[19][9]  = - (sqrt(230945.))/1920.;
        mat1d[19][7]  =  (sqrt(230945.))/320.;
        mat1d[19][5]  = - (3.*sqrt(230945.))/320.;
        mat1d[19][3]  =  (3.*sqrt(230945.))/256.;
        mat1d[19][1]  = - (sqrt(230945.))/256.;
    }

    if (degree >= 20) {
        mat1d[20][20] = (sqrt(46189.))/335221286400.;
        mat1d[20][18] = - (sqrt(46189.))/1764322560.;
        mat1d[20][16] =  (sqrt(46189.))/23063040.;
        mat1d[20][14] = - (sqrt(46189.))/576576.;
        mat1d[20][12] =  (sqrt(46189.))/25344.;
        mat1d[20][10] = - (sqrt(46189.))/1920.;
        mat1d[20][8]  =  (sqrt(46189.))/256.;
        mat1d[20][6]  = - (sqrt(46189.))/64.;
        mat1d[20][4]  =  (15.*sqrt(46189.))/512.;
        mat1d[20][2]  = - (5.*sqrt(46189.))/256.;
        mat1d[20][0]  =  sqrt(46189.)/512.;
    }

    if (degree > 20) {
        cout << "Degree too high" << endl; exit(0);
    }

    int nb = bin(degree + n_vars, n_vars);
    vector< vector<double> > matnd(nb, vector<double>(nb,0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i,degree,n_vars);
        for (int j = 0; j < nb; ++j) {
            vector<int> m2 = ind2mult(j,degree,n_vars);
            matnd[i][j] = 1.;
            for (int k = 0; k < n_vars; ++k) {
                matnd[i][j] *= mat1d[m1[k]][m2[k]];
            }
        }
    }

    this->hermiteCoeffs_1d = mat1d;
    this->hermiteCoeffs_nd = matnd;
}

int Solver_spectral::mult2ind(vector<int> m, int d) {
    int l = m.size() - 1; int i;
    for(i = l; (i > 0) & (m[i] == 0); i--);

    if ((i == 0) & (m[0] == 0))
        return 0;

    int s = 0;
    for (unsigned int j = 0; j < m.size(); ++j)
        s += m[j];

    int dr = d - s;
    int vr = l - i;
    m[i] = m[i] - 1;
    return bin(dr + vr + 1, vr) + mult2ind(m, d);
}

vector<int> Solver_spectral::ind2mult(int ind, int d, int n) {
    vector<int> m(n,0); int s = 0;
    for (int i = 0; i < ind; ++i) {
        if (s < d) {
            m[n-1] ++; s++;
        } else {
            int j; for(j = n-1; m[j] == 0; j--);
            s -= m[j]; m[j] = 0; m[j-1] ++; s ++;
        }
    }
    return m;
}

vector<double> Solver_spectral::mon2herm (vector<double> mcoeffs, int n, int d) {
    vector<double> result(mcoeffs.size(), 0.);
    for (unsigned int i = 0; i < mcoeffs.size(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            result[i] += hermiteCoeffs_nd[i][j] * mcoeffs[j];
        }
    }
    return result;
}
