#include "aux.hpp"

using namespace std;

// Delta function
double delta(int a, int b) {
    if (fabs(b-a) < 0.1)
        return 1.;
    else return 0.;
}

// Symmetric part of a matrix
vector< vector<double> > symmetric(vector< vector<double> > A) {
    int n = A.size();
    vector< vector <double> > result(n,vector<double>(n,0.));

    for (int i = 0.; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            result[i][j] = 0.5*(A[i][j] + A[j][i]);
            result[j][i] = result[i][j];
        }
    }
    return result;
}

// Cholesky factorization of a matrix
vector< vector<double> > cholesky(vector< vector<double> > A) {
    int n = A.size();
    vector< vector <double> > L(n,vector<double>(n,0.));

    for (int i = 0.; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.;
            for (int k = 0; k < j; k++) {
                sum += L[i][k]*L[j][k];
            }
            if (i == j) {
                L[i][i] = sqrt(A[i][i]-sum);
            }
            else {
                if(fabs(L[j][j])>1e-14) {
                    L[i][j] = (A[i][j] - sum)/L[j][j];
                }
                else {
                    L[i][j] = 0.;
                }
            }
        }
    }
    return L;
}

// Binomial coefficients
int bin(int n, int k) {
    int res = 1;

    if ( k > n - k )
        k = n - k;

    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

int mult2ind(vector<int> m, int d) {
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

vector<int> ind2mult(int ind, int d, int n) {
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

double ipow(double x, int e) {
    if(e == 0) return 1;
    if (e == 1) return x;
    double aux = ipow(x,e/2);
    if (e%2 == 0) return aux*aux;
    return x*aux*aux;
}

double monomial(vector<int> mult, vector<double> x, vector<double> sigmas) {
    double result = 1.;
    for (unsigned int i = 0; i < mult.size(); ++i) {
        result *= ipow(x[i]/sigmas[i], mult[i]);
    }
    return result;
}

vector< vector<double> > hermiteCoeffs(int degree) {

    vector< vector<double> > coefficients(degree + 1, vector<double>(degree + 1,0.));

    coefficients[0][0] = 1.0;
    if (degree == 0) return coefficients;

    coefficients[1][0] = 0.0;
    coefficients[1][1] = 1.0;
    if (degree == 1) return coefficients;

    coefficients[2][0] = sqrt(2.0)/2.0*(-1.0);
    coefficients[2][1] = sqrt(2.0)/2.0*( 0.0);
    coefficients[2][2] = sqrt(2.0)/2.0*( 1.0);
    if (degree == 2) return coefficients;

    coefficients[3][0] = sqrt(6.0)/6.0*( 0.0);
    coefficients[3][1] = sqrt(6.0)/6.0*(-3.0);
    coefficients[3][2] = sqrt(6.0)/6.0*( 0.0);
    coefficients[3][3] = sqrt(6.0)/6.0*( 1.0);
    if (degree == 3) return coefficients;

    coefficients[4][0] = sqrt(6.0)/12.0*( 3.0);
    coefficients[4][1] = sqrt(6.0)/12.0*( 0.0);
    coefficients[4][2] = sqrt(6.0)/12.0*(-6.0);
    coefficients[4][3] = sqrt(6.0)/12.0*( 0.0);
    coefficients[4][4] = sqrt(6.0)/12.0*( 1.0);
    if (degree == 4) return coefficients;

    coefficients[5][0] = sqrt(30.0)/60.0*(  0.0);
    coefficients[5][1] = sqrt(30.0)/60.0*( 15.0);
    coefficients[5][2] = sqrt(30.0)/60.0*(  0.0);
    coefficients[5][3] = sqrt(30.0)/60.0*(-10.0);
    coefficients[5][4] = sqrt(30.0)/60.0*(  0.0);
    coefficients[5][5] = sqrt(30.0)/60.0*(  1.0);
    if (degree == 5) return coefficients;

    coefficients[6][0] = sqrt(5.0)/60.0*(-15.0);
    coefficients[6][1] = sqrt(5.0)/60.0*(  0.0);
    coefficients[6][2] = sqrt(5.0)/60.0*( 45.0);
    coefficients[6][3] = sqrt(5.0)/60.0*(  0.0);
    coefficients[6][4] = sqrt(5.0)/60.0*(-15.0);
    coefficients[6][5] = sqrt(5.0)/60.0*(  0.0);
    coefficients[6][6] = sqrt(5.0)/60.0*(  1.0);
    if (degree == 6) return coefficients;

    coefficients[7][0] = sqrt(35.0)/420.0*(   0.0);
    coefficients[7][1] = sqrt(35.0)/420.0*(-105.0);
    coefficients[7][2] = sqrt(35.0)/420.0*(   0.0);
    coefficients[7][3] = sqrt(35.0)/420.0*( 105.0);
    coefficients[7][4] = sqrt(35.0)/420.0*(   0.0);
    coefficients[7][5] = sqrt(35.0)/420.0*( -21.0);
    coefficients[7][6] = sqrt(35.0)/420.0*(   0.0);
    coefficients[7][7] = sqrt(35.0)/420.0*(   1.0);
    if (degree == 7) return coefficients;

    coefficients[8][0] = sqrt(70.0)/1680.0*( 105.0);
    coefficients[8][1] = sqrt(70.0)/1680.0*(   0.0);
    coefficients[8][2] = sqrt(70.0)/1680.0*(-420.0);
    coefficients[8][3] = sqrt(70.0)/1680.0*(   0.0);
    coefficients[8][4] = sqrt(70.0)/1680.0*( 210.0);
    coefficients[8][5] = sqrt(70.0)/1680.0*(   0.0);
    coefficients[8][6] = sqrt(70.0)/1680.0*( -28.0);
    coefficients[8][7] = sqrt(70.0)/1680.0*(   0.0);
    coefficients[8][8] = sqrt(70.0)/1680.0*(   1.0);
    if (degree == 8) return coefficients;

    coefficients[9][0] = sqrt(70.0)/5040.0*(    0.0);
    coefficients[9][1] = sqrt(70.0)/5040.0*(  945.0);
    coefficients[9][2] = sqrt(70.0)/5040.0*(    0.0);
    coefficients[9][3] = sqrt(70.0)/5040.0*(-1260.0);
    coefficients[9][4] = sqrt(70.0)/5040.0*(    0.0);
    coefficients[9][5] = sqrt(70.0)/5040.0*(  378.0);
    coefficients[9][6] = sqrt(70.0)/5040.0*(    0.0);
    coefficients[9][7] = sqrt(70.0)/5040.0*(  -36.0);
    coefficients[9][8] = sqrt(70.0)/5040.0*(    0.0);
    coefficients[9][9] = sqrt(70.0)/5040.0*(    1.0);
    if (degree == 9) return coefficients;

    coefficients[10][0]  = sqrt(7.0)/5040.0*( -945.0);
    coefficients[10][1]  = sqrt(7.0)/5040.0*(    0.0);
    coefficients[10][2]  = sqrt(7.0)/5040.0*( 4725.0);
    coefficients[10][3]  = sqrt(7.0)/5040.0*(    0.0);
    coefficients[10][4]  = sqrt(7.0)/5040.0*(-3150.0);
    coefficients[10][5]  = sqrt(7.0)/5040.0*(    0.0);
    coefficients[10][6]  = sqrt(7.0)/5040.0*(  630.0);
    coefficients[10][7]  = sqrt(7.0)/5040.0*(    0.0);
    coefficients[10][8]  = sqrt(7.0)/5040.0*(  -45.0);
    coefficients[10][9]  = sqrt(7.0)/5040.0*(    0.0);
    coefficients[10][10] = sqrt(7.0)/5040.0*(    1.0);
    if (degree == 10) return coefficients;

    cout << "Degree too high" << endl;
    exit(0);
}

vector<vector<double> > hermiteCoeffs_nd(int d, int n) {
    int nb = bin(d + n, n);
    vector< vector<double> > mat1d = hermiteCoeffs(d);
    vector< vector<double> > matnd(nb, vector<double>(nb,0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i,d,n);
        for (int j = 0; j < nb; ++j) {
        vector<int> m2 = ind2mult(j,d,n);
            matnd[i][j] = 1.;
            for (int k = 0; k < n; ++k) {
                matnd[i][j] *= mat1d[m1[k]][m2[k]];
            }
        }
    }
    return matnd;
}

vector<double> hcoeffs (vector<double> mcoeffs, int n, int d) {
    vector<double> result(mcoeffs.size(), 0.);
    vector< vector<double> > mat = hermiteCoeffs_nd(d, n);
    for (unsigned int i = 0; i < mcoeffs.size(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            result[i] += mat[i][j] * mcoeffs[j];
        }
    }
    return result;
}

// Canonical integer associated with a multindex
int canonicalInd(vector<int> alpha, int n, int degree) {
    int toReturn = 0;
    for (int j = 0; j < n; ++j) {
        toReturn += alpha[j]*pow(degree + 1, n -j - 1);
    }
    return toReturn;
}
