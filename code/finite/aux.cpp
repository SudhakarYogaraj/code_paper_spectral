#include "aux.hpp"

using namespace std;

// Standard normal gaussian
double gaussian(vector<double> y, vector<double> sigmas) {
    double result = 1.;
    for (unsigned int i = 0; i < y.size(); ++i) {
        double s = sigmas[i];
        result *= exp(-y[i]*y[i]/(2*s*s))/(sqrt(2*PI)*s);
    }
    return result;
}

double ipow(double x, int e) {
    if(e == 0) return 1;
    if (e == 1) return x;
    double aux = ipow(x,e/2);
    if (e%2 == 0) return aux*aux;
    return x*aux*aux;
}

double laplacian(vector<double> mult, vector<double> y) {
    double result = 0.;
    for (unsigned int i = 0; i < mult.size(); ++i) {
        if (mult[i] >= 2) {
            double result_i = 1.;
            vector<double> newMult = mult;
            newMult[i] -= 2;
            for (unsigned int j = 0; j < mult.size(); ++j) {
                result_i *= ipow(y[i], mult[i]);
            }
            result += result_i;
        }
    }
    return result;
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

vector< vector<double> > transpose(vector< vector<double> > A)
{
    int n = A.size();
    vector< vector <double> > T(n,vector<double>(n,0.));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            T[i][j] = A[j][i];
        }
    }
    return T;
}

// Cholesky factorization of a matrix
vector<double> lsolve(vector< vector<double> > U, vector<double> b) {
    int n = U.size();
    vector<double> result(n, 0.);

    for (int i = 0.; i < n; i++) {
        double tmp = b[i];
        for (int j = 0; j < i; ++j) {
            tmp -= U[i][j]*result[j];
        }
        result[i] = tmp/U[i][i];
    }
    return result;
}

// Cholesky factorization of a matrix
vector<double> usolve(vector< vector<double> > U, vector<double> b) {
    int n = U.size();
    vector<double> result(n, 0.);

    for (int i = n-1; i >= 0; i--) {
        double tmp = b[i];
        for (int j = n-1; j > i; --j) {
            tmp -= U[i][j]*result[j];
        }
        result[i] = tmp/U[i][i];
    }
    return result;
}

void niceMat(vector< vector<double> > a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void niceVec(vector<double> a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
            cout << a[i] << " ";
            cout << endl;
    }
}

vector<double> solve(vector< vector<double> > A, vector<double> b) {
    vector<double> result = b;
    vector< vector<double> > L = cholesky(A);
    vector< vector<double> > U = transpose(L);
    result = lsolve(L, result);
    result = usolve(U, result);
    return result;
}

// Binomial coefficients
int bin(int n, int k) {
    int res = 1;
    if ( k > n - k )
        k = n - k;
    for (int i = 0; i < k; ++i) {
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

double monomial(vector<int> mult, vector<double> x, vector<double> sigmas) {
    double result = 1.;
    for (unsigned int i = 0; i < mult.size(); ++i) {
        result *= ipow(x[i]/sigmas[i], mult[i]);
    }
    return result;
}

double monomial(vector<int> mult, vector<double> x) {
    double result = 1.;
    for (unsigned int i = 0; i < mult.size(); ++i) {
        result *= ipow(x[i], mult[i]);
    }
    return result;
}

vector< vector<double> > hermiteCoeffs(int degree) {

    vector< vector<double> > coefficients(degree + 1, vector<double>(degree + 1,0.));

    coefficients[0][0] = 1.0;
    if (degree == 0) return coefficients;

    coefficients[1][1] = 1.0;
    if (degree == 1) return coefficients;

    coefficients[2][0] = sqrt(2.0)/2.0*(-1.0);
    coefficients[2][2] = sqrt(2.0)/2.0*( 1.0);
    if (degree == 2) return coefficients;

    coefficients[3][1] = sqrt(6.0)/6.0*(-3.0);
    coefficients[3][3] = sqrt(6.0)/6.0*( 1.0);
    if (degree == 3) return coefficients;

    coefficients[4][0] = sqrt(6.0)/12.0*( 3.0);
    coefficients[4][2] = sqrt(6.0)/12.0*(-6.0);
    coefficients[4][4] = sqrt(6.0)/12.0*( 1.0);
    if (degree == 4) return coefficients;

    coefficients[5][1] = sqrt(30.0)/60.0*( 15.0);
    coefficients[5][3] = sqrt(30.0)/60.0*(-10.0);
    coefficients[5][5] = sqrt(30.0)/60.0*(  1.0);
    if (degree == 5) return coefficients;

    coefficients[6][0] = sqrt(5.0)/60.0*(-15.0);
    coefficients[6][2] = sqrt(5.0)/60.0*( 45.0);
    coefficients[6][4] = sqrt(5.0)/60.0*(-15.0);
    coefficients[6][6] = sqrt(5.0)/60.0*(  1.0);
    if (degree == 6) return coefficients;

    coefficients[7][1] = sqrt(35.0)/420.0*(-105.0);
    coefficients[7][3] = sqrt(35.0)/420.0*( 105.0);
    coefficients[7][5] = sqrt(35.0)/420.0*( -21.0);
    coefficients[7][7] = sqrt(35.0)/420.0*(   1.0);
    if (degree == 7) return coefficients;

    coefficients[8][0] = sqrt(70.0)/1680.0*( 105.0);
    coefficients[8][2] = sqrt(70.0)/1680.0*(-420.0);
    coefficients[8][4] = sqrt(70.0)/1680.0*( 210.0);
    coefficients[8][6] = sqrt(70.0)/1680.0*( -28.0);
    coefficients[8][8] = sqrt(70.0)/1680.0*(   1.0);
    if (degree == 8) return coefficients;

    coefficients[9][1] = sqrt(70.0)/5040.0*(  945.0);
    coefficients[9][3] = sqrt(70.0)/5040.0*(-1260.0);
    coefficients[9][5] = sqrt(70.0)/5040.0*(  378.0);
    coefficients[9][7] = sqrt(70.0)/5040.0*(  -36.0);
    coefficients[9][9] = sqrt(70.0)/5040.0*(    1.0);
    if (degree == 9) return coefficients;

    coefficients[10][0]  = sqrt(7.0)/5040.0*( -945.0);
    coefficients[10][2]  = sqrt(7.0)/5040.0*( 4725.0);
    coefficients[10][4]  = sqrt(7.0)/5040.0*(-3150.0);
    coefficients[10][6]  = sqrt(7.0)/5040.0*(  630.0);
    coefficients[10][8]  = sqrt(7.0)/5040.0*(  -45.0);
    coefficients[10][10] = sqrt(7.0)/5040.0*(    1.0);
    if (degree == 10) return coefficients;

    coefficients[11][11] = (sqrt(77.))/55440.;
    coefficients[11][9]  = - (sqrt(77.))/1008.;
    coefficients[11][7]  =  (sqrt(77.))/56.;
    coefficients[11][5]  = - (sqrt(77.))/8.;
    coefficients[11][3]  =  (5.*sqrt(77.))/16.;
    coefficients[11][1]  = - (3.*sqrt(77.))/16.;
    if (degree == 11) return coefficients;

    coefficients[12][12] = (sqrt(231.))/332640.;
    coefficients[12][10] = - (sqrt(231.))/5040.;
    coefficients[12][8]  =  (sqrt(231.))/224.;
    coefficients[12][6]  = - (sqrt(231.))/24.;
    coefficients[12][4]  =  (5.*sqrt(231.))/32.;
    coefficients[12][2]  = - (3.*sqrt(231.))/16.;
    coefficients[12][0]  =  sqrt(231.)/32.;
    if (degree == 12) return coefficients;

    coefficients[13][13] = (sqrt(3003.))/4324320.;
    coefficients[13][11] = - (sqrt(3003.))/55440.;
    coefficients[13][9]  =  (sqrt(3003.))/2016.;
    coefficients[13][7]  = - (sqrt(3003.))/168.;
    coefficients[13][5]  =  (sqrt(3003.))/32.;
    coefficients[13][3]  = - (sqrt(3003.))/16.;
    coefficients[13][1]  =  (sqrt(3003.))/32.;
    if (degree == 13) return coefficients;

    coefficients[14][14] = (sqrt(858.))/8648640.;
    coefficients[14][12] = - (sqrt(858.))/95040.;
    coefficients[14][10] =  (sqrt(858.))/2880.;
    coefficients[14][8]  = - (sqrt(858.))/192.;
    coefficients[14][6]  =  (7.*sqrt(858.))/192.;
    coefficients[14][4]  = - (7.*sqrt(858.))/64.;
    coefficients[14][2]  =  (7.*sqrt(858.))/64.;
    coefficients[14][0]  = - sqrt(858.)/64.;
    if (degree == 14) return coefficients;

    coefficients[15][15] = (sqrt(1430.))/43243200.;
    coefficients[15][13] = - (sqrt(1430.))/411840.;
    coefficients[15][11] =  (sqrt(1430.))/10560.;
    coefficients[15][9]  = - (sqrt(1430.))/576.;
    coefficients[15][7]  =  (sqrt(1430.))/64.;
    coefficients[15][5]  = - (21.*sqrt(1430.))/320.;
    coefficients[15][3]  =  (7.*sqrt(1430.))/64.;
    coefficients[15][1]  = - (3.*sqrt(1430.))/64.;
    if (degree == 15) return coefficients;

    coefficients[16][16] = (sqrt(1430.))/172972800.;
    coefficients[16][14] = - (sqrt(1430.))/1441440.;
    coefficients[16][12] =  (sqrt(1430.))/31680.;
    coefficients[16][10] = - (sqrt(1430.))/1440.;
    coefficients[16][8]  =  (sqrt(1430.))/128.;
    coefficients[16][6]  = - (7.*sqrt(1430.))/160.;
    coefficients[16][4]  =  (7.*sqrt(1430.))/64.;
    coefficients[16][2]  = - (3.*sqrt(1430.))/32.;
    coefficients[16][0]  =  (3.*sqrt(1430.))/256.;
    if (degree == 16) return coefficients;

    coefficients[17][17] = (sqrt(24310.))/2940537600.;
    coefficients[17][15] = - (sqrt(24310.))/21621600.;
    coefficients[17][13] =  (sqrt(24310.))/411840.;
    coefficients[17][11] = - (sqrt(24310.))/15840.;
    coefficients[17][9]  =  (sqrt(24310.))/1152.;
    coefficients[17][7]  = - (sqrt(24310.))/160.;
    coefficients[17][5]  =  (7.*sqrt(24310.))/320.;
    coefficients[17][3]  = - (sqrt(24310.))/32.;
    coefficients[17][1]  =  (3.*sqrt(24310.))/256.;
    if (degree == 17) return coefficients;

    coefficients[18][18] = (sqrt(12155.))/8821612800.;
    coefficients[18][16] = - (sqrt(12155.))/57657600.;
    coefficients[18][14] =  (sqrt(12155.))/960960.;
    coefficients[18][12] = - (sqrt(12155.))/31680.;
    coefficients[18][10] =  (sqrt(12155.))/1920.;
    coefficients[18][8]  = - (3.*sqrt(12155.))/640.;
    coefficients[18][6]  =  (7.*sqrt(12155.))/320.;
    coefficients[18][4]  = - (3.*sqrt(12155.))/64.;
    coefficients[18][2]  =  (9.*sqrt(12155.))/256.;
    coefficients[18][0]  = - sqrt(12155.)/256.;
    if (degree == 18) return coefficients;

    coefficients[19][19] = (sqrt(230945.))/167610643200.;
    coefficients[19][17] = - (sqrt(230945.))/980179200.;
    coefficients[19][15] =  (sqrt(230945.))/14414400.;
    coefficients[19][13] = - (sqrt(230945.))/411840.;
    coefficients[19][11] =  (sqrt(230945.))/21120.;
    coefficients[19][9]  = - (sqrt(230945.))/1920.;
    coefficients[19][7]  =  (sqrt(230945.))/320.;
    coefficients[19][5]  = - (3.*sqrt(230945.))/320.;
    coefficients[19][3]  =  (3.*sqrt(230945.))/256.;
    coefficients[19][1]  = - (sqrt(230945.))/256.;
    if (degree == 19) return coefficients;

    coefficients[20][20] = (sqrt(46189.))/335221286400.;
    coefficients[20][18] = - (sqrt(46189.))/1764322560.;
    coefficients[20][16] =  (sqrt(46189.))/23063040.;
    coefficients[20][14] = - (sqrt(46189.))/576576.;
    coefficients[20][12] =  (sqrt(46189.))/25344.;
    coefficients[20][10] = - (sqrt(46189.))/1920.;
    coefficients[20][8]  =  (sqrt(46189.))/256.;
    coefficients[20][6]  = - (sqrt(46189.))/64.;
    coefficients[20][4]  =  (15.*sqrt(46189.))/512.;
    coefficients[20][2]  = - (5.*sqrt(46189.))/256.;
    coefficients[20][0]  =  sqrt(46189.)/512.;
    if (degree == 20) return coefficients;

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

vector<double> mon2herm (vector<double> mcoeffs, int n, int d) {
    vector<double> result(mcoeffs.size(), 0.);
    vector< vector<double> > mat = hermiteCoeffs_nd(d, n);
    for (unsigned int i = 0; i < mcoeffs.size(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            result[i] += mat[i][j] * mcoeffs[j];
        }
    }
    return result;
}
