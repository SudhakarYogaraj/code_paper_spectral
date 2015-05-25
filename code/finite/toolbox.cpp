#include "toolbox.hpp"

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
                if (A[i][i] - sum < 0) {
                    cout << "Warning: matrix is not positive definite: (" << A[i][i] - sum << ")." << endl;
                    if (A[i][i] - sum < -1e-12) {
                        cout << "Error occured during Cholesky factorization of line " << i << "." << endl;
                        exit(0);
                    }
                }
                else {
                    /* cout << i << ". " << A[i][i] - sum << endl; */
                    L[i][i] = sqrt(A[i][i]-sum);
                }
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
    vector< vector<double> > At = transpose(A);
    for (unsigned int i = 0; i < A.size(); ++i) {
        for (unsigned int j = 0; j < A.size(); j++) {
            if ( fabs(A[i][j] - A[j][i]) > 1e-10) {
                cout << "Matrix must be symmetric, but A[" << i << "][" << j << "] - A[" << j << "][" << i << "] = " << A[i][j] - A[j][i] << "." << endl;
                exit(0);
            }
        }
    }
    /* for (unsigned int i = 0; i < A.size(); ++i) { */
    /*     for (unsigned int j = 0; j < i; j++) { */
    /*         double tmp = (A[i][j] + A[j][i])/2.; */
    /*         A[i][j] = tmp; */
    /*         A[j][i] = tmp; */
    /*     } */
    /* } */
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
