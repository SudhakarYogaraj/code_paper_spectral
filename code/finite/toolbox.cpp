#include "toolbox.hpp"
#include "templates.hpp"

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
                if (A[i][i] - sum < -1e-14) {
                    cout << "Warning: matrix is not positive definite: (" << A[i][i] - sum << ")." << endl;
                    if (A[i][i] - sum < -1e-12) {
                        cout << "Error occured during Cholesky factorization of line " << i << "." << endl;
                        exit(0);
                    }
                }
                else {
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

void lu(vector< vector<double> > a, vector< vector<double> >& l, vector< vector<double> >& u) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < a.size(); i++)
    {
        for (j = 0; j < a.size(); j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < a.size(); j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
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

vector<double> solve(vector< vector<double> > A, vector<double> b) {
    vector<double> result = b;
    vector< vector<double> > At = transpose(A);
    for (unsigned int i = 0; i < A.size(); ++i) {
        for (unsigned int j = 0; j < A.size(); j++) {
            if ( fabs(A[i][j] - A[j][i]) > 1e-14) {
                cout << "Warning: matrix must be symmetric, but A[" << i << "][" << j << "] - A[" << j << "][" << i << "] = " << A[i][j] - A[j][i] << "." << endl;
            }
        }
    }

    vector< vector<double> > L(A.size(), vector<double> (A.size(), 0.));
    vector< vector<double> > U(A.size(), vector<double> (A.size(), 0.));
    lu(A,L,U);
    /* vector< vector<double> > L = cholesky(A); */
    /* vector< vector<double> > U = transpose(L); */
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

void qr(vector< vector<double> > a, vector< vector<double> >& q, vector< vector<double> >& r) {
    vector< vector<double> > at = transpose(a);
    vector< vector<double> > qt (a.size(), vector<double> (a.size(), 0.));
    for (unsigned int i = 0; i < a.size(); ++i) {
        vector<double> u = at[i];
        for (unsigned int j = 0; j < i; ++j) {
            u = u - qt[j] * (at[i] * qt[j]);
            r[i][j] = qt[j] * at[i];
        }
        qt[i] = u * (1/fabs(u));
        r[i][i] = qt[i] * at[i];

    }
    q = transpose(qt);

    if ( fabs(a - q*r) > 1e-10) {
        cout << "Error in qr decomposition" << endl;
        exit(0);
    }
}

void eig_qr(vector< vector<double> > a, vector< vector<double> >& v, vector<double>& l) {
    vector< vector<double> > q (a.size(), vector<double> (a.size(), 0.));
    vector< vector<double> > r (a.size(), vector<double> (a.size(), 0.));

    int n_iter = 100;
    for (int n = 0; n < n_iter; ++n) {
        qr(a,q,r);
        a = r*q;
    }

    for (int i = 0; i < a.size(); ++i) {
        l[i] = a[i][i];
    }
    v = q;
}
