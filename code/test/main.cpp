#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <vector>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <time.h>
#include <functional>
#include <stack>
#include <typeinfo>

#define PI 3.141592653589793238463

// Namespace
using namespace std;

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

void printMat(vector< vector<double> > a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void printVec(vector<double> a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
            cout << a[i] << " ";
            cout << endl;
    }
}

vector<double> solve(vector< vector<double> > A, vector<double> b) {
    vector<double> result = b;
    vector< vector<double> > L = cholesky(A);
    printMat(L);
    vector< vector<double> > U = transpose(L);
    printMat(U);
    result = lsolve(L, result);
    printVec(result);
    result = usolve(U, result);
    printVec(result);
    return result;
}

int main(int argc, char *argv[])
{
    vector<double> a = {10, 3, 1};
    vector<double> b = {0, 15, -1};
    vector<double> c = {4, 0, 20};
    vector< vector<double> > A = {a,b,c};

    vector<double> r = {3, 2, 1};
    vector<double> result = solve(A,r);

    printMat(A);

    cout << result[0] << " " << result[1] << " " << result[2] << endl;

    return 0;
}
