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

// Normalized Hermite polynomials
double hermite(int n, double x, double sigma) {
    double toReturn = 0.;

/* FIXME: this is a trick (urbain, Wed 29 Apr 2015 11:29:39 BST) */

    if (abs(sigma) > 1E-12)
        x = x/sigma;
    else {
        cout << "Variance has to be positive, but sigma = " << sigma <<  endl;
        exit(0);
    }
    switch (n) {
        case 0: toReturn = 1.0; break;
        case 1: toReturn = x; break;
        case 2: toReturn = sqrt(2.0)*(x*x-1.0)*(1.0/2.0); break;
        case 3: toReturn = sqrt(6.0)*x*(x*x-3.0)*(1.0/6.0); break;
        case 4: toReturn = sqrt(6.0)*((x*x)*-6.0+x*x*x*x+3.0)*(1.0/1.2E1); break;
        case 5: toReturn = sqrt(3.0E1)*x*((x*x)*-1.0E1+x*x*x*x+1.5E1)*(1.0/6.0E1); break;
        case 6: toReturn = sqrt(5.0)*((x*x)*4.5E1-(x*x*x*x)*1.5E1+x*x*x*x*x*x-1.5E1)*(1.0/6.0E1); break;
        case 7: toReturn = sqrt(3.5E1)*x*((x*x)*1.05E2-(x*x*x*x)*2.1E1+x*x*x*x*x*x-1.05E2)*(1.0/4.2E2); break;
        case 8: toReturn = sqrt(7.0E1)*((x*x)*-4.2E2+(x*x*x*x)*2.1E2-(x*x*x*x*x*x)*2.8E1+x*x*x*x*x*x*x*x+1.05E2)*5.952380952380952E-4; break;
        case 9: toReturn = sqrt(7.0E1)*x*((x*x)*-1.26E3+(x*x*x*x)*3.78E2-(x*x*x*x*x*x)*3.6E1+x*x*x*x*x*x*x*x+9.45E2)*1.984126984126984E-4; break;
        case 10: toReturn = sqrt(7.0)*((x*x)*4.725E3-(x*x*x*x)*3.15E3+(x*x*x*x*x*x)*6.3E2-(x*x*x*x*x*x*x*x)*4.5E1+pow(x,1.0E1)-9.45E2)*1.984126984126984E-4; break;
        case 11: toReturn = sqrt(7.7E1)*x*((x*x)*1.7325E4-(x*x*x*x)*6.93E3+(x*x*x*x*x*x)*9.9E2-(x*x*x*x*x*x*x*x)*5.5E1+pow(x,1.0E1)-1.0395E4)*1.803751803751804E-5; break;
        case 12: toReturn = sqrt(2.31E2)*((x*x)*-6.237E4+(x*x*x*x)*5.1975E4-(x*x*x*x*x*x)*1.386E4+(x*x*x*x*x*x*x*x)*1.485E3-pow(x,1.0E1)*6.6E1+pow(x,1.2E1)+1.0395E4)*3.006253006253006E-6; break;
        case 13: toReturn = 5.479963503528103E1*x*((x*x)*-2.7027E5+(x*x*x*x)*1.35135E5-(x*x*x*x*x*x)*2.574E4+(x*x*x*x*x*x*x*x)*2.145E3-pow(x,1.0E1)*7.8E1+pow(x,1.2E1)+1.35135E5)*2.312502312502313E-7; break;
        case 14: toReturn = sqrt(8.58E2)*((x*x)*9.45945E5-(x*x*x*x)*9.45945E5+(x*x*x*x*x*x)*3.15315E5-(x*x*x*x*x*x*x*x)*4.5045E4+pow(x,1.0E1)*3.003E3-pow(x,1.2E1)*9.1E1+pow(x,1.4E1)-1.35135E5)*1.156251156251156E-7; break;
        case 15: toReturn = 3.781534080237807E1*x*((x*x)*4.729725E6-(x*x*x*x)*2.837835E6+(x*x*x*x*x*x)*6.75675E5-(x*x*x*x*x*x*x*x)*7.5075E4+pow(x,1.0E1)*4.095E3-pow(x,1.2E1)*1.05E2+pow(x,1.4E1)-2.027025E6)*2.312502312502313E-8; break;
        default: cout << "Degree too high" << endl; exit(0);
    }
    return toReturn;
}


// Multidimensional Hermite
double hermiteM(vector<int> multIndex, vector<double> x, vector<double> sigmas) {
    double h_eval = 1.;
    for (unsigned int i = 0; i < multIndex.size(); ++i) {
        h_eval *= hermite(multIndex[i],x[i],sigmas[i]);
    }
    return h_eval;
}

// Canonical integer associated with a multindex
int canonicalInd(vector<int> alpha, int n, int degree) {
    int toReturn = 0;
    for (int j = 0; j < n; ++j) {
        toReturn += alpha[j]*pow(degree + 1, n -j - 1);
    }
    return toReturn;
}

