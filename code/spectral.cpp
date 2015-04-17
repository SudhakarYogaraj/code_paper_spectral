#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <time.h>

// Namespace
using namespace std;

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

// Functions of the problem
double slow_drift(double x, double y) {
    return sin(x)*cos(y);
}
double h(double x, double y) {
    return cos(x)*cos(y);
}

// Physicists' Hermite polynomials
double hermite(int n, double x, double sigma) {    
    double toReturn = 0.;

    // Scaling in case the variance is different from 1.
    if (abs(sigma) > 1E-16)
        x = x/sigma;
    else {
        cout << "Please provide a non-negative variance" << endl;
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
        case 16: toReturn = 3.781534080237807E1*((x*x)*-1.62162E7+(x*x*x*x)*1.89189E7-(x*x*x*x*x*x)*7.56756E6+(x*x*x*x*x*x*x*x)*1.35135E6-pow(x,1.0E1)*1.2012E5+pow(x,1.2E1)*5.46E3-pow(x,1.4E1)*1.2E2+pow(x,1.6E1)+2.027025E6)*5.781255781255781E-9; break;
        case 17: toReturn = 1.559166443969341E2*x*((x*x)*-9.18918E7+(x*x*x*x)*6.432426E7-(x*x*x*x*x*x)*1.837836E7+(x*x*x*x*x*x*x*x)*2.55255E6-pow(x,1.0E1)*1.8564E5+pow(x,1.2E1)*7.14E3-pow(x,1.4E1)*1.36E2+pow(x,1.6E1)+3.4459425E7)*3.400738694856342E-10; break;
        case 18: toReturn = 1.102497165529236E2*((x*x)*3.10134825E8-(x*x*x*x)*4.135131E8+(x*x*x*x*x*x)*1.9297278E8-(x*x*x*x*x*x*x*x)*4.135131E7+pow(x,1.0E1)*4.59459E6-pow(x,1.2E1)*2.7846E5+pow(x,1.4E1)*9.18E3-pow(x,1.6E1)*1.53E2+pow(x,1.8E1)-3.4459425E7)*1.133579564952114E-10; break;
        case 19: toReturn = 4.805673730081975E2*x*((x*x)*1.964187225E9-(x*x*x*x)*1.57134978E9+(x*x*x*x*x*x)*5.2378326E8-(x*x*x*x*x*x*x*x)*8.729721E7+pow(x,1.0E1)*7.93611E6-pow(x,1.2E1)*4.0698E5+pow(x,1.4E1)*1.1628E4-pow(x,1.6E1)*1.71E2+pow(x,1.8E1)-6.54729075E8)*5.966208236590074E-12; break;
        case 20: toReturn = 2.149162627629654E2*((x*x)*-6.54729075E9+(x*x*x*x)*9.820936125E9-(x*x*x*x*x*x)*5.2378326E9+(x*x*x*x*x*x*x*x)*1.30945815E9-pow(x,1.0E1)*1.7459442E8+pow(x,1.2E1)*1.322685E7-pow(x,1.4E1)*5.814E5+pow(x,1.6E1)*1.4535E4-pow(x,1.8E1)*1.9E2+pow(x,2.0E1)+6.54729075E8)*2.983104118295037E-12; break;
        default: cout << "Degree too high" << endl; exit(0);
    }
    return toReturn;
}

// Main method
int main(int argc, char *argv[]) {

    // Degree of polynomials approximation.
    int degree = 10;

    // Macro time step
    double Dt = .1;

    // Number of fast variables
    int nf = 3;

    // Eigenvalues
    vector<double> lambdas = {1.,2.,4.};
    vector<double> qs      = {2.,4.,2.};
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        sigmas[i] = sqrt(0.5*qs[i]*qs[i]/lambdas[i]);
    }

    // Number of polynomials in the basis
    int nBasis = bin(degree + nf, nf);

    // Relation linear index - multiindex
    vector< vector<int> > indexMap(nBasis, vector<int>(nf,0));

    vector<int> currentMult(nf,0);
    currentMult[nf-1] = 0;
    int movingIndex;
    for (int i = 0; i < nBasis; ++i) {
        for (int j = 0; j < nf; ++j) {
            indexMap[i][j] = currentMult[j];
        }
        int sum = 0;
        for (int j = 0; j < nf; ++j) {
            sum += currentMult[j];
        }
        if (sum < degree) {
            currentMult[nf-1] += 1;
        } else {
            int auxIndex = nf - 1;
            while (currentMult[auxIndex] == 0) {
                auxIndex --;
            }
            currentMult[auxIndex] = 0;
            currentMult[auxIndex-1] ++;
        } 
    }
    for (int i = 0; i < nBasis; ++i) {
       cout << "Element number " << i << ": ";
       for (int j = 0; j < nf; ++j) {
           cout << indexMap[i][j] << ", ";
       }
       cout << endl;
    }

    // Final time
    double T = 1;

    // Vector of coefficients
    vector<double> coefficients(degree+1, 0.);

    // Vector of times
    int sizet = (int) (T/Dt) + 1;
    vector<double> t(sizet,0.);

    // Construction of the vector containing the times at which
    // the slow variable is calculated
    for (int i = 0; i < sizet; i++) {
        t[i] = i*Dt;
    }

    // Random variable for generator
    int seed = time(NULL);

    // Random number generator
    default_random_engine generator; generator.seed(seed);
    normal_distribution<double> distribution(0.0,1.0);

    // Slow variable
    vector<double> x(sizet,0.2);

    // Macro loop
    for (int i = 0; i < sizet; ++i) {
        
        cout << "Time of simulation: " << t[i] << endl;

        // Solution of the cell problem

        // Computation of the coefficients
        for (int n = 0; n <= nBasis; ++n) {

            // Monte Carlo
            int N_mc = 100000;
            double sum = 0.;
            for (int k = 0; k < N_mc; ++k) {
                double randn = distribution(generator);
                sum += hermite(n,randn,1.)*slow_drift(x[i],randn);   
            }
            coefficients[n] = sum/N_mc;
            cout << "Value of coefficient " << coefficients[n] << endl;
        }
        cout << endl;
    }
}
