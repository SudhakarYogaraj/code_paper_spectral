/* #include <iostream> */
/* #include <vector> */
/* #include <cmath> */
/* #include <random> */
/* #include <algorithm> */
/* #include <iterator> */
/* #include <fstream> */
/* #include <time.h> */

/* /1* TODO: Take problem as input (urbain, Sun 19 Apr 2015 13:45:33 BST) *1/ */

/* // Namespace */
/* using namespace std; */

/* // Canonical integer associated with a multindex */
/* int canonicalInd(vector<int> alpha, int n, int degree) { */
/*     int toReturn; */ 
/*     for (int j = 0; j < n; ++j) { */
/*         toReturn += alpha[j]*pow(degree + 1, n -j - 1); */
/*     } */
/*     return toReturn; */
/* } */

/* // Main method */
/* int main(int argc, char *argv[]) { */

/*     // Degree of polynomials approximation. */
/*     int degree = 6; */

/*     // Macro time step */
/*     double Dt = .1; */

/*     // Number of fast variables */
/*     int nf = 3; */

/*     // Eigenvalues */
/*     vector<double> lambdas = {1.,2.,4.}; */
/*     vector<double> qs      = {2.,4.,2.}; */
/*     vector<double> sigmas(nf, 0.); */
/*     for (int i = 0; i < nf; ++i) { */
/*         sigmas[i] = sqrt(0.5*qs[i]*qs[i]/lambdas[i]); */
/*     } */

/*     // Number of polynomials in the basis */
/*     int nBasis = bin(degree + nf, nf); */

/*     // Relation linear index - multiindex */
/*     vector< vector<int> > ind2mult(nBasis, vector<int>(nf,0)); */
/*     vector<int> mult2ind(pow(degree + 1, nf), -1); */
/*     vector<int> currentMult(nf,0); */
/*     for (int i = 0; i < nBasis; ++i) { */
/*         for (int j = 0; j < nf; ++j) { */
/*             ind2mult[i][j] = currentMult[j]; */
/*         } */
/*         mult2ind[canonicalInd(currentMult, nf, degree)] = i; */
/*         int sum = 0; */
/*         for (int j = 0; j < nf; ++j) { */
/*             sum += currentMult[j]; */
/*         } */
/*         if (sum < degree) { */
/*             currentMult[nf-1] += 1; */
/*         } else { */
/*             int auxIndex = nf - 1; */
/*             while (currentMult[auxIndex] == 0) { */
/*                 auxIndex --; */
/*             } */
/*             currentMult[auxIndex] = 0; */
/*             currentMult[auxIndex-1] ++; */
/*         } */ 
/*     } */

/*     // Final time */
/*     double T = 1; */


/*     // Vector of times */
/*     int sizet = (int) (T/Dt) + 1; */
/*     vector<double> t(sizet,0.); */
/*     for (int i = 0; i < sizet; i++) { */
/*         t[i] = i*Dt; */
/*     } */

/*     // Parameters for random numbers */
/*     int seed = time(NULL); */
/*     default_random_engine generator; generator.seed(seed); */
/*     normal_distribution<double> distribution(0.0,1.0); */

/*     // Slow variable */
/*     vector<double> x(sizet,0.2); */

/*     for (int i = 0; i < sizet; ++i) { */
/*         cout << "Time of simulation: " << t[i] << endl; */

/*         // Expansion of right-hand side of the Poisson equation */
/*         vector<double> coefficients(nBasis, 0.); */
/*         vector<double> coefficients_dx(nBasis, 0.); */
/*         vector< vector<double> > coefficients_h(nBasis, vector<double>(nf,0.)); */
/*         for (int j = 0; j < nBasis; ++j) { */
/*             vector<int> multIndex = ind2mult[j]; */
/*             int N_mc = 100000; */
/*             double sum = 0.; */

/*             // Monte-Carlo to compute the coefficients */
/*             for (int k = 0; k < N_mc; ++k) { */
/*                 vector<double> randn(nf, 0.); */
/*                 double h_eval = 1.; */
/*                 for (int l = 0; l < nf; ++l) { */
/*                     randn[l] = distribution(generator); */
/*                     h_eval *= hermite(multIndex[l],randn[l],sigmas[l]); */
/*                 } */
/*                 coefficients[j] += h_eval*slow_drift(x[i],randn); */
/*                 coefficients_dx[j] += h_eval*slow_drift_dx(x[i],randn); */
/*                 vector<double> fast_drift_aux = fast_drift_h(x[i],randn); */
/*                 for (int l = 0; l < nf; ++l) { */
/*                     coefficients_h[j][l] += h_eval*fast_drift_aux[l]; */
/*                 } */
/*             } */
/*             coefficients[j] = coefficients[j]/N_mc; */
/*             coefficients_dx[j] = coefficients_dx[j]/N_mc; */
/*             for (int l = 0; l < nf; ++l) { */
/*                 coefficients_h[j][l] = coefficients_h[j][l]/N_mc; */
/*             } */
/*         } */

/*         // Solution of the Poisson equation */
/*         vector<double> solution(nBasis, 0.); */
/*         vector<double> solution_dx(nBasis, 0.); */
/*         vector< vector<double> > solution_dy(nBasis, vector<double>(nf,0.)); */
/*         for (int j = 0; j < nBasis; ++j) { */
/*             double eig = 0.; */
/*             for (int k = 0; k < nf; ++k) { */
/*                 eig += ind2mult[j][k]*lambdas[k]; */
/*             } */
/*             solution[j] = coefficients[j]/eig; */
/*             solution_dx[j] = coefficients_dx[j]/eig; */
/*             vector<int> thisMult = ind2mult[j]; */
/*             int sum = 0; */
/*             for (int l = 0; l < nf; ++l) { */
/*                 sum += thisMult[l]; */
/*             } */
/*             if (sum < degree) { */
/*                 for (int l = 0; l < nf; ++l) { */
/*                     vector<int> newMult(nf, 0); */
/*                     for (int m = 0; m < nf; ++m) { */
/*                         newMult[m] = thisMult[m]; */
/*                     } */
/*                     newMult[l] = newMult[l] + 1; */
/*                     int newInd = mult2ind[canonicalInd(newMult, nf, degree)]; */
/*                     solution_dy[j][l] = solution[newInd]*sqrt(newMult[l])/sigmas[l]; */
/*                 } */
/*             } */
/*         } */


/*         // Calculation of the coefficients of the simplified equation */
/*         double F1 = 0., F2 = 0., A0 = 0.; */
/*         for (int j = 0; j < nBasis; ++j) { */

/*             F1 += solution_dx[j]*coefficients[j]; */

/*             for (int k = 0; k < nf; ++k) { */
/*                 F2 += solution_dy[j][k]*coefficients_h[j][k]; */
/*             } */

/*             A0 += solution[j]*coefficients[j]; */
/*         } */
/*         cout << endl; */
/*     } */
/* } */
