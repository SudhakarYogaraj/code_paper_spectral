// Includes
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <time.h>

#define PI 3.141592653589793238463

// Namespace
using namespace std;

class Problem {
    public:
        
        // Final time
        double t_end;

        // Number of fast processes
        int nf;

        // Number of slow processes
        int d;

        // Initial condition for the slow process
        vector<double> x0;
        
        // Drift coefficient of the slow process
        vector<double> a(vector<double> x, vector<double> y); // f0 in Pavliotis-Stuart
        
        // Derivatives of a(,)
        vector< vector<double> > dax(vector<double> x, vector<double> y);
        vector< vector<double> > day(vector<double> x, vector<double> y);

        // FOR SPECTRAL
        // Non-leading order part of drift in the fast process
        vector<double>  fast_drift_h(vector<double> x, vector<double> y);
        vector<double> lambdas;
        vector<double> betas;

        // FOR HMM
        // Drift and diffusion terms of the fast system, in its transformed 
        // version. The first nf components correspond to the initial 
		// variables, whereas the last nf components correspend to the
		// auxiliary variables.
        vector<double> drif(vector<double> x, vector<double> y);
        vector<double> diff(vector<double> x, vector<double> y);

		// Drift and diffusion coefficients of the solution
		vector<double> soldrif(vector<double> x);
		vector< vector<double> > soldiff(vector<double> x);

		// Initialization of the problem
		void init();
};

class Solver {
    public: 

        // Macro and micro time-steps 
        double macro_dt; 
        double micro_dt;

        // Parameters of the estimator
        int n; 
        int nt; 
        int np;
        int M;

        // Order of accuracy of the micro-solver
        double l; 

        // Precision parameter of the solver
        double p;

		void set(double p, int M);
};

void solve_hmm(Problem &problem, Solver &solver, vector<double> xt, vector<double>& yInit, vector<double>& fi, vector< vector<double> >& hi, int seed, double t);
void solve_spectral(Problem &problem, Solver &solver, vector<double> xt, vector<double>& fi, vector< vector<double> >& hi, int seed, double t);

// Function to write vector to a file
void writeToFile(string s, vector<double> x);

// Kronecker delta
double delta(int a, int b);

// Symmetric part
vector< vector<double> > symmetric(vector< vector<double> > A);

// Cholesky factorization
vector< vector<double> > cholesky(vector< vector<double> > A);

// Printing vector
void print2Vecs (vector<double> x, vector<double> y);
void printVec (vector<double> x);

// Printing matrix
void print2Mats (vector< vector<double> > x, vector< vector<double> > y);
void printMat (vector< vector<double> > x);

// Norm of a vector
double normVec (vector<double> x);

// Norm of a matrix
double normMat (vector< vector<double> > x);

// Write matrix to a file
void writeMatToFile(string s, vector< vector<double> > x);

// Normalized Hermite polynomials
double hermite(int n, double x, double sigma);    
double hermiteM(vector<int> multIndex, vector<double> x, vector<double> sigmas);

// Binomial coefficients
int bin(int n, int k);

// Function that associates a multiIndex with a "canonical linear index"
int canonicalInd(vector<int> alpha, int n, int degree);
