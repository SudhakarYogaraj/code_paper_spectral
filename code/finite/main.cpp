#include "header.h"

/* TODO: Integrate spectral method in the code (urbain, Mon 20 Apr 2015 14:40:06 BST) */

// Main function
int main(int argc, char* argv[])
{
    // Initialization of the problem
    Problem problem; 
    problem.init();

    // Values of the precision parameter
    vector<double> p_values = {5.};

    // Vector of the log of the error
    vector<double> errors(p_values.size(), 0.);

    // Number of replicas of the fast process
    int M = 1;

    // Random variable for generator
    int seed = time(NULL);

    // Initialization of the solver
    Solver solver;

    // Random numbers generator
    default_random_engine generator; generator.seed(seed);
    normal_distribution<double> distribution(0.0,1.0);
    
    // Precision of the cout command
    cout.precision(8);
    cout << scientific;

    // Macro time-step
    double Dt = .01;

    // Initialization of variables and vectors
    // Size of the time vector for the macro scheme.
    int sizet = (int) (problem.t_end/Dt) + 1;

    // Vector of times of the macro-simulation
    vector<double> t(sizet,0.);
    for (int i = 0; i < sizet; i++) {
        t[i] = i*solver.macro_dt;
    }

    // Vector of random variables used to simulate the brownian 
    // motion for the evolution of the slow variable.
    vector< vector<double> > dWs(sizet,vector<double>(problem.d, 0.));
    for (int i = 0; i < sizet; i++) {
        for (int j = 0; j < problem.d ; j++) {
            dWs[i][j] = distribution(generator);
        }
    }

    for (int j = 0; j < p_values.size(); ++j) {

        // Approximate and exact solutions
        vector< vector<double> > xt(sizet,vector<double>(problem.d,0.)); 
        vector< vector<double> > x_exact(sizet,vector<double>(problem.d,0.)); 

        // Initial condition
        xt[0] = problem.x0;
        x_exact[0] = problem.x0;

        // Setting to solver to match the precision parameter p_value[i];
        solver.set(p_values[j],M);

        // Error due to the estimation of the coefficients
        double error = 0.; 
        
        for (int i = 0; i < sizet - 1; i++) {
            cout << "<< i = " << i << ", p = " << solver.p <<" >>" << endl << endl;


            // Numerical estimation of the drift and diffusion coefficients
            // at time step i (initialized at 0).
            vector<double> fi(problem.d, 0.);
            vector< vector<double> > hi(problem.d, vector<double>(problem.d,0.));

            // Initial value for the fast process at each macro time-step
            vector<double> yInit(2*problem.nf, 0.); 

            int seed = (int) abs(1000*distribution(generator));

            // Solution of the problem using the HMM method
            solve_hmm(problem, solver, xt[i], yInit, fi, hi, seed, t[i]);


            // Output to terminal
            cout << "   Drift coefficient: " << endl; printVec(fi); cout << endl;
            cout << "   Diffusion coefficient: " << endl; printMat(hi); cout << endl;
            cout << "   Exact value of the drift coefficient: " << endl;
            printVec(problem.soldrif(xt[i])); cout << endl;
            cout << "   Exact value of the diffusion coefficient: " << endl;
            printMat(problem.soldiff(xt[i])); cout << endl;
            cout << "   Value of the slow variable x: " << endl;
            printVec(xt[i]); cout << endl;
            cout << "   Exact value of the slow variable: " << endl;
            printVec(x_exact[i]); cout << endl;

            // Exact drift and diffusion coefficients
            vector<double> exact_drif = problem.soldrif(xt[i]);
            vector< vector<double> > exact_diff = problem.soldiff(xt[i]);

            vector<double> Ddrif(problem.d, 0.);
            vector< vector<double> > Ddiff(problem.d, vector<double>(problem.d,0.));

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    Ddiff[i1][i2] = hi[i1][i2] - exact_diff[i1][i2];
                }
                Ddrif[i1] = fi[i1] - exact_drif[i1];
            }

            error += 1./sizet*(normVec(Ddrif) + normMat(Ddiff));
            double errorDrift = normVec(Ddrif)/normVec(exact_drif);
            double errorDiff  = normMat(Ddiff)/normMat(exact_diff);

            cout << endl << "   Error in the drift term: " << errorDrift << endl;
            cout << "   Error in the diffusion term: " << errorDiff << endl;
            cout << "   Total error up to the current iteration: " << error << endl << endl;

            // Computation of the exact coefficients based on the exact solution
            exact_drif = problem.soldrif(x_exact[i]);
            exact_diff = problem.soldiff(x_exact[i]);

            x_exact[i+1] = x_exact[i];
            xt[i+1] = xt[i];

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    x_exact[i+1][i1] += exact_diff[i1][i2]*sqrt(solver.macro_dt)*dWs[i][i2]; 
                    xt[i+1][i1] += hi[i1][i2]*sqrt(solver.macro_dt)*dWs[i][i2];
                }
                x_exact[i+1][i1] += solver.macro_dt*exact_drif[i1];
                xt[i+1][i1] += solver.macro_dt*fi[i1];
            }
        }
        writeToFile("time.dat",t); int p_aux = (int) (10*solver.p + 0.0001);
        writeMatToFile("xt" + to_string(p_aux) + ".dat", xt);
        writeMatToFile("x_exact.dat", x_exact);

        cout << "Total error: " << error << endl; 

        // log2 of the error, used to produce a plot
        errors[j] = log2(error);
    }

    for (unsigned int i = 0; i < p_values.size(); i++) {
        cout << "Error for p = " << p_values[i] << ": " << errors[i] << endl; 
    }
    cout << endl << endl;

    // writeToFile("errors.dat", errors);
    // writeToFile("p_values.dat", p_values); 
}
