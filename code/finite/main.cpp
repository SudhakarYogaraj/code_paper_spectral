#include "header.h"

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
    cout.precision(2);
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
        double error_hmm = 0.; 
        double error_spectral = 0.; 
        
        for (int i = 0; i < sizet - 1; i++) {

            // Numerical estimation of the drift and diffusion coefficients
            // at time step i (initialized at 0).
            vector<double> fi_hmm(problem.d, 0.);
            vector<double> fi_spectral(problem.d, 0.);
            vector< vector<double> > hi_hmm(problem.d, vector<double>(problem.d,0.));
            vector< vector<double> > hi_spectral(problem.d, vector<double>(problem.d,0.));

            // Initial value for the fast process at each macro time-step
            vector<double> yInit(2*problem.nf, 0.); 

            int seed = (int) abs(1000*distribution(generator));

            // Solution of the problem using the HMM method
            solve_hmm(problem, solver, xt[i], yInit, fi_hmm, hi_hmm, seed, t[i]);
            solve_spectral(problem, solver, xt[i], fi_spectral, hi_spectral, seed, t[i]);


            // Exact drift and diffusion coefficients
            vector<double> exact_drif = problem.soldrif(xt[i]);
            vector< vector<double> > exact_diff = problem.soldiff(xt[i]);

            vector<double> Ddrif(problem.d, 0.);
            vector< vector<double> > Ddiff(problem.d, vector<double>(problem.d,0.));

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    Ddiff[i1][i2] = hi_hmm[i1][i2] - exact_diff[i1][i2];
                }
                Ddrif[i1] = fi_hmm[i1] - exact_drif[i1];
            }

            error_hmm += 1./sizet*(normVec(Ddrif) + normMat(Ddiff));
            double errorDrift_hmm = normVec(Ddrif)/normVec(exact_drif);
            double errorDiff_hmm  = normMat(Ddiff)/normMat(exact_diff);

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    Ddiff[i1][i2] = hi_spectral[i1][i2] - exact_diff[i1][i2];
                }
                Ddrif[i1] = fi_spectral[i1] - exact_drif[i1];
            }

            error_spectral += 1./sizet*(normVec(Ddrif) + normMat(Ddiff));
            double errorDrift_spectral = normVec(Ddrif)/normVec(exact_drif);
            double errorDiff_spectral  = normMat(Ddiff)/normMat(exact_diff);

            // Output to terminal
            cout << "o-----------------------------------------------------------------------------------------------------o" << endl;
            cout << "|----------------- Iteration " << setw(3) <<  i+1 << "/" << sizet-1 << ". Time: " << t[i] << ". Precision parameter: " << solver.p <<". -----------------|" << endl;
            cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "|" << left << setw(50) <<  " 1) The HMM method:" << "|";
            cout << setw(50) <<  " 2) The Hermite spectral method:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "|" << setw(50) << "   -Drift coefficient:" << "|"; 
            cout << setw(50) <<  "   -Drift coefficient:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            print2Vecs(fi_hmm, fi_spectral);
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "|" << setw(50) << "   -Diffusion coefficient:" << "|"; 
            cout << setw(50) <<  "   -Diffusion coefficient:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            print2Mats(hi_hmm, hi_spectral);
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "|" << setw(50) <<  "   -Value of the slow variable x:" << "|";
            cout << setw(50) <<  "   -Value of the slow variable x:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            print2Vecs(xt[i],xt[i]);
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;

            cout << "|" << left << setw(50) <<  " 3) Error for the HMM method:" << "|";
            cout << setw(50) <<  " 4) Error for the Hermite spectral method:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "|" << setw(50) << "   -Error for drift coefficient:" << "|"; 
            cout << setw(50) <<  "   -Error for drift coefficient:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            cout << "|" << "    " << setw(46) <<  errorDrift_hmm  << "|" << "    " <<  setw(46) << errorDrift_spectral <<  "|" <<  endl; 
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            cout << "|" << setw(50) << "   -Error for diffusion coefficient:" << "|"; 
            cout << setw(50) <<  "   -Error for diffusion coefficient:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            cout << "|" << "    " << setw(46) <<  errorDiff_hmm << "|" << "    " << setw(46) << errorDiff_spectral <<  "|" <<  endl; 
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            cout << "|" << setw(50) << "   -Total error up to current iteration:" << "|"; 
            cout << setw(50) <<  "   -Total error up to current iteration:" << "|" << endl;
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 
            cout << "|" << "    " << setw(46) <<  error_hmm << "|" << "    " << setw(46) << error_spectral <<  "|" <<  endl; 
            cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl; 

            cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;

            cout << "|" << setw(101) <<  "5) Exact solution" << "|" << endl; 
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            cout << "|" << setw(101) << "   -Drift coefficient" << "|" << endl;
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            printVec(problem.soldrif(xt[i]));
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            cout << "|" << setw(101) << "   -Diffusion coefficient" << "|" << endl;
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            printMat(problem.soldiff(xt[i]));
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            cout << "|" << setw(101) << "   -Value of the slow variable x" << "|" << endl;
            cout << "|" << setw(101) <<  " " << "|" << endl; 
            printVec(x_exact[i]);
            cout << "|" << setw(101) <<  " " << "|" << endl; 

            cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;
            cout << endl << endl;

            // Computation of the exact coefficients based on the exact solution
            exact_drif = problem.soldrif(x_exact[i]);
            exact_diff = problem.soldiff(x_exact[i]);

            x_exact[i+1] = x_exact[i];
            xt[i+1] = xt[i];

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    x_exact[i+1][i1] += exact_diff[i1][i2]*sqrt(solver.macro_dt)*dWs[i][i2]; 
                    xt[i+1][i1] += hi_hmm[i1][i2]*sqrt(solver.macro_dt)*dWs[i][i2];
                }
                x_exact[i+1][i1] += solver.macro_dt*exact_drif[i1];
                xt[i+1][i1] += solver.macro_dt*fi_hmm[i1];
            }
        }
        writeToFile("time.dat",t); int p_aux = (int) (10*solver.p + 0.0001);
        writeMatToFile("xt" + to_string(p_aux) + ".dat", xt);
        writeMatToFile("x_exact.dat", x_exact);

        /* cout << "Total error: " << error << endl; */ 

        // log2 of the error, used to produce a plot
        /* errors[j] = log2(error); */
    }

    /* for (unsigned int i = 0; i < p_values.size(); i++) { */
    /*     cout << "Error for p = " << p_values[i] << ": " << errors[i] << endl; */ 
    /* } */
    /* cout << endl << endl; */

    // writeToFile("errors.dat", errors);
    // writeToFile("p_values.dat", p_values); 
}
