#include "header.h"

// Main function
int main(int argc, char* argv[])
{
    // Initialization of the problem
    Problem problem;
    problem.init();

    // Values of the precision parameter
    vector<double> p_values = {4.};

    // Vector of the log of the error
    vector<double> errors_hmm(p_values.size(), 0.);
    vector<double> errors_spectral(p_values.size(), 0.);

    // Number of replicas of the fast process
    int M = 1;

    // Random variable for generator
    int seed = time(NULL);

    // Initialization of the solver
    Solver solver;
    Solver_spectral solver_spectral;

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
        t[i] = i*Dt;
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
        vector< vector<double> > xt_hmm(sizet,vector<double>(problem.d,0.));
        vector< vector<double> > xt_spectral(sizet,vector<double>(problem.d,0.));
        vector< vector<double> > x_exact(sizet,vector<double>(problem.d,0.));

        // Initial condition
        xt_hmm[0] = problem.x0;
        xt_spectral[0] = problem.x0;
        x_exact[0] = problem.x0;

        // Setting to solver to match the precision parameter p_value[i];
        solver.set(p_values[j],M);
        solver_spectral.set(p_values[j], problem.nf);

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
            solve_hmm(problem, solver, xt_hmm[i], yInit, fi_hmm, hi_hmm, seed, t[i]);
            solve_spectral(problem, solver_spectral, xt_spectral[i], fi_spectral, hi_spectral, t[i]);

            // Exact drift and diffusion coefficients
            vector<double> exact_drif = problem.soldrif(xt_hmm[i]);
            vector< vector<double> > exact_diff = problem.soldiff(xt_hmm[i]);

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

            exact_drif = problem.soldrif(xt_spectral[i]);
            exact_diff = problem.soldiff(xt_spectral[i]);

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    Ddiff[i1][i2] = hi_spectral[i1][i2] - exact_diff[i1][i2];
                }
                Ddrif[i1] = fi_spectral[i1] - exact_drif[i1];
            }

            error_spectral += 1./sizet*(normVec(Ddrif) + normMat(Ddiff));
            double errorDrift_spectral = normVec(Ddrif)/normVec(exact_drif);
            double errorDiff_spectral  = normMat(Ddiff)/normMat(exact_diff);

            // Computation of the exact coefficients based on the exact solution
            exact_drif = problem.soldrif(x_exact[i]);
            exact_diff = problem.soldiff(x_exact[i]);

            x_exact[i+1] = x_exact[i];
            xt_hmm[i+1] = xt_hmm[i];
            xt_spectral[i+1] = xt_hmm[i];

            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    x_exact[i+1][i1] += exact_diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    xt_hmm[i+1][i1] += hi_hmm[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    xt_spectral[i+1][i1] += hi_spectral[i1][i2]*sqrt(Dt)*dWs[i][i2];
                }
                x_exact[i+1][i1] += Dt*exact_drif[i1];
                xt_hmm[i+1][i1] += Dt*fi_hmm[i1];
                xt_spectral[i+1][i1] += Dt*fi_spectral[i1];
            }

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
            print2Vecs(xt_hmm[i],xt_spectral[i]);
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
            printVec(exact_drif);
            cout << "|" << setw(101) <<  " " << "|" << endl;
            cout << "|" << setw(101) << "   -Diffusion coefficient" << "|" << endl;
            cout << "|" << setw(101) <<  " " << "|" << endl;
            printMat(exact_diff);
            cout << "|" << setw(101) <<  " " << "|" << endl;
            cout << "|" << setw(101) << "   -Value of the slow variable x" << "|" << endl;
            cout << "|" << setw(101) <<  " " << "|" << endl;
            printVec(x_exact[i]);
            cout << "|" << setw(101) <<  " " << "|" << endl;
            cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;
            cout << endl << endl;
        }
        writeToFile("time.dat",t); int p_aux = (int) (10*solver.p + 0.0001);
        writeMatToFile("xt_hmm" + to_string(p_aux) + ".dat", xt_hmm);
        writeMatToFile("xt_spectral" + to_string(p_aux) + ".dat", xt_spectral);
        writeMatToFile("x_exact.dat", x_exact);
        errors_hmm[j] = log2(error_hmm);
        errors_spectral[j] = log2(error_spectral);
    }

    cout << endl << "Summary for the hmm method" << endl;
    for (int i = 0; i < p_values.size(); i++) {
        cout << "Error for p = " << p_values[i] << ": " << errors_hmm[i] << endl;
    }

    cout << endl << "Summary for the spectral method" << endl;
    for (int i = 0; i < p_values.size(); i++) {
        cout << "Error for p = " << p_values[i] << ": " << errors_spectral[i] << endl;
    }
    writeToFile("errors_hmm.dat", errors_hmm);
    writeToFile("errors_spectral.dat", errors_spectral);
    writeToFile("p_values.dat", p_values);
}
