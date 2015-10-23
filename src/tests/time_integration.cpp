// Time tracking
#include "io/tictoc.hpp"

// Structures
#include "structures.hpp"

// Problem and analyser
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"

// Solvers
#include "solvers/Solver_spectral.hpp"
#include "solvers/Solver_hmm.hpp"
#include "solvers/Solver_exact.hpp"

// Own header
#include "tests/time_integration.hpp"

using namespace std;

namespace tests {

    void integrate(Problem *problem, Solver *se, Solver_spectral *ss, Solver_hmm *sh) {

        // Precision of the cout command
        cout.precision(6);
        cout << scientific;

        // Macro time-step
        double Dt = .01;
        int sizet = (int) (problem->t_end/Dt) + 1;

        // Vector of times of the macro-simulation
        vector<double> t(sizet,0.);
        for (unsigned int i = 0; i < sizet; i++) {
            t[i] = i*Dt;
        }

        // Random numbers generator
        default_random_engine generator; generator.seed(time(NULL));
        normal_distribution<double> distribution(0.0,1.0);

        // Vector of random variables used to simulate the brownian
        // motion for the evolution of the slow variable.
        vector< vector<double> > dWs(sizet,vector<double>(problem->ns, 0.));
        for (unsigned int i = 0; i < sizet; i++) {
            for (int j = 0; j < problem->ns ; j++) {
                dWs[i][j] = distribution(generator);
            }
        }

        // Approximate and exact solutions
        vector< vector<double> > x_exact(sizet,vector<double>(problem->ns,0.));
        vector< vector<double> > x_spectral(sizet,vector<double>(problem->ns,0.));
        vector< vector<double> > x_hmm(sizet,vector<double>(problem->ns,0.));

        // Initial condition
        x_exact[0] = problem->x0;
        x_spectral[0] = problem->x0;
        x_hmm[0] = problem->x0;

        // Streams
        ofstream out_time("out/time_integration/time");
        ofstream out_exact("out/time_integration/exact");
        ofstream out_spectral("out/time_integration/spectral");
        ofstream out_hmm("out/time_integration/hmm");

        for (unsigned int i = 0; i < sizet - 1; i++) {

            struct SDE_coeffs c_exact;
            struct SDE_coeffs c_spectral;
            struct SDE_coeffs c_hmm;

            // Solution of the problem using the HMM method
            c_exact = se->estimator(x_exact[i], t[i]);
            c_spectral = ss->estimator(x_spectral[i], t[i]);
            c_hmm = sh->estimator(x_hmm[i], t[i]);

            // Initial value for new time step
            x_exact[i+1] = x_exact[i];
            x_hmm[i+1] = x_hmm[i];
            x_spectral[i+1] = x_spectral[i];

            for (int i1 = 0; i1 < problem->ns; i1++) {
                for (int i2 = 0; i2 < problem->ns; i2++) {
                    x_exact[i+1][i1] += c_exact.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    x_spectral[i+1][i1] += c_spectral.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    x_hmm[i+1][i1] += c_hmm.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                }
                x_exact[i+1][i1] += Dt*c_exact.drif[i1];
                x_spectral[i+1][i1] += Dt*c_spectral.drif[i1];
                x_hmm[i+1][i1] += Dt*c_hmm.drif[i1];
            }

            out_time << t[i] << endl;
            for (int j = 0; j < problem->ns; ++j) {
                out_exact << x_exact[i][j];
                out_spectral << x_spectral[i][j];
                out_hmm << x_hmm[i][j];
                if (j != problem->ns -1) {
                    out_exact << " ";
                    out_spectral << " ";
                    out_hmm << " ";
                }
                else {
                    out_exact << endl;
                    out_spectral << endl;
                    out_hmm << endl;
                }
            }
        }
    }
}
