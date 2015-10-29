// Time tracking
#include "io/tictoc.hpp"

// Structures
#include "structures.hpp"

// Solvers
#include "solvers/Solver_spectral.hpp"
#include "solvers/Solver_exact.hpp"

// Own header
#include "tests/error_spectral.hpp"

using namespace std;

namespace tests {

    void error_spectral(vector<double> x, Problem *problem, Analyser *analyser) {

        // Maximal degree of Hermite polynomials
        int degree_min = 2;
        int degree_max = 30;

        // Number of integration points
        int ni = 100;

        vector<int> degrees(degree_max - degree_min + 1);
        for (int i = 0; i < degrees.size(); ++i)
            degrees[i] = i + degree_min;

        // Update of the statistics of the invariant measure
        analyser->update_stats(x);

        // Computation of the exact solution
        Solver_exact solver_exact(problem, analyser);
        vector<double> exact_drift = solver_exact.soldrif(x);
        vector< vector<double> > exact_diff = solver_exact.soldiff(x);

        vector<double> estimator_time(degrees.size());
        vector<double> estimator_error(degrees.size());

        // Files to write to
        ofstream out_degree("out/error_spectral/degree");
        ofstream out_time("out/error_spectral/time");
        ofstream out_errs("out/error_spectral/error");

        for (unsigned int i = 0; i < degrees.size(); ++i) {

            // Create new solvers
            Solver_spectral solver_spectral(problem, analyser, degrees[i], ni);

            // Measure time of execution
            tic(); SDE_coeffs c = solver_spectral.estimator(x, 0.); estimator_time[i] = toc();

            // Error
            vector<double> Ddrif = (c.drif - exact_drift);
            vector< vector<double> > Ddiff = (c.diff - exact_diff);
            estimator_error[i] = fabs(Ddrif)/fabs(exact_drift) + fabs(Ddiff)/fabs(exact_diff);

            // Write to file
            out_degree << degrees[i] << endl;
            out_time << estimator_time[i] << endl;
            out_errs << estimator_error[i] << endl;

            // Summary of iteration
            cout << "> Degree" << degrees[i] << ". Time: " << estimator_time[i] << ". Error: " << estimator_error[i] << endl;
        }
    }
}
