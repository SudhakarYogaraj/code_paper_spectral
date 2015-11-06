// Time tracking
#include "io/tictoc.hpp"
#include "templates.hpp"

// Problem and analyser
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"

// Solvers
#include "solvers/Solver_hmm.hpp"
#include "solvers/Solver_exact.hpp"

// Own header
#include "tests/error_hmm.hpp"

using namespace std;

namespace tests {

    void error_hmm(vector<double> x, Problem *problem, Analyser *analyser) {

        // Minimal and maximal values of the precision parameter
        int p_min = 3;
        int p_max = 6;

        vector<int> p_values(p_max - p_min + 1);
        for (int i = 0; i < p_values.size(); ++i)
            p_values[i] = p_min + i;

        // Update of the statistics of the invariant measure
        analyser->update_stats(x);

        // Computation of the exact solution
        Solver_exact solver_exact(problem, analyser);
        vector<double> exact_drift = solver_exact.soldrif(x);
        vector< vector<double> > exact_diff = solver_exact.soldiff(x);

        vector<double> estimator_time(p_values.size());
        vector<double> estimator_error(p_values.size());

        // Files to write to
        ofstream out_time("out/hmm_time");
        ofstream out_errs("out/hmm_error");

        for (unsigned int i = 0; i < p_values.size(); ++i) {

            // Create new solvers
            config_hmm config = Solver_hmm::sensible_conf(p_values[i], 1);
            Solver_hmm solver_hmm(problem, &config);

            // Measure time of execution
            tic(); SDE_coeffs c = solver_hmm.estimator(x, 0.); estimator_time[i] = toc();

            // Error
            vector<double> Ddrif = (c.drif - exact_drift);
            vector< vector<double> > Ddiff = (c.diff - exact_diff);
            cout << "fabs(Ddrif)" << fabs(Ddrif) << endl;
            cout << "fabs(Ddiff)" << fabs(Ddiff) << endl;
            estimator_error[i] = fabs(Ddrif)/fabs(exact_drift) + fabs(Ddiff)/fabs(exact_diff);

            // Write to file
            out_time << estimator_time[i] << endl;
            out_errs << estimator_error[i] << endl;

            // Summary of iteration
            cout << "> Precision :" << p_values[i] << ". Time: " << estimator_time[i] << ". Error: " << estimator_error[i] << endl;
        }

    }
}
