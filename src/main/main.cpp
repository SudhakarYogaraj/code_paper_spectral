#include "main/main.hpp"

using namespace std;

namespace tests {

    void error_hmm(vector<double> x, Problem *problem, Analyser *analyser) {

        // Minimal and maximal values of the precision parameter
        int p_min = 4;
        int p_max = 6;

        vector<int> p_values(p_max - p_min + 1);
        for (int i = 0; i < p_values.size(); ++i)
            p_values[i] = p_min + i;

        // Update of the statistics of the invariant measure
        analyser->update_stats(x);

        // Computation of the exact solution
        Solver_exact solver_exact(problem, analyser);
        SDE_coeffs c_exact = solver_exact.estimator(x, 0.);
        vector<double> exact_drift = c_exact.drif;
        vector< vector<double> > exact_diff = c_exact.diff;

        vector<double> estimator_time(p_values.size());
        vector<double> estimator_error(p_values.size());

        // Files to write to
        ofstream out_time("hmm_time");
        ofstream out_errs("hmm_error");

        for (unsigned int i = 0; i < p_values.size(); ++i) {

            // Create new solvers
            config_hmm config = Solver_hmm::sensible_conf(p_values[i], 1);
            Solver_hmm solver_hmm(problem, &config);

            // Measure time of execution
            tic(); SDE_coeffs c = solver_hmm.estimator(x, 0.); estimator_time[i] = toc();

            // Error
            vector<double> Ddrif = (c.drif - exact_drift);
            vector< vector<double> > Ddiff = (c.diff - exact_diff);

            cout << "fabs(exact_drift)" << fabs(exact_drift) << endl;
            cout << "fabs(exact_diff)" << fabs(exact_diff) << endl;
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

    void error_spectral(vector<double> x, Problem *problem, Analyser *analyser) {

        // Maximal degree of Hermite polynomials
        int degree_min = 5;
        int degree_max = 30;

        vector<int> degrees(degree_max - degree_min + 1);
        for (int i = 0; i < degrees.size(); ++i)
            degrees[i] = i + degree_min;

        // Update of the statistics of the invariant measure
        analyser->update_stats(x);
        // Computation of the exact solution
        Solver_exact solver_exact(problem, analyser);
        SDE_coeffs c_exact = solver_exact.estimator(x, 0.);
        vector<double> exact_drift = c_exact.drif;
        vector< vector<double> > exact_diff = c_exact.diff;

        vector<double> estimator_time(degrees.size());
        vector<double> estimator_error(degrees.size());

        // Files to write to
        ofstream out_degree("out/error_spectral/degree");
        ofstream out_time("out/error_spectral/time");
        ofstream out_errs("out/error_spectral/error");

        // Configuration for the spectral solver
        config_spectral conf_spectral; {
            conf_spectral.n_nodes = 100;
            conf_spectral.degree = degrees[degrees.size()-1];
            conf_spectral.scaling = vector<double> (problem->nf, 0.5);
        }

        // Create new solver
        Solver_spectral solver_spectral(problem, analyser, &conf_spectral);

        // Find full solution
        tic(); vector<SDE_coeffs> full_c = solver_spectral.full_estimator(x, 0.,degrees);

        for (unsigned int i = 0; i < degrees.size(); ++i) {

            // Error
            vector<double> Ddrif = (full_c[i].drif - exact_drift);
            vector< vector<double> > Ddiff = (full_c[i].diff - exact_diff);
            estimator_error[i] = fabs(Ddrif)/fabs(exact_drift) + fabs(Ddiff)/fabs(exact_diff);

            // Summary of iteration
            cout << "> Degree" << degrees[i] << ". Error: " << estimator_error[i] << endl;
        }

        cout << "Time for full solution: " << toc() << endl << endl;

        for (unsigned int i = 0; i < degrees.size(); ++i) {

            // Configuration for the spectral solver
            config_spectral conf_spectral; {
                conf_spectral.n_nodes = 100;
                conf_spectral.degree = degrees[i];
                conf_spectral.scaling = vector<double> (problem->nf, 0.5);
            }

            // Create new solvers
            Solver_spectral solver_spectral(problem, analyser, &conf_spectral);

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

        /* Do you notice a difference between the two solutions?
         * This is because the analyser updates the statistics each time.
         */
    }
}

int main(int argc, char* argv[]) {

    // Initialization of the problem and helper analyser
    Problem problem;
    problem.init();
    Analyser analyser(&problem);

    // Test error_spectral
    tests::error_hmm(problem.x0, &problem, &analyser);
    tests::error_spectral(problem.x0, &problem, &analyser);
}
