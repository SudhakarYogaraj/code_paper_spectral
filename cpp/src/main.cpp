#include "main.hpp"

/* TODO: Pass to general LU solver (urbain, Sat 23 May 2015 20:12:16 BST) */
/* TODO: Test gradient structures in several dimensions (urbain, Thu 21 May 2015 12:12:14 BST) */
/* TODO: automatic problem generation for hmm (urbain, Thu 21 May 2015 18:23:04 BST) */

using namespace std;

int main(int argc, char* argv[]) {

    int output1 = 0;

    // Initialization of the problem
    Problem problem;
    problem.init();

    // Values of the precision parameter
    vector<double> p_values = {6.};

    // Vector of the log of the error
    vector<double> errors_hmm(p_values.size(), 0.);
    vector<double> errors_spectral(p_values.size(), 0.);

    // Random numbers generator
    default_random_engine generator; generator.seed(time(NULL));
    normal_distribution<double> distribution(0.0,1.0);

    // Precision of the cout command
    cout.precision(6);
    cout << scientific;

    // Macro time-step
    double Dt = .01;

    // Initialization of variables and vectors
    // Size of the time vector for the macro scheme.
    unsigned int sizet = (int) (problem.t_end/Dt) + 1;

    // Vector of times of the macro-simulation
    vector<double> t(sizet,0.);
    for (unsigned int i = 0; i < sizet; i++) {
        t[i] = i*Dt;
    }

    // Vector of random variables used to simulate the brownian
    // motion for the evolution of the slow variable.
    vector< vector<double> > dWs(sizet,vector<double>(problem.ns, 0.));
    for (unsigned int i = 0; i < sizet; i++) {
        for (int j = 0; j < problem.ns ; j++) {
            dWs[i][j] = distribution(generator);
        }
    }

    // OUTPUT 1: GRAPH TIME - PRECISION FOR SPECTRAL METHOD

    if (output1) {

        vector<int> estimator_degrees = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
        vector<double> estimator_time(estimator_degrees.size());
        vector<double> estimator_error(estimator_degrees.size());

        ofstream out_time("out/spectral_time");
        ofstream out_errs("out/spectral_error");

        cout << "* Generating output for graph precision-time (spectral)" << endl;
        for (unsigned int i = 0; i < estimator_degrees.size(); ++i) {

            // Create new solver
            Solver_spectral solver_spectral = Solver_spectral(estimator_degrees[i], 100, problem.nf);

            // Update of the statistics of the invariant measure
            problem.update_stats(problem.x0);

            // Measure time of execution
            tic(); SDE_coeffs ci = solver_spectral.estimator(problem, problem.x0, 0.);
            estimator_time[i] = toc();

            // Erorr
            vector<double> Ddrif = (ci.drif - problem.soldrif(problem.x0));
            vector< vector<double> > Ddiff = (ci.diff - problem.soldiff(problem.x0));
            estimator_error[i] = fabs(Ddrif);// + fabs(Ddiff);

            // Summary of iteration
            cout << "    Iteration " << i << ". Time: " << estimator_time[i] << ". Error: " << estimator_error[i] << endl;

            // Write to file
            out_time << estimator_time[i] << endl;
            out_errs << estimator_error[i] << endl;
        }

        writeToFile("spectral_time", estimator_time);
        writeToFile("spectral_error", estimator_error);
    }

    for (unsigned int j = 0; j < p_values.size(); ++j) {

        Solver_hmm solver_hmm = Solver_hmm(p_values[j], 1);
        Solver_spectral solver_spectral = Solver_spectral(10, 100, problem.nf);

        // Approximate and exact solutions
        vector< vector<double> > xt_hmm(sizet,vector<double>(problem.ns,0.));
        vector< vector<double> > xt_spectral(sizet,vector<double>(problem.ns,0.));
        vector< vector<double> > x_exact(sizet,vector<double>(problem.ns,0.));

        // Initial condition
        xt_hmm[0] = problem.x0;
        xt_spectral[0] = problem.x0;
        x_exact[0] = problem.x0;

        // Error due to the estimation of the coefficients
        double error_hmm = 0.;
        double error_spectral = 0.;

        for (unsigned int i = 0; i < sizet - 1; i++) {

            struct SDE_coeffs c_hmm;
            struct SDE_coeffs c_spectral;
            vector<SDE_coeffs> vec_spectral;

            // Initial value for the fast process at each macro time-step
            vector<double> yInit(2*problem.nf, 0.);

            int seed = (int) abs(1000*distribution(generator));

            // Update of the statistics of the invariant measure
            problem.update_stats(x_exact[i]);

            // Solution of the problem using the HMM method
            solver_hmm.estimator(problem, xt_hmm[i], yInit, c_hmm, seed, t[i]);
            c_spectral = solver_spectral.estimator(problem, xt_spectral[i], t[i]);

            // Exact drift and diffusion coefficients
            vector<double> Ddrif = (c_hmm.drif - problem.soldrif(xt_hmm[i]));
            vector< vector<double> > Ddiff = (c_hmm.diff - problem.soldiff(xt_hmm[i]));
            error_hmm += 1./sizet*(fabs(Ddrif) + fabs(Ddiff));
            double errorDrift_hmm = fabs(Ddrif);
            double errorDiff_hmm  = fabs(Ddiff);

            /* cout << "--Drift--" << endl; */
            /* vector<double> exact_drift = problem.soldrif(xt_spectral[i]); */
            /* for (unsigned int m = 0; m < vec_spectral.size(); ++m) { */
            /*     Ddrif = vec_spectral[m].drif - exact_drift; */
            /*     double errorDrift_spectral = fabs(Ddrif); */
            /*     cout << errorDrift_spectral << endl; */
            /* } */
            /* cout << "--Diffusion--" << endl; */
            /* vector< vector<double> > exact_diffu = problem.soldiff(xt_spectral[i]); */
            /* for (unsigned int m = 0; m < vec_spectral.size(); ++m) { */
            /*     Ddiff = vec_spectral[m].diff - exact_diffu; */
            /*     double errorDiff_spectral  = fabs(Ddiff); */
            /*     cout << errorDiff_spectral << endl; */
            /* } */
            /* cout << "----" << endl; */

            Ddrif = c_spectral.drif - problem.soldrif(xt_spectral[i]);
            Ddiff = c_spectral.diff - problem.soldiff(xt_spectral[i]);
            error_spectral += 1./sizet*(fabs(Ddrif) + fabs(Ddiff));
            double errorDrift_spectral = fabs(Ddrif);
            double errorDiff_spectral  = fabs(Ddiff);

            // Computation of the exact coefficients based on the exact solution
            vector<double> exact_drif = problem.soldrif(x_exact[i]);
            vector< vector<double> > exact_diff = problem.soldiff(x_exact[i]);

            x_exact[i+1] = x_exact[i];
            xt_hmm[i+1] = xt_hmm[i];
            xt_spectral[i+1] = xt_spectral[i];

            for (int i1 = 0; i1 < problem.ns; i1++) {
                for (int i2 = 0; i2 < problem.ns; i2++) {
                    x_exact[i+1][i1] += exact_diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    xt_hmm[i+1][i1] += c_hmm.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                    xt_spectral[i+1][i1] += c_spectral.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                }
                x_exact[i+1][i1] += Dt*exact_drif[i1];
                xt_hmm[i+1][i1] += Dt*c_hmm.drif[i1];
                xt_spectral[i+1][i1] += Dt*c_spectral.drif[i1];
            }

            // Output to terminal
            if (SUMMARY) {
                cout << "o-----------------------------------------------------------------------------------------------------o" << endl;
                cout << "|------------- Iteration " << setw(3) <<  i+1 << "/" << sizet-1 << ". Time: " << t[i] << ". Precision parameter: " << solver_hmm.p <<". -------------|" << endl;
                cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << left << setw(50) <<  " 1) The HMM method:" << "|";
                cout << setw(50) <<  " 2) The Hermite spectral method:" << "|" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << setw(50) << "   -Drift coefficient:" << "|";
                cout << setw(50) <<  "   -Drift coefficient:" << "|" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                print2Vecs(c_hmm.drif, c_spectral.drif);
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << setw(50) << "   -Diffusion coefficient:" << "|";
                cout << setw(50) <<  "   -Diffusion coefficient:" << "|" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                print2Mats(c_hmm.diff, c_spectral.diff);
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
                cout << "|" << setw(50) << "   -Error for the slow process:" << "|";
                cout << setw(50) <<  "   -Error for the slow process:" << "|" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << "    " << setw(46) <<  fabs(x_exact[i] - xt_hmm[i]) << "|" << "    " << setw(46) << fabs(x_exact[i] - xt_spectral[i]) <<  "|" <<  endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << setw(50) << "   -Total error up to current iteration:" << "|";
                cout << setw(50) <<  "   -Total error up to current iteration:" << "|" << endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "|" << "    " << setw(46) <<  error_hmm << "|" << "    " << setw(46) << error_spectral <<  "|" <<  endl;
                cout << "|" << setw(50) <<  " " << "|" << setw(50) << " " <<  "|" <<  endl;
                cout << "o--------------------------------------------------o--------------------------------------------------o" << endl;
                cout << "|" << setw(101) <<  " 5) Exact solution" << "|" << endl;
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
                cout << "o-----------------------------------------------------------------------------------------------------o" << endl;
                cout << endl << endl;
            }
        }
        writeToFile("time.dat",t); int p_aux = (int) (10*solver_hmm.p + 0.0001);
        writeMatToFile("xt_hmm" + to_string(p_aux) + ".dat", xt_hmm);
        writeMatToFile("xt_spectral" + to_string(p_aux) + ".dat", xt_spectral);
        writeMatToFile("x_exact.dat", x_exact);
        errors_hmm[j] = log2(error_hmm);
        errors_spectral[j] = log2(error_spectral);
    }

    cout << endl << "Summary for the hmm method" << endl;
    for (unsigned int i = 0; i < p_values.size(); i++) {
        cout << "Error for p = " << p_values[i] << ": " << errors_hmm[i] << endl;
    }

    cout << endl << "Summary for the spectral method" << endl;
    for (unsigned int i = 0; i < p_values.size(); i++) {
        cout << "Error for p = " << p_values[i] << ": " << errors_spectral[i] << endl;
    }
    writeToFile("errors_hmm.dat", errors_hmm);
    writeToFile("errors_spectral.dat", errors_spectral);
    writeToFile("p_values.dat", p_values);
}
