#include "tests/time_integration.hpp"

using namespace std;

namespace tests {

    void integrate(Problem *problem, Solver *solver, int seed, vec& time, mat& solution) {

        // Macro time-step
        double Dt = .01;
        int nSteps = 100;

        // Create vector of times
        vec t(nSteps + 1,0.);
        for (int i = 0; i < nSteps + 1; i++)
        {
            t[i] = i*Dt;
        }

        // Create brownian motion
        mat dWs(nSteps + 1,vec(problem->ns, 0.));
        default_random_engine generator; generator.seed(seed);
        normal_distribution<double> distribution(0.0,1.0);

        for (int i = 0; i < nSteps; i++)
        {
            for (int j = 0; j < problem->ns ; j++)
            {
                dWs[i][j] = distribution(generator);
            }
        }

        // Calculate approximate solution
        mat x(nSteps + 1, vec(problem->ns,0.));
        x[0] = problem->x0;

        for (int i = 0; i < nSteps; i++) {

            struct SDE_coeffs c = solver->estimator(x[i], t[i]);
            x[i+1] = x[i];

            for (int i1 = 0; i1 < problem->ns; i1++) {
                for (int i2 = 0; i2 < problem->ns; i2++) {
                    x[i+1][i1] += c.diff[i1][i2]*sqrt(Dt)*dWs[i][i2];
                }
                x[i+1][i1] += Dt*c.drif[i1];
            }
        }

        time = t;
        solution = x;
    }
}

int main(int argc, char* argv[]) {

    // Initialization of the problem and helper analyser
    Problem problem; problem.init();
    Analyser analyser(&problem);

    // Degrees for spectral methmd
    int degree_min = 5;
    int degree_max = 30;

    vec degrees(degree_max - degree_min + 1);
    for (unsigned int i = 0; i < degrees.size(); ++i)
    {
        degrees[i] = degree_min + i;
    }

    // Define variables
    vec time;
    mat sol_exact;
    cube sol_spectral(degrees.size());

    // Integrate in time using exact and spectral solvers
    Solver_exact solver_exact(&problem, &analyser);
    tests::integrate(&problem, &solver_exact, 0, time, sol_exact);

    for (unsigned int i = 0; i < degrees.size(); ++i)
    {
        config_spectral conf_spectral; {
            conf_spectral.n_nodes = 100;
            conf_spectral.degree = degrees[i];
            conf_spectral.scaling = vec(problem.nf, problem.sigma);
        }

        Solver_spectral solver_spectral = Solver_spectral(&problem, &analyser, &conf_spectral);
        tests::integrate(&problem, &solver_spectral, 0, time, sol_spectral[i]);
    }

    // Calculate errors
    cube error(degrees.size());
    mat error_abs(degrees.size(), vec(time.size()));
    vec error_sup(degrees.size(), 0.);

    for (unsigned int i = 0; i < degrees.size(); i++)
    {
        error[i] = sol_spectral[i] - sol_exact;
        for (unsigned int j = 0; j < time.size(); j++)
        {
            error_abs[i][j] = fabs(error[i][j]);
            if (error_abs[i][j] > error_sup[i])
            {
                error_sup[i] = error_abs[i][j];
            }
        }
    }

    ofstream out_degree("degree");
    ofstream out_error("error");
    for (unsigned int i = 0; i < degrees.size(); i++)
    {
        out_degree << degrees[i] << endl;
        out_error << error_sup[i] << endl;
    }
}
