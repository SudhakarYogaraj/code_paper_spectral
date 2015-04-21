#include "header.h"

// Function implementing the numerical method
double solve(int seed, Problem &problem, Solver &solver) { 
    
    // ----- Random numbers generator ----- //

    default_random_engine generator; generator.seed(seed);
    normal_distribution<double> distribution(0.0,1.0);

    // ----- Precision of the cout command ----- //
    cout.precision(2);
    cout << scientific;
    
    // ----- Implementation of the HMM ----- // 

    // Initialization of variables and vectors
    // Size of the time vector for the macro scheme.
    int sizet = (int) (problem.t_end/solver.macro_dt) + 1;
        
    // Initialization of the arrays that will used to store 
    // the solution of the problem.
    vector<double> t(sizet,0.);
    vector< vector<double> > xt(sizet,vector<double>(problem.d,0.)); 
            
    // Vector that will contain the exact solution
    vector< vector<double> > x_exact(sizet,vector<double>(problem.d,0.)); 
    
    // Initial condition
    xt[0] = problem.x0;
    x_exact[0] = problem.x0;

    // Vector of random variables used to simulate the brownian 
    // motion for the evolution of the slow variable.
    vector< vector<double> > dWs(sizet,vector<double>(problem.d, 0.));
    for (int i = 0; i < sizet; i++) {
        for (int j = 0; j < problem.d ; j++) {
            dWs[i][j] = distribution(generator);
        }
    }

    // The brownian motion for the slow process is taken to be the same
    // for every value of the precision parameter p, so that a proper 
    // comparison can be made. 
    generator.seed(time(NULL));

    // Construction of the vector containing the times at which
    // the slow variable is calculated
    for (int i = 0; i < sizet; i++) {
        t[i] = i*solver.macro_dt;
    }

    // Construction of the array that will contain the solution for the 
    // fast process at each macro time-step. 
    vector< vector<double> > yAux(solver.nt + solver.n + solver.np, \
            vector<double>(2*problem.nf,0.)); 

    // Initial value for the fast process at each macro time-step
    vector<double> yInit(2*problem.nf, 0.); 
   
    // Error due to the estimation of the coefficients
    double error = 0.; 

    // Implementation of the numerical method. 
    for (int i = 0; i < sizet - 1; i++) {
        
        cout << "<< i = " << i << ", p = " << solver.p <<" >>" << endl << endl;

        // Numerical estimation of the drift and diffusion coefficients
        // at time step i (initialized at 0).
        vector<double> fi(problem.d, 0.);
        vector< vector<double> > hi(problem.d, vector<double>(problem.d,0.));
       
        // Loop for ensemble average
        for (int m = 0; m < solver.M; m++) {
            
            // Initialization of the fast variables
            yAux[0] = yInit;
            
            vector<double> drift;
            vector<double> diffu;

            
            // Euler-Maruyama method for the fast processes:
            for (unsigned int j = 0; j < yAux.size() - 1; j++) {

                drift = problem.drif(xt[i], yAux[j]);
                diffu = problem.diff(xt[i], yAux[j]);

                for (int k = 0; k < 2*problem.nf; k++)
                {
                    yAux[j+1][k] = yAux[j][k] + drift[k]*solver.micro_dt + 
                    diffu[k]*sqrt(solver.micro_dt)*distribution(generator);
                }
            }
            
            // Construction of auxiliary vector for efficiency.
            
            // sumsAux1 =approx= dax 
            vector< vector< vector<double> > > sumsAux1(solver.nt + solver.n, \
                    vector< vector<double> >(problem.d,vector<double>(problem.d,0.)));

            // sumsAux2 =approx= a
            vector< vector<double> > sumsAux2(solver.nt + solver.n, \
                    vector<double>(problem.d, 0.));

            // First component of each vector
            for (int index = 0; index <= solver.np; index++) {
                for (int i1 = 0; i1 < problem.d; i1++) {
                    for (int i2 = 0; i2 < problem.d; i2++) {
                        sumsAux1[0][i1][i2] += problem.dax(xt[i], yAux[index])[i1][i2];
                    }
                    sumsAux2[0][i1] += problem.a(xt[i], yAux[index])[i1];
                }
            }
            
            // Recurrence to obtain the other components
            for (int index = 1; index < solver.nt + solver.n; index++) {
                for (int i1 = 0; i1 < problem.d; i1++) {
                    for (int i2 = 0; i2 < problem.d; i2++) {
                        sumsAux1[index][i1][i2] = sumsAux1[index-1][i1][i2] + problem.dax(xt[i], \
                        yAux[index + solver.np])[i1][i2] - problem.dax(xt[i], yAux[index-1])[i1][i2];
                    }
                    sumsAux2[index][i1] = sumsAux2[index-1][i1] + problem.a(xt[i], \
                    yAux[index + solver.np])[i1] - problem.a(xt[i], yAux[index-1])[i1];
                }
            }

            // Drift term by the HMM
            vector<double> fim(problem.d, 0.);
            vector<double> fim1(problem.d, 0.);
            vector<double> fim2(problem.d, 0.);

            for (int j = solver.nt; j < solver.nt + solver.n; j++) {

                // evaluation of data
                vector< vector<double> > dya_j = problem.day(xt[i], yAux[j]);             

                // first term
                for (int i1 = 0; i1 < problem.d; i1++) {
                    for (int k = 0; k < problem.nf; k++) 
                        fim1[i1] += dya_j[i1][k]*yAux[j][problem.nf+k];
                }
                
                // second term: improved 
                for (int i1 = 0; i1 < problem.d; i1++) {
                    for (int i2 = 0; i2 < problem.d; i2++) {
                        fim2[i1] += solver.micro_dt*problem.a(xt[i], \
                                yAux[j])[i2]*sumsAux1[j][i1][i2];
                    }
                }
            }

            for (int i1 = 0; i1 < problem.d; i1++) {
                fim[i1] = fim1[i1] + fim2[i1];
                fim[i1] = fim[i1]/solver.n;
            }
            
            // Diffusion term by the HMM
            vector< vector<double> > him(problem.d, vector<double>(problem.d,0.));

            for (int i1 = 0; i1 < problem.d; i1++) {
                // TODO: Be sure that the indices are good here
                for (int i2 = 0; i2 < problem.d; i2++) { 
                    for (int j = solver.nt; j < solver.nt + solver.n; j++) {
                        him[i1][i2] += solver.micro_dt*problem.a(xt[i], \
                                yAux[j])[i1]*sumsAux2[j][i2];
                    }
                    him[i1][i2] = 2*him[i1][i2]/solver.n;
                }
            }
            
            // Incremental mean over the samples
            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    hi[i1][i2] = (hi[i1][i2]*m + him[i1][i2])/(m+1);
                }
                fi[i1] = (fi[i1]*m + fim[i1])/(m+1);
            }
        } 
        
        // Approximate diffusion coefficient
        cout << "   Diffusion coefficient before cholesky: " << endl;
        printMat(hi); cout << endl;

        hi = symmetric(hi);

        cout << "   Symmetric part of the diffusion" << endl;
        printMat(hi); cout << endl;

        hi = cholesky(hi);
        
        // Initial condition for next iteratation
        // and storage of y
        yInit = yAux[yAux.size()-1]; 

        // cout << "\ec";
        cout << "   Drift coefficient: " << endl;
        printVec(fi);
        cout << endl;

        cout << "   Diffusion coefficient: " << endl;
        printMat(hi);
        cout << endl;

        cout << "   Exact value of the drift coefficient: " << endl;
        printVec(problem.soldrif(xt[i]));
        cout << endl;

        cout << "   Exact value of the diffusion coefficient: " << endl;
        printMat(problem.soldiff(xt[i]));
        cout << endl;
        
        cout << "   Value of the slow variable x: " << endl;
        printVec(xt[i]);
        cout << endl;

        cout << "   Exact value of the slow variable: " << endl;
        printVec(x_exact[i]);
        cout << endl;
        
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
        cout << "   Total error up to the current iteration: " \
             << error << endl << endl;

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
     
    return error;
}
