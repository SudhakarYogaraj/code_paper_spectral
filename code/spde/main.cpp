#include "header.h"

// Main function
int main(int argc, char* argv[])
{
	// Initialization of the problem
	Problem problem; 
	problem.init();

    // Values of the precision parameter
    vector<double> p_values = {2,2.5,3,3.5,4,4.5,5,5.5,6};

	// Vector of the log of the error
    vector<double> errors(p_values.size(), 0.);
        
    // Number of replicas of the fast process
    int M = 1;
    
    // Random variable for generator
    int seed = time(NULL);

	// Initialization of the solver
	Solver solver;

	for (unsigned int i = 0; i < p_values.size(); i++) {

		// Setting to solver to match the precision parameter p_value[i];
		solver.set(p_values[i],M);
		
		// Error corresponding to the current value of p
		double error = solve(seed, problem,solver);

		// log2 of the error, used to produce a plot
		errors[i] = log2(error);
	}

	for (unsigned int i = 0; i < p_values.size(); i++) {
		cout << "Error for p = " << p_values[i] << ": " << errors[i] << endl; 
	}
	cout << endl << endl;

    // writeToFile("errors.dat", errors);
    // writeToFile("p_values.dat", p_values); 
}
