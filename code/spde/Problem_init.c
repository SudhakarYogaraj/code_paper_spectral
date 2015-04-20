#include "header.h"

void Problem::init() {

    // Final time
    this->t_end = 1.;

    // Constant apppearing in the SPDE
    this->nu = 0.0;
    
    // Number of fast processes
    this->nf = 4;

    // Number of slow processes
    this->d = 2;

    // Initial condition
    this->x0 = vector<double>(d,1.2);
}
