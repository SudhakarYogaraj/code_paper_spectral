#include "header.h"

void Solver::set(double p, int M)
{
    // Precision parameter
    this->p = p;

    // Order of the micro solver
    this->l = 1;

    // Macro time-step
    this->macro_dt = .01;

    // Micro time-step
    
    // If KS
    // this->micro_dt = 0.05*pow(2.,-this->p/this->l);
    
    // If Burgers
    this->micro_dt = pow(2.,-this->p/this->l);

    // Number of micro time-steps taken into account
    // in the average to obtain the coefficients of the 
    // effective equation at each macro time-step.
    this->n = (int) 10*pow(2,this->p*(2+1./this->l));

    // Number of micro time-steps that are not taken 
    // into account in the averages
    this->nt = 16;

    // Number of micre time-steps used for the discretization
    // of the integrals in time. 
    this->np = (int) pow(2.,this->p/this->l)*this->p;

    // Number of replicas of the fast process
    this->M = M;
}
