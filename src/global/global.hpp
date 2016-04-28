#ifndef GLOBAL_H
#define GLOBAL_H

#define VERBOSE 1
#define SUMMARY 1
#define DEBUG 0

#include<vector>
#include<armadillo>

// Definition of types
namespace std {
    typedef std::vector< std::vector< std::vector<double> > > cube;
    typedef std::vector< std::vector<double> > mat;
    typedef std::vector<double> vec;
}

// Structures
struct SDE_coeffs {
    std::vec drif;
    std::mat diff;
};

#endif
