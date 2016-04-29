#ifndef IO_H
#define IO_H

#include <iostream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
void writeToFile(std::string s, std::vector<double> x);
void writeMatToFile(std::string s, std::vector< std::vector<double> > x);
void niceMat(std::vector< std::vector<double> > a);
void niceVec(std::vector<double> a);
#endif
