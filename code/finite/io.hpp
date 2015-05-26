#ifndef IO_H
#define IO_H

#include <iostream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>

void progress_bar(double progress);
void end_progress_bar();
void printVec (std::vector<double> x);
void printMat (std::vector< std::vector<double> > x);
void print2Vecs (std::vector<double> x, std::vector<double> y);
void print2Mats (std::vector< std::vector<double> > x, std::vector< std::vector<double> > y);
void writeToFile(std::string s, std::vector<double> x);
void writeMatToFile(std::string s, std::vector< std::vector<double> > x);
#endif
