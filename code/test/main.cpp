// Includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <time.h>

#define PI 3.141592653589793238463

// Namespace
using namespace std;

vector<double> oneverynicetestfunction(vector<double> arg) {
    vector<double> barg = {4.,5.,6.};
    vector<double> carg(3,0.);
    arg = barg;
    /* arg[0] ++; arg[1] ++; arg[2] ++; */
    return arg;
}

int main(int argc, char *argv[])
{
    stringstream s; 
    s << "a" << "b" << "c" << "\0";
    cout << setw(60) <<  s.str() << endl;
    /* vector<double> a = {1.,2.,3.}; */
    /* vector<double> b; */
    /* b = oneverynicetestfunction(a); */

    /* cout << b[0] << b[1] << b[2] << endl; */
    /* cout << a[0] << a[1] << a[2] << endl; */
    
    return 0;
}
