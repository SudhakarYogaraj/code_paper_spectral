#include "io/io.hpp"
#include "toolbox/linear_algebra.hpp"

using namespace std;

const string output_path = "out/";

void writeToFile(string s, vec x) {
    std::ofstream output_file(output_path + s);
    ostream_iterator<double> out_it (output_file,"\n");
    copy (x.begin(), x.end(), out_it);
}

void writeMatToFile(string s, mat x) {
    std::ofstream fout(output_path + s);
    fout.precision(2);
    fout << scientific;

    if (fout.is_open()) {
        for (unsigned int i = 0; i < x.size(); i++) {
            for (unsigned int j = 0; j < x[0].size(); j++) {
                fout.width(15); fout <<  x[i][j];
            }
            fout << endl;
        }
    }
    else {
        cout << "There was a mistake while writing the output file" << endl;
        exit(0);
    }
}

void niceMat(mat a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void niceVec(vec a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
            cout << a[i] << " ";
            cout << endl;
    }
}
