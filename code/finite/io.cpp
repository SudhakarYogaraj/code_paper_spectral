#include "io.hpp"

using namespace std;

const string output_path = "/home/urbain/output/";

void writeToFile(string s, vector<double> x) {
    std::ofstream output_file(output_path + s);
    ostream_iterator<double> out_it (output_file,"\n");
    copy (x.begin(), x.end(), out_it);
}

void writeMatToFile(string s, vector< vector<double> > x) {
    std::ofstream fout(output_path + s);
    fout.precision(5);
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

void printVec (vector<double> x) {
    stringstream vec;
    vec << "    ";
    for (unsigned int i = 0; i < x.size(); i++) {
        vec << x[i] << "  ";
    }
    vec << "\0";
    cout << "|" << setw(101) << vec.str() << "|" <<  endl;
}

void print2Vecs (vector<double> x, vector<double> y) {
    stringstream vec1;
    vec1 << "    ";
    for (unsigned int i = 0; i < x.size(); i++) {
        vec1 << x[i] << "  ";
    }
    vec1 << "\0";

    stringstream vec2;
    vec2 << "    ";
    for (unsigned int i = 0; i < y.size(); i++) {
        vec2 << y[i] << "  ";
    }
    vec2 << "\0";

    cout << "|" << setw(50) << vec1.str() << "|" << setw(50) << vec2.str() << "|" <<  endl;
}

void printMat (vector< vector<double> > x) {
    for (unsigned int i = 0; i < x.size(); i++) {
        printVec(x[i]);
    }
}

void print2Mats (vector< vector<double> > x, vector< vector<double> > y) {
    for (unsigned int i = 0; i < x.size(); i++) {
        print2Vecs(x[i], y[i]);
    }
}
