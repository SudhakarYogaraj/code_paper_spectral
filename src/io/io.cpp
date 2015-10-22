#include "io/io.hpp"
#include "linear_algebra.hpp"

using namespace std;

const string output_path = "out/";

void writeToFile(string s, vector<double> x) {
    std::ofstream output_file(output_path + s);
    ostream_iterator<double> out_it (output_file,"\n");
    copy (x.begin(), x.end(), out_it);
}

void writeMatToFile(string s, vector< vector<double> > x) {
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

void printVec (vector<double> x) {
    stringstream vec;
    vec.precision(2);
    vec << scientific;
    vec << "    ";
    for (unsigned int i = 0; i < x.size(); i++) {
        vec << x[i] << "  ";
    }
    vec << "\0";
    cout << "|" << setw(101) << vec.str() << "|" <<  endl;
}

void print2Vecs (vector<double> x, vector<double> y) {
    stringstream vec1;
    vec1.precision(2);
    vec1 << scientific;
    vec1 << "    ";
    for (unsigned int i = 0; i < x.size(); i++) {
        vec1 << x[i] << "  ";
    }
    vec1 << "\0";

    stringstream vec2;
    vec2.precision(2);
    vec2 << scientific;
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

void niceMat(vector< vector<double> > a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void niceVec(vector<double> a) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
            cout << a[i] << " ";
            cout << endl;
    }
}


void progress_bar(double progress) {
    int width = 50;
    int position = width * progress;
    cout << "    Progress: [ " << setw(3) << int(progress * 100.0) << " % ] [";
    for (int i = 0; i < width; ++i) {
        if (i < position) cout << "=";
        else if (i == position) {
            cout << "x";
        }
        else cout << "·";
    }
    cout << "] \r";
    cout.flush();
}

void end_progress_bar() {
    progress_bar(1.0);
    cout << endl << endl;
}
