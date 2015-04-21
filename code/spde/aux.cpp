#include "header.h"

const string output_path = "/home/urbain/mres/output/spde/";

// Function to write vector to file
void writeToFile(string s, vector<double> x) {
    ofstream output_file(output_path + s);
    ostream_iterator<double> out_it (output_file,"\n");
    copy (x.begin(), x.end(), out_it);
}

// Write matrix to a file
void writeMatToFile(string s, vector< vector<double> > x) {
    ofstream fout(output_path + s);
    fout.precision(5);
    fout << scientific;

    if (fout.is_open()) {
        for (int i = 0; i < x.size(); i++) {
            for (int j = 0; j < x[0].size(); j++) {
                fout.width(15); fout <<  x[i][j];
            }
            fout << endl;
        }
    }
    // else {
    //     cout << "There was a mistake while writing the output file" << endl;
    //     exit(0);
    // }
}

// Delta function
double delta(int a, int b) {
	if (fabs(b-a) < 0.1)
		return 1.;
	else return 0.;
}

// Symmetric part of a matrix
vector< vector<double> > symmetric(vector< vector<double> > A) {
    int n = A.size();
    vector< vector <double> > result(n,vector<double>(n,0.));

    for (int i = 0.; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            result[i][j] = 0.5*(A[i][j] + A[j][i]);
            result[j][i] = result[i][j];
        }
    }
    return result;
}

// Cholesky factorization of a matrix
vector< vector<double> > cholesky(vector< vector<double> > A) {
    int n = A.size();
    vector< vector <double> > L(n,vector<double>(n,0.));

    for (int i = 0.; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.;
            for (int k = 0; k < j; k++) {
                sum += L[i][k]*L[j][k];
            }
            if (i == j) {
                L[i][i] = sqrt(A[i][i]-sum);
            }
            else {
                if(fabs(L[j][j])>1e-14) {
                    L[i][j] = (A[i][j] - sum)/L[j][j];
                }
                else {
                    L[i][j] = 0.;
                }
            }
        }
    }
    return L;
}

// Print vec
void printVec (vector<double> x) {
    cout << "   ";
    for (int i = 0; i < x.size(); i++) {
        cout << "\t" << x[i];
    }
    cout << endl;
}

// Print mat
void printMat (vector< vector<double> > x) {
    for (int i = 0; i < x.size(); i++) {
        cout << "   ";
        for (int j = 0; j < x[0].size(); j++) {
            cout << "\t" << x[i][j];
        }
        cout << endl;
    }
}

// Norm vec
double normVec (vector<double> x) {
    double result = 0.;
    for (int i = 0; i < x.size(); i++) {
        result += x[i]*x[i];
    }
    return sqrt(result);
}

// Norm mat
double normMat (vector< vector<double> > x) {
    double result = 0.;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[0].size(); j++) {
            result += x[i][j]*x[i][j];
        }
    }
    return sqrt(result);
}
