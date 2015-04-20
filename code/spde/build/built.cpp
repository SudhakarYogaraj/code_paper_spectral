vector<double> Problem::soldrif(vector<double> x) {
    vector<double> result(this->d,0.);
    result[0] = x[0]*((x[0]*x[0])*1.188E5+(x[1]*x[1])*1.188E5-7.87E2)*(-7.014590347923681E-7);
    result[1] = x[1]*5.520482603815937E-4-(x[0]*x[0])*x[1]*(1.0/1.2E1)-(x[1]*x[1]*x[1])*(1.0/1.2E1);
    return result;
}

vector< vector<double> > Problem::soldiff(vector<double> x) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = sqrt((x[0]*x[0])*(1.0/3.24E2)+(x[1]*x[1])*(1.0/5.76E2)+2.104377104377104E-6);
    result[1][0] = x[0]*x[1]*1.0/sqrt((x[0]*x[0])*(1.0/3.24E2)+(x[1]*x[1])*(1.0/5.76E2)+2.104377104377104E-6)*(-1.350308641975309E-3);
    result[1][1] = 1.0/sqrt((x[0]*x[0])*4.4E3+(x[1]*x[1])*2.475E3+3.0)*sqrt((x[0]*x[0])*(x[1]*x[1])*1.787712191358025E1-pow(fabs(x[0]*x[1]),2.0)*2.599344135802469+(x[0]*x[0])*(1.0/9.6E1)+(x[1]*x[1])*1.218894675925926E-2+(x[0]*x[0]*x[0]*x[0])*(2.75E2/3.6E1)+(x[1]*x[1]*x[1]*x[1])*(2.75E2/3.6E1)+3.551136363636364E-6);
    return result;
}

vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = x[0]*y[0]*(-1.0/2.0)-x[1]*y[1]*(1.0/2.0)-y[0]*y[2]*(1.0/2.0)-y[1]*y[3]*(1.0/2.0);
    result[1] = x[0]*y[1]*(-1.0/2.0)+x[1]*y[0]*(1.0/2.0)-y[0]*y[3]*(1.0/2.0)+y[1]*y[2]*(1.0/2.0);
    return result;
}

vector<double> Problem::a_nu(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = nu*x[0];
    result[1] = nu*x[1];
    return result;
}

vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) { 
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = y[0]*(-1.0/2.0);
    result[0][1] = y[1]*(-1.0/2.0);
    result[1][0] = y[1]*(-1.0/2.0);
    result[1][1] = y[0]*(1.0/2.0);
    return result;
}

vector< vector<double> > Problem::day(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->nf,0.));
    result[0][0] = x[0]*(-1.0/2.0)-y[2]*(1.0/2.0);
    result[0][1] = x[1]*(-1.0/2.0)-y[3]*(1.0/2.0);
    result[0][2] = y[0]*(-1.0/2.0);
    result[0][3] = y[1]*(-1.0/2.0);
    result[1][0] = x[1]*(1.0/2.0)-y[3]*(1.0/2.0);
    result[1][1] = x[0]*(-1.0/2.0)+y[2]*(1.0/2.0);
    result[1][2] = y[1]*(1.0/2.0);
    result[1][3] = y[0]*(-1.0/2.0);
    return result;
}

vector<double> Problem::drif(vector<double> x, vector<double> y) {
    vector<double> result(2*this->nf,0.);
    result[0] = y[0]*-3.0;
    result[1] = y[1]*-3.0;
    result[2] = y[2]*-8.0;
    result[3] = y[3]*-8.0;
    result[4] = y[4]*-3.0-x[0]*y[2]-x[1]*y[3]+(x[0]*x[0])*(1.0/2.0)-(x[1]*x[1])*(1.0/2.0);
    result[5] = y[5]*-3.0+x[0]*x[1]-x[0]*y[3]+x[1]*y[2];
    result[6] = y[6]*-8.0+x[0]*y[0]*(3.0/2.0)-x[1]*y[1]*(3.0/2.0);
    result[7] = y[7]*-8.0+x[0]*y[1]*(3.0/2.0)+x[1]*y[0]*(3.0/2.0);
    return result;
}

vector<double> Problem::diff(vector<double> x, vector<double> y) {
    vector<double> result(2*nf,0.);
    result[0] = 1.0/3.0;
    result[1] = 1.0/4.0;
    result[2] = 1.0/5.0;
    return result;
}
