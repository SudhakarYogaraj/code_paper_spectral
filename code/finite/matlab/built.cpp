
vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = -sin(x[0])*(sin(y[0])+(y[0]*y[0]*y[0])*cos(y[0])-y[0]*cos(y[0]));
    return result;
}

vector< vector<double> > Problem::grad_h(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->nf, vector<double>(this->nf, 0.));
    result[0][0] = -cos(x[0])*sin(y[0]);
    return result;
}

vector<double> Problem::grad(vector<double> x, vector<double> y){
    vector<double> result(this->nf);
    result[0] = -y[0]+y[0]*y[0]*y[0];
    return result;
}

vector<double> Problem::phi(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = cos(x[0])*sin(y[0]);
    return result;
}

vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = cos(x[0])*(sin(y[0])+(y[0]*y[0]*y[0])*cos(y[0])-y[0]*cos(y[0]));
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = cos(x[0])*cos(y[0]);
    return result;
}

double Problem::potential(vector<double> x, vector<double> y) {
    double result = 0.;

    return result;
}

double Problem::linearTerm(vector<double> x, vector<double> y){
    double result;
    result = pow(y[0]-y[0]*y[0]*y[0],2.0)*(-1.0/4.0)+(y[0]*y[0])*(3.0/2.0)-1.0/2.0;
    return result;
}

double Problem::rho(vector<double> x, vector<double> y) {
    double result = 0.;
    result = exp((y[0]*y[0])*(1.0/2.0)-(y[0]*y[0]*y[0]*y[0])*(1.0/4.0))*2.560729512196574E-1;
    return result;
}
