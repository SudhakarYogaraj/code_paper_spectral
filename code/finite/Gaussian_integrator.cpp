#include "Gaussian_integrator.hpp"
#include "templates.hpp"

using namespace std;

static vector<double> nodes_4 = { 5.24647623275290317900e-01,    1.65068012388578455585e+00 };
static vector<double> weights_4 = { 8.04914090005512836482e-01,    8.13128354472451771398e-02 };
static vector<double> nodes_6 = { 2.35060497367449222281e+00,    4.36077411927616508688e-01, 1.33584907401369694976e+00 };
static vector<double> weights_6 = { 4.53000990550884564102e-03,    7.24629595224392524086e-01, 1.57067320322856643914e-01 };
static vector<double> nodes_10 = { 3.42901327223704608781e-01,    1.03661082978951365422e+00, 1.75668364929988177344e+00,    2.53273167423278979644e+00, 3.43615911883773760341e+00 };
static vector<double> weights_10 = { 6.10862633735325798764e-01,    2.40138611082314686412e-01, 3.38743944554810631376e-02,    1.34364574678123269223e-03, 7.64043285523262062930e-06 };
static vector<double> nodes_20 = { 2.45340708300901249903e-01,    7.37473728545394358719e-01, 1.23407621539532300786e+00,    1.73853771211658620678e+00, 2.25497400208927552311e+00,    2.78880605842813048055e+00, 3.34785456738321632688e+00,    3.94476404011562521040e+00, 4.60368244955074427298e+00,    5.38748089001123286199e+00 };
static vector<double> weights_20 = { 4.62243669600610089640e-01,    2.86675505362834129720e-01, 1.09017206020023320014e-01,    2.48105208874636108814e-02, 3.24377334223786183217e-03,    2.28338636016353967260e-04, 7.80255647853206369398e-06,    1.08606937076928169398e-07, 4.39934099227318055366e-10,    2.22939364553415129254e-13 };
static vector<double> nodes_30 = { 0.2011285765488714855458, 0.6039210586255523077782, 1.008338271046723461805, 1.415527800198188511941, 1.826741143603688038836, 2.243391467761504072473, 2.667132124535617200571, 3.099970529586441748689, 3.544443873155349886925, 4.003908603861228815228, 4.483055357092518341887, 4.988918968589943944486, 5.533147151567495725118, 6.138279220123934620395, 6.863345293529891581061};
static vector<double> weights_30 = { 0.3863948895418138625556, 0.2801309308392126674135, 0.1467358475408900997517, 0.05514417687023425116808, 0.01470382970482668351528, 0.00273792247306765846299, 3.48310124318685523421E-4, 2.9387252289229876415E-5, 1.57909488732471028835E-6, 5.10852245077594627739E-8, 9.178580424378528209E-10, 8.10618629746304420399E-12, 2.87860708054870606219E-14, 2.8103336027509037088E-17, 2.90825470013122622941E-21};

Gaussian_integrator::Gaussian_integrator(int nNodes, int nVars) {

    vector<double> nodes_1d(nNodes);
    vector<double> weights_1d(nNodes);

    switch (nNodes) {
        case 4:
            nodes_1d = nodes_4;
            weights_1d = weights_4;
           break;
        case 6:
            nodes_1d = nodes_6;
            weights_1d = weights_6;
            break;
        case 10:
            nodes_1d = nodes_10;
            weights_1d = weights_10;
            break;
        case 20:
            nodes_1d = nodes_20;
            weights_1d = weights_20;
            break;
        case 30:
            nodes_1d = nodes_30;
            weights_1d = weights_30;
            break;
        default: cout << "Invalid number of nodes for Gauss-hermite integration" << endl;
                 exit(0);
    }

    for (int i = 0; i < nNodes/2; ++i) {
        nodes_1d[i + nNodes/2] = -nodes_1d[i];
        weights_1d[i + nNodes/2] = weights_1d[i];
    }

    int nPoints = pow(nNodes, nVars);
    vector< vector<double> > x(nPoints, vector<double>(nVars,0.));
    vector<double> w(nPoints, 1.);
    for (int i = 0; i < nPoints; ++i) {
        int tmp = i;
        for (int j = 0; j < nVars; ++j) {
            x[i][j] = nodes_1d[tmp%nNodes];
            w[i] *= weights_1d[tmp%nNodes];
            tmp = tmp/nNodes;
        }
    }

    this->nodes = x;
    this->weights = w;
}
