#include "thinfilm.hh"
#include <iostream>

using namespace std;

int main() {

    thinfilm::complex incidentCosTheta;
    thinfilm::complex nIncident, nExit;
    vector<thinfilm::Layer> layers;

    cin >> incidentCosTheta;
    cin >> nIncident;
    cin >> nExit;

    size_t n;
    cin >> n;

    layers.resize(n);
    for (size_t i = 0; i < n; ++i) {
        cin >> layers[i].refractiveIndex;
        cin >> layers[i].thickness;
    }

    double reflectanceP;
    double reflectanceS;
    double transmittanceP;
    double transmittanceS;
    thinfilm::reflectance_transmittance(incidentCosTheta, nIncident,
                                        nExit, layers,
                                        &reflectanceP, &reflectanceS,
                                        &transmittanceP, &transmittanceS);

    cout.precision(15);
    cout << reflectanceP << " " << reflectanceS << " " << transmittanceP << " " << transmittanceS << endl;

    return 0;
}
