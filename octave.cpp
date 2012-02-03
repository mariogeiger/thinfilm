#include <QString>
#include <QTime>

#include "thinfilm.hpp"

#include <iostream>

using namespace blitz;
using namespace std;

int main(int argc, char *argv[])
{
    double theta;
    double lamda;
    double polar;
    thinfilm::complex nInc;
    thinfilm::complex nExi;
    vector<TinyVector<double, 3> > layers;

    if (argc >= 8 && (argc - 8) % 3 == 0) {
        theta = QString(argv[1]).toDouble();
        lamda = QString(argv[2]).toDouble();
        polar = QString(argv[3]).toDouble();
        nInc.real() = QString(argv[4]).toDouble();
        nInc.imag() = QString(argv[5]).toDouble();
        nExi.real() = QString(argv[6]).toDouble();
        nExi.imag() = QString(argv[7]).toDouble();

        for (int i = 8; i < argc; ++i) {
            TinyVector<double, 3> layer;
            layer[(i-8) % 3] = QString(argv[i]).toDouble(); i++;
            layer[(i-8) % 3] = QString(argv[i]).toDouble(); i++;
            layer[(i-8) % 3] = QString(argv[i]).toDouble();
            layers.push_back(layer);
        }
    } else {
        return 1;
    }

    // n-ik
    nInc = conj(nInc);
    nExi = conj(nExi);

    double reflectance;
    double transmittance;
    double absorptance;
//    double psi, delta;

    thinfilm::simulate(cos(theta * (M_PI / 180.0)), lamda, polar * (M_PI / 180.0), nInc, nExi, layers,
                       &reflectance, &transmittance, &absorptance);

    cout << reflectance << " " << transmittance << " " << absorptance << endl;

    return 0;
}
