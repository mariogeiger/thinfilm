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

    if (argc >= 2 && QString(argv[1]).compare("--help") == 0) {
        cout << "theta(deg) lamda(nm) pol(deg0S) ninc kinc nexit kexit [ dlayer(nm) nlayer klayer ...]" << endl;
        return 0;
    }

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
        cout << "theta in deg :";
        cin >> theta;

        cout << "lamda in nm :";
        cin >> lamda;

        cout << "polarization in deg 0isS :";
        cin >> polar;

        cout << "n incident :";
        cin >> nInc.real();

        cout << "k incident :";
        cin >> nInc.imag();

        cout << "n exit :";
        cin >> nExi.real();

        cout << "k exit :";
        cin >> nExi.imag();

        char yes;
        do {
            TinyVector<double, 3> layer;

            cout << "d layer in nm :";
            cin >> layer[0];

            cout << "n layer :";
            cin >> layer[1];

            cout << "k layer :";
            cin >> layer[2];

            layers.push_back(layer);

            cout << "add an other layer ? y/n :";
            cin >> yes;
        } while (yes == 'y');
    }

    // n-ik
    nInc = conj(nInc);
    nExi = conj(nExi);

    cout << nInc << endl;
    for (uint i = 0; i < layers.size(); ++i) {
        cout << layers[i] << endl;
    }
    cout << nExi << endl;

    double reflectance;
    double transmittance;
    double absorptance;
    double psi, delta;

    //    QTime time;
    //    time.start();
    //for (int i = 0; i < 100000; ++i)
    thinfilm::simulate(cos(theta * (M_PI / 180.0)), lamda, polar * (M_PI / 180.0), nInc, nExi, layers,
                       &reflectance, &transmittance, &absorptance, &psi, &delta);

    //    cout << "execution time = " << time.elapsed()/100.0 << "us" << endl;

    cout << "reflectance  =  " << reflectance * 100.0 << "%" << endl;
    cout << "transmittance = " << transmittance * 100.0 << "%" << endl;
    cout << "absorptance  =  " << absorptance * 100.0 << "%" << endl;

    cout << "psi  =  " << psi * (180.0 / M_PI) << "°" << endl;
    cout << "delta = " << delta * (180.0 / M_PI) << "°" << endl;

    return 0;
}
