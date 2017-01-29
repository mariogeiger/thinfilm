#include "../thinfilm2.hh"

#include <QString>
#include <QTime>


#include <iostream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    double theta;
    double lamda;
    thinfilm::complex nInc;
    thinfilm::complex nExi;
    vector<struct thinfilm::Layer> layers;

    if (argc >= 2 && QString(argv[1]).compare("--help") == 0) {
		cout << "theta(deg) lamda(nm) ninc kinc nexit kexit [ dlayer(nm) nlayer klayer ...]" << endl;
        return 0;
    }

    if (argc >= 8 && (argc - 8) % 3 == 0) {
        theta = QString(argv[1]).toDouble();
        lamda = QString(argv[2]).toDouble();
		nInc.real(QString(argv[3]).toDouble());
		nInc.imag(QString(argv[4]).toDouble());
		nExi.real(QString(argv[5]).toDouble());
		nExi.imag(QString(argv[6]).toDouble());

		for (int i = 7; i < argc; i += 3) {
            struct thinfilm::Layer layer;
            layer.thickness = QString(argv[i]).toDouble();
            layer.refractiveIndex = thinfilm::complex(QString(argv[i+1]).toDouble(),
                                                      -QString(argv[i+2]).toDouble());
            layers.push_back(layer);
        }
    } else {
        cout << "theta in deg :";
        cin >> theta;

        cout << "lamda in nm :";
        cin >> lamda;

		double x;
        cout << "n incident :";
		cin >> x;
		nInc.real(x);

        cout << "k incident :";
		cin >> x;
		nInc.imag(x);

        cout << "n exit :";
		cin >> x;
		nExi.real(x);

        cout << "k exit :";
		cin >> x;
		nExi.imag(x);

        char yes;
        do {
            struct thinfilm::Layer layer;

            cout << "d layer in nm :";
            cin >> layer.thickness;

            cout << "n layer :";
            double n;
            cin >> n;

            cout << "k layer :";
            double k;
            cin >> k;

            layer.refractiveIndex = thinfilm::complex(n, -k);

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
        cout << layers[i].thickness << " " << layers[i].refractiveIndex << endl;
    }
    cout << nExi << endl;

	double reflectanceP, reflectanceS;
	double transmittanceP, transmittanceS;

    QTime time;
    time.start();
    for (int i = 0; i < 100000; ++i)
		thinfilm::compute(cos(theta * (M_PI / 180.0)), lamda, nInc, nExi, layers,
						   &reflectanceP, &reflectanceS, &transmittanceP, &transmittanceS);

    cout << "execution time = " << time.elapsed()/100.0 << "us" << endl;

    cout << setprecision(30);
	cout << "                    P / S" << endl;
	cout << "reflectance    = " << reflectanceP * 100.0 << " / " << reflectanceS * 100.0 << " %" << endl;
	cout << "transmittance  = " << transmittanceP * 100.0 << " / " << transmittanceS * 100.0 << " %" << endl;

    return 0;
}
