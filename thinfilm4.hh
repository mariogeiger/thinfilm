#ifndef THINFILM_H
#define THINFILM_H

#include <complex>
#include <cmath>
#include <vector>

// This version is based on https://arxiv.org/abs/1603.02720

namespace thinfilm {

typedef std::complex<double> complex;

struct Layer {
    // Thickness of the layer in unit of wavelength
    // thickness == 1  <=> thickness = 1 wavelength
    double thickness;

    // with k non-negative
    // example (1.5, +0.001)
    complex refractiveIndex;
};

struct Matrix22 {
    Matrix22() {
    }
    Matrix22(complex m11, complex m12, complex m21, complex m22)
        : m11(m11), m12(m12), m21(m21), m22(m22) {
    }

    Matrix22 &operator*=(const Matrix22 &b) {
        complex x = m11 * b.m11 + m12 * b.m21;
        m12 = m11 * b.m12 + m12 * b.m22;
        m11 = x;
        x = m21 * b.m11 + m22 * b.m21;
        m22 = m21 * b.m12 + m22 * b.m22;
        m21 = x;
        return *this;
    }

    Matrix22 &operator/=(complex b) {
        m11 /= b;
        m12 /= b;
        m21 /= b;
        m22 /= b;
        return *this;
    }

    /*
        ( m11   m12 )
        (           )
        ( m21   m22 )
     */
    complex m11;
    complex m12;
    complex m21;
    complex m22;
};

const Matrix22 operator*(Matrix22 a, const Matrix22 &b) {
    return a *= b;
}

std::pair<Matrix22, Matrix22> transfer_matrix_PS(double incidentCosTheta, complex nIncident,
                                              complex nExit, const std::vector<Layer> &layers) {
    // nIncident * incidentSinTheta MUST be real, see
    // https://arxiv.org/abs/1603.02720

    Matrix22 MatP(nIncident, incidentCosTheta, nIncident, -incidentCosTheta);
    MatP /= 2.0 * nIncident * incidentCosTheta;
    Matrix22 MatS(nIncident * incidentCosTheta, 1.0, nIncident * incidentCosTheta,
                  -1.0);
    MatS /= 2.0 * nIncident * incidentCosTheta;

    for (std::size_t i = 0; i < layers.size(); ++i) {
        complex n = layers[i].refractiveIndex;
        complex cosTheta = std::sqrt(1.0 -
                                     (1.0 - incidentCosTheta * incidentCosTheta) *
                                     (nIncident * nIncident) / (n * n));

        complex deltaLayer = 2.0 * M_PI * n * layers[i].thickness * cosTheta;
        complex c = std::cos(deltaLayer);
        complex s = std::sin(deltaLayer);

        const complex j(0.0, 1.0);
        MatP *= Matrix22(c, -j * s * cosTheta / n, -j * s * n / cosTheta, c);
        MatS *= Matrix22(c, -j * s / cosTheta / n, -j * s * n * cosTheta, c);
    }

    complex exitCosTheta =
        std::sqrt(1.0 -
                  (1.0 - incidentCosTheta * incidentCosTheta) *
                  (nIncident * nIncident) / (nExit * nExit));

    MatP *= Matrix22(exitCosTheta, exitCosTheta, nExit, -nExit);
    MatS *= Matrix22(1.0, 1.0, nExit * exitCosTheta, -nExit * exitCosTheta);

    return std::pair<Matrix22, Matrix22>(MatP, MatS);
}

void reflectance_transmittance(double incidentCosTheta, complex nIncident,
                               complex nExit, const std::vector<Layer> &layers,
                               double *reflectanceP, double *reflectanceS,
                               double *transmittanceP, double *transmittanceS) {

    std::pair<Matrix22, Matrix22> matricies = transfer_matrix_PS(incidentCosTheta, nIncident, nExit, layers);
    Matrix22 MatP = matricies.first;
    Matrix22 MatS = matricies.second;

    complex rP = MatP.m21 / MatP.m11;
    complex rS = MatS.m21 / MatS.m11;

    complex tP = 1.0 / MatP.m11;
    complex tS = 1.0 / MatS.m11;

    *reflectanceP = std::norm(rP);
    *reflectanceS = std::norm(rS);

    *transmittanceP = std::norm(tP) * std::real(nExit * std::conj(exitCosTheta)) /
                      std::real(nIncident * std::conj(incidentCosTheta));
    *transmittanceS = std::norm(tS) * std::real(nExit * exitCosTheta) /
                      std::real(nIncident * incidentCosTheta);
}

void reflectance(double incidentCosTheta, complex nIncident, complex nExit,
                 const std::vector<Layer> &layers, double *reflectanceP,
                 double *reflectanceS) {

    std::pair<Matrix22, Matrix22> matricies = transfer_matrix_PS(incidentCosTheta, nIncident, nExit, layers);
    Matrix22 MatP = matricies.first;
    Matrix22 MatS = matricies.second;

    complex rP = MatP.m21 / MatP.m11;
    complex rS = MatS.m21 / MatS.m11;

    *reflectanceP = std::norm(rP);
    *reflectanceS = std::norm(rS);
}

}
#endif

/*
 \               ^
 \             /  Reflectance
 \           /
 \         /
 \       /
 \     /
        v   /

        Incident medium
            n + i*k
   -------------------------------
        Layer 0 medium
         d0 (n0 + i*k0)
   -------------------------------
             ...

   -------------------------------
        Layer j medium
         dj (nj + i*kj)
   -------------------------------
             ...

   -------------------------------
        Layer n medium
         dn (nn + i*kn)
   -------------------------------
         Exit medium
            n + i*k
 \
 \
 \
                       v  Transmittance

 */
