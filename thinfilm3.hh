#ifndef THINFILM_H
#define THINFILM_H

#include <vector>
#include <complex>
#include <cmath>

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
    Matrix22(complex m11, complex m12, complex m21, complex m22) :
        m11(m11), m12(m12), m21(m21), m22(m22)
    {
    }

    Matrix22& operator*=(const Matrix22& b) {
        complex x = m11 * b.m11 + m12 * b.m21;
        m12 = m11 * b.m12 + m12 * b.m22;
        m11 = x;
        x = m21 * b.m11 + m22 * b.m21;
        m22 = m21 * b.m12 + m22 * b.m22;
        m21 = x;
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

const Matrix22 operator*(Matrix22 a, const Matrix22& b)
{
    return a *= b;
}

void compute(double incidentCosTheta, complex nIncident, complex nExit, const std::vector<Layer>& layers,
             double* reflectanceP, double* reflectanceS, double* transmittanceP, double* transmittanceS) {
    Matrix22 productMatrixP(1.0, 0.0, 0.0, 1.0);
    Matrix22 productMatrixS(1.0, 0.0, 0.0, 1.0);

    for (std::size_t i = 0; i < layers.size(); ++i) {
        complex n = layers[i].refractiveIndex;
        complex cosTheta = std::sqrt(1.0 - (1.0 - incidentCosTheta*incidentCosTheta) * (nIncident*nIncident) / (n*n));
        complex deltaLayer = 2.0 * M_PI * n * layers[i].thickness * cosTheta;
        complex c = std::cos(deltaLayer);
        complex s = std::sin(deltaLayer);

        const complex j(0.0, 1.0);
        productMatrixP *= Matrix22(c, -j * s * cosTheta / n, -j * s * n / cosTheta, c);
        productMatrixS *= Matrix22(c, -j * s / cosTheta / n, -j * s * n * cosTheta, c);
    }

    complex exitCosTheta = std::sqrt(1.0 - (1.0 - incidentCosTheta*incidentCosTheta) * (nIncident*nIncident) / (nExit*nExit));

    complex bP = productMatrixP.m11 * exitCosTheta + productMatrixP.m12 * nExit; // has been multiplied by exitCosTheta to avoid division by 0
    complex cP = productMatrixP.m21 * exitCosTheta + productMatrixP.m22 * nExit;

    complex bS = productMatrixS.m11 + productMatrixS.m12 * nExit * exitCosTheta;
    complex cS = productMatrixS.m21 + productMatrixS.m22 * nExit * exitCosTheta;

    double numerator = 4.0 * std::abs(nExit * exitCosTheta * nIncident * incidentCosTheta);
    double denomP = std::norm(bP * nIncident + cP * incidentCosTheta);
    double denomS = std::norm(bS * nIncident * incidentCosTheta + cS);

    if (denomP != 0.0) {
        *reflectanceP = std::norm(bP * nIncident - cP * incidentCosTheta) / denomP;
        *transmittanceP = numerator / denomP;
    }
    if (denomS != 0.0) {
        *reflectanceS = std::norm(bS * nIncident * incidentCosTheta - cS) / denomS;
        *transmittanceS = numerator / denomS;
    }
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
