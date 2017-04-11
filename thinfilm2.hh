/***********************************************************
   License: GNU Public License

   (c) Copyright LESO-PB 2011

   Address:
    EPFL ENAC IIC LESO-PB
    CH-1015 Lausanne

   Author: Mario Geiger
   E-Mail: geiger.mario@gmail.com

   Description: Header file to calculate the
                reflectance, transmittance,
                psi and delta of a
                multilayer coating.

   Limitations: Tested with Film Wizard and
                Tfcalc. absorbing incident
                medium, psi and delta are
                not verified.
***********************************************************/

#ifndef THINFILM_H
#define THINFILM_H

#include <complex>
#include <cstdio>
#include <vector>


// ----------------------------------------------------------------------------
namespace thinfilm {

typedef std::complex<double> complex;

// ----------------------------------------------------------------------------
/**
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

**/
// ----------------------------------------------------------------------------


struct Layer {
  // express in the same unit as the wavelength
  double thickness;

  // with k non-negative
  // example (1.5, +0.001)
  complex refractiveIndex;
};

struct Matrix22 {
  Matrix22() {}
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

// ----------------------------------------------------------------------------
inline void compute(
    // cosine of insident angle
    complex incidentCosTheta,
    // wavelength of light (same unit as layers thickness)
    double lambda,
    // complex index of refraxion of insident medium,
    // k value must be non-negative
    complex nIncident,
    // complex index of refraxion of exit medium,
    // k value must be non-negative
    complex nExit,

    // array of layers
    // in the array the layers are presented from
    //           the insident one to the exit one
    std::vector<Layer> layers,

    double *reflectanceP,
    double *reflectanceS,
    double *transmittanceP,
    double *transmittanceS)
{
  int layersAmount = layers.size();

  /**
      always for each P and S polarization

      1.1 : we calculate the admittance of the incident and exit medium

      1.2 : the big unique loop
            calculate :
            admittance of the layer
            shift phase of the layer
            characteristic matrix
            product of the layer matrix and the main matrix

      2   : we do strange equations with the main matrix to extract
            reflectance transmittance and absorptance
      **/

  // ------------------------------------------------------------------------
  // ---------- 1.1 ----------


  // calculate admittance of the incident medium for P and S polarization
  const complex admittanceIncidentP = nIncident / incidentCosTheta;
  const complex admittanceIncidentS = nIncident * incidentCosTheta;


  // calculate cosine of exit medium with snell law optimized
  //  sin(theta0)    n0  =   sin(theta1)    n1
  //       s0        n0  =        s1        n1
  // sqrt( 1 - c0² ) n0  =  sqrt( 1 - c1² ) n1
  // and solve for c1 :
  // c1 = sqrt( c0² n0² - n0² + n1² ) / n1
  const complex squareIncidentCosTheta = incidentCosTheta * incidentCosTheta;
  const complex ratioIncidentExit = nIncident / nExit;
  const complex exitCosTheta =
      std::sqrt(1. - (1. - squareIncidentCosTheta) * ratioIncidentExit*ratioIncidentExit);

  // calculate the exit medium admitance P and S
  const complex admittanceExitP = nExit / exitCosTheta;
  const complex admittanceExitS = nExit * exitCosTheta;








  // ------------------------------------------------------------------------
  // ---------- 1.2 ----------

  // mains matrix, the product of each layer matrix
  Matrix22 productMatrixP(1.0, 0.0, 0.0, 1.0);
  Matrix22 productMatrixS(1.0, 0.0 ,0.0, 1.0);
  // initialized to Identity Matrix
  /*
      ( 1    0 )
      (        )
      ( 0    1 )
     */



  // the matrix of the layer
  Matrix22 layerMatrixP;
  Matrix22 layerMatrixS;



  // i : for each layer
  for (int i = 0; i < layersAmount; ++i) {

    // layerRefractiveIndex complex number : n + ik

    // snell law
    const complex ratio = nIncident / layers[i].refractiveIndex;
    const complex layerCosTheta =
        std::sqrt(1. - (1. - squareIncidentCosTheta) * ratio*ratio);



    // then calculate the admittance of the layer P and S
    const complex admittanceLayerP =
        layers[i].refractiveIndex / layerCosTheta;
    const complex admittanceLayerS =
        layers[i].refractiveIndex * layerCosTheta;

    // now the delta dephasing of the layer
    const complex deltaLayer =
        2.0 * M_PI * layers[i].refractiveIndex * layers[i].thickness * layerCosTheta / lambda;
    // juillet 2015 : il faut bien multiplier par le cos => il faut regarder les front d'ondes !
    // the thickness layer, is need to be the same unit of lambda


    // create the matrix layer
    const complex c = cos(deltaLayer);
    const complex s = sin(deltaLayer);

    layerMatrixP.m11 = c;
    layerMatrixP.m12 = -complex(0.0, 1.0) * s / admittanceLayerP;
    layerMatrixP.m21 = -complex(0.0, 1.0) * s * admittanceLayerP;
    layerMatrixP.m22 = c;


    layerMatrixS.m11 = c;
    layerMatrixS.m12 = -complex(0.0, 1.0) * s / admittanceLayerS;
    layerMatrixS.m21 = -complex(0.0, 1.0) * s * admittanceLayerS;
    layerMatrixS.m22 = c;


    // now we need to make the product of the main matrix with the layer matrix

    // do the matrix product
    productMatrixP *= layerMatrixP;
    productMatrixS *= layerMatrixS;
  }







  // ------------------------------------------------------------------------
  // ----------  2  ----------

  const complex bP = productMatrixP.m11 + productMatrixP.m12 * admittanceExitP;
  const complex cP = productMatrixP.m21 + productMatrixP.m22 * admittanceExitP;



  const complex bS = productMatrixS.m11 + productMatrixS.m12 * admittanceExitS;
  const complex cS = productMatrixS.m21 + productMatrixS.m22 * admittanceExitS;



  // calculate the reflectance
  const complex reflectionCoefficientP =
      (admittanceIncidentP - cP / bP) / (admittanceIncidentP + cP / bP);

  const complex reflectionCoefficientS =
      (admittanceIncidentS - cS / bS) / (admittanceIncidentS + cS / bS);


  // norm returns the norm value of the complex number : norm(3+4i) = 25
  *reflectanceP = std::norm(reflectionCoefficientP);
  *reflectanceS = std::norm(reflectionCoefficientS);

  // warning ! : the transmittance is correct only if kIncident == 0
  if (imag(nIncident) != 0.0) {
    fprintf(stderr, "%s:%d : warning ! the transmittance"
                    " is maybe false (incident k != 0)", __FILE__, __LINE__);
  }

  // the transmittance
  // pour les std::abs je n'ai aucune idée.
  // Dans le cas où l'indice de refraction incident et sortie
  // sont réels j'ai fais le calcul dans doc/calculs.pdf.
  *transmittanceP = std::abs(admittanceExitP / admittanceIncidentP) * std::norm(2.0 / (bP + cP / admittanceIncidentP));
  *transmittanceS = std::abs(admittanceExitS / admittanceIncidentS) * std::norm(2.0 / (bS + cS / admittanceIncidentS));
}


} // end namespace

#endif // THINFILM_H
