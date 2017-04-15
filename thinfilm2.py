# pylint: disable=C,R
import cmath
import numpy as np

# layers = [(index of refraxion, thickness / lambda in vaccum), second layer, ...]
def compute(incidentCosTheta, nIncident, nExit, layers):

    # sin(Theta)*nIncident must be purely real

    MatrixP = np.array([[nIncident, incidentCosTheta], [nIncident, -incidentCosTheta]]) / (2 * nIncident * incidentCosTheta)
    MatrixS = np.array([[nIncident * incidentCosTheta, 1], [nIncident * incidentCosTheta, -1]]) / (2 * nIncident * incidentCosTheta)

    for layer in layers:
        n = layer[0]
        thickness = layer[1] # express in mutiple of the wavelength

        cosTheta = cmath.sqrt(1. - (1. - incidentCosTheta ** 2) * (nIncident / n) ** 2)

        deltaLayer = 2.0 * cmath.pi * n * thickness * cosTheta
        c = cmath.cos(deltaLayer)
        s = cmath.sin(deltaLayer)

        MatrixP = np.dot(MatrixP, np.array([[c, -1j * s * cosTheta / n], [-1j * s * n / cosTheta, c]]))
        MatrixS = np.dot(MatrixS, np.array([[c, -1j * s / cosTheta / n], [-1j * s * n * cosTheta, c]]))
    
    exitCosTheta = cmath.sqrt(1. - (1. - incidentCosTheta ** 2) * (nIncident / nExit) ** 2)
    
    MatrixP = np.dot(MatrixP, np.array([[exitCosTheta, exitCosTheta], [nExit, -nExit]]))
    MatrixS = np.dot(MatrixS, np.array([[1, 1], [nExit * exitCosTheta, -nExit * exitCosTheta]]))
    
    rP = MatrixP[1,0] / MatrixP[0,0]
    rS = MatrixS[1,0] / MatrixS[0,0]
    
    tP = 1 / MatrixP[0,0]
    tS = 1 / MatrixS[0,0]
    
    reflectanceP = abs(rP)**2
    reflectanceS = abs(rS)**2
    
    transmittanceP = abs(tP)**2 * (nExit * np.conj(exitCosTheta)).real / (nIncident * np.conj(incidentCosTheta)).real
    transmittanceS = abs(tS)**2 * (nExit * exitCosTheta).real / (nIncident * incidentCosTheta).real
    
    return reflectanceP, reflectanceS, transmittanceP, transmittanceS
