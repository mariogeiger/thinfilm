# pylint: disable=C,R
import cmath
import numpy as np

# layers = [(index of refraxion, thickness / lambda in vaccum), second layer, ...]
def compute(incidentCosTheta, nIncident, nExit, layers):

    productMatrixP = np.eye(2)
    productMatrixS = np.eye(2)

    for layer in layers:
        n = layer[0]
        thickness = layer[1] # express in mutiple of the wavelength

        cosTheta = cmath.sqrt(1. - (1. - incidentCosTheta ** 2) * (nIncident / n) ** 2)

        deltaLayer = 2.0 * cmath.pi * n * thickness * cosTheta
        c = cmath.cos(deltaLayer)
        s = cmath.sin(deltaLayer)

        productMatrixP = np.dot(productMatrixP, np.array([[c, -1j * s * cosTheta / n], [-1j * s * n / cosTheta, c]]))
        productMatrixS = np.dot(productMatrixS, np.array([[c, -1j * s / cosTheta / n], [-1j * s * n * cosTheta, c]]))

    exitCosTheta = cmath.sqrt(1. - (1. - incidentCosTheta ** 2) * (nIncident / nExit) ** 2)

    bP, cP = np.dot(productMatrixP, [exitCosTheta, nExit]) # has been multiplied by exitCosTheta to avoid division by 0
    bS, cS = np.dot(productMatrixS, [1.0, nExit * exitCosTheta])

    numerator = 4 * abs(nExit * exitCosTheta * nIncident * incidentCosTheta)
    denomP = abs(bP * nIncident + cP * incidentCosTheta) ** 2
    denomS = abs(bS * nIncident * incidentCosTheta + cS) ** 2

    if denomP != 0.0:
        reflectanceP = abs(bP * nIncident - cP * incidentCosTheta) ** 2 / denomP
        transmittanceP = numerator / denomP
    else:
        reflectanceP = -1
        transmittanceP = -1

    if denomS != 0.0:
        reflectanceS = abs(bS * nIncident * incidentCosTheta - cS) ** 2 / denomS
        transmittanceS = numerator / denomS
    else:
        reflectanceS = -1
        transmittanceS = -1

    return reflectanceP, reflectanceS, transmittanceP, transmittanceS
