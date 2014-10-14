#include <mex.h>
#include "thinfilm.hh"

std::vector<thinfilm::complex> readComplexVector(const mxArray* array);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nrhs < 7)
		mexErrMsgTxt("input : incident angle, wavelength, polarization (0=P,pi/2=S), incident and exit complex ref. index, [thicknesses], [ref. indexes]");
	if (nlhs == 0)
		mexErrMsgTxt("output : power ratios : reflectivity, transmittance, absorptance");

	const mxArray* theta = prhs[0];
	int thetaLength = mxGetNumberOfElements(theta);
	
	const mxArray* lambdaArray = prhs[1];
	int lambdaLength = mxGetNumberOfElements(lambdaArray);
	
	double* reflectance = 0;
	if (nlhs >= 1) {
		plhs[0] = mxCreateDoubleMatrix(thetaLength,lambdaLength,mxREAL);
		reflectance = mxGetPr(plhs[0]);
	}
	double* transmittance = 0;
	if (nlhs >= 2) {
		plhs[1] = mxCreateDoubleMatrix(thetaLength,lambdaLength,mxREAL);
		transmittance = mxGetPr(plhs[1]);
	}
	double* absorptance = 0;
	if (nlhs >= 3) {
		plhs[2] = mxCreateDoubleMatrix(thetaLength,lambdaLength,mxREAL);
		absorptance = mxGetPr(plhs[2]);
	}
	double* psi = 0;
	double* delta = 0;
	if (nlhs >= 4) {
		mexErrMsgTxt("only 3 outputs are supported");
	}
	
	double polarization = mxGetScalar(prhs[2]);
	
	std::vector<thinfilm::complex> nIncident = readComplexVector(prhs[3]);
	if (nIncident.size() != 1 && nIncident.size() != lambdaLength)
		mexErrMsgTxt("the incident index of refraction must either be a complex scalar or a vector of same length of lambda");

	std::vector<thinfilm::complex> nExit     = readComplexVector(prhs[4]);
	if (nExit.size() != 1 && nExit.size() != lambdaLength)
		mexErrMsgTxt("the exit index of refraction must either be a complex scalar or a vector of same length of lambda");
	
	const mxArray* thicknesses = prhs[5];
	const mxArray* indexes     = prhs[6];
	
	if (mxGetNumberOfElements(thicknesses) != mxGetNumberOfElements(indexes))
		mexErrMsgTxt("thicknesses and indexes must contain the same amount of elements");
	
	int nLayers = mxGetNumberOfElements(thicknesses);
	std::vector<thinfilm::Layer> layers;
	for (int i = 0; i < nLayers; ++i) {
		thinfilm::complex index(*mxGetPr(indexes), mxIsComplex(indexes) ? *mxGetPi(indexes) : 0.0);
		if (imag(index) > 0.0) mexErrMsgTxt("the complex part of the index of refraction must be negative or null");
		thinfilm::Layer layer;
		layer.thickness = mxGetPr(thicknesses)[i];
		layer.refractiveIndex = index;
		layers.push_back(layer);
	}
	
	for (int j = 0; j < lambdaLength; ++j) {
		thinfilm::complex nInc = nIncident.size() == lambdaLength ? nIncident[j] : nIncident[0];
		thinfilm::complex nExi = nExit.size() == lambdaLength ? nExit[j] : nExit[0];
		for (int i = 0; i < thetaLength; ++i) {
			thinfilm::complex incidentCosTheta = cos(mxGetPr(theta)[i]);

			int k = j * thetaLength + i;
			thinfilm::simulate(incidentCosTheta, mxGetPr(lambdaArray)[j], polarization, nInc, nExi, layers, 
				reflectance   ? &reflectance[k]   : 0, 
				transmittance ? &transmittance[k] : 0,
				absorptance   ? &absorptance[k]   : 0, 
				psi           ? &psi[k]           : 0, 
				delta         ? &delta[k]         : 0);
		}
	}
}

std::vector<thinfilm::complex> readComplexVector(const mxArray* array)
{
	std::vector<thinfilm::complex> result;
	int n = mxGetNumberOfElements(array);
	
	double* re = mxGetPr(array);
	double* im = mxGetPi(array);
	if (mxIsComplex(array)) {
		for (int i = 0; i < n; ++i) {
			result.push_back(thinfilm::complex(re[i], im[i]));
		}
	} else {
		for (int i = 0; i < n; ++i) {
			result.push_back(re[i]);
		}
	}
	return result;
}
