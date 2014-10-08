#include <mex.h>
#include "thinfilm.hh"

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
	
	thinfilm::complex nIncident(*mxGetPr(prhs[3]), mxIsComplex(prhs[3]) ? *mxGetPi(prhs[3]) : 0.0);
	thinfilm::complex nExit    (*mxGetPr(prhs[4]), mxIsComplex(prhs[4]) ? *mxGetPi(prhs[4]) : 0.0);
	
	const mxArray* thicknesses = prhs[5];
	const mxArray* indexes   = prhs[6];
	
	if (mxGetNumberOfElements(thicknesses) != mxGetNumberOfElements(indexes))
		mexErrMsgTxt("thicknesses and indexes must contain the same amount of elements");
	
	std::vector<thinfilm::Layer> layers;
	for (int i = 0; i < mxGetNumberOfElements(thicknesses); ++i) {
		thinfilm::complex index(*mxGetPr(indexes), mxIsComplex(indexes) ? *mxGetPi(indexes) : 0.0);
		if (imag(index) > 0.0) mexErrMsgTxt("the complex part of the index of refraction must be negative or null");
		thinfilm::Layer layer;
		layer.thickness = mxGetPr(thicknesses)[i];
		layer.refractiveIndex = index;
		layers.push_back(layer);
	}
	
	for (int j = 0; j < lambdaLength; ++j) {
		for (int i = 0; i < thetaLength; ++i) {
			thinfilm::complex incidentCosTheta = cos(mxGetPr(theta)[i]);

			int k = j * thetaLength + i;
			thinfilm::simulate(incidentCosTheta, mxGetPr(lambdaArray)[j], polarization, nIncident, nExit, layers, 
				reflectance   ? &reflectance[k]   : 0, 
				transmittance ? &transmittance[k] : 0,
				absorptance   ? &absorptance[k]   : 0, 
				psi           ? &psi[k]           : 0, 
				delta         ? &delta[k]         : 0);
		}
	}
}
