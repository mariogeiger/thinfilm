mkoctfile -mex -v -o thinfilm.mex tf-octave.cc

angles = linspace(0, 0.8*pi/2, 100);
lambdas = linspace(200, 1000, 110);

[r,t,a] = thinfilm(angles, lambdas, pi/4, 1, 2, [200], [1.5]);

mesh(lambdas, angles, r);
