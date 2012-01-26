## usage: thinfilm_spectrum need a array of lamda, n and k values

function rta = thinfilm_spectrum (theta, lamda, pol, ninc, kinc, nexit, kexit, varargin)

    w = length (lamda);
    
    for i = 1:w
        theta(i) = theta(1);
        pol(i) = pol(1);
        ninc
    endfor
    
    arrayfun(@thinfilm, [0,0], [550,600], [45,45], [1,1], [0,0], [1,1], [0,0], [30,30], [1.54,1.56], [0.001,0.001], "UniformOutput", false);
    
    for i = 1:w
    
        ## create layers for one wavelength
        clear layers;
        for j = 1:length (varargin)
            if mod (j, 3) == 1
                layers{j} = varargin{j};
            else
                layers{j} = varargin{j}(i);
            endif
        endfor

        ## run thinfilm for one wavelength       
        ret = thinfilm (theta, lamda(i), pol, ninc(i), kinc(i), nexit(i), kexit(i), layers)
        
        rta(i,1) = ret(1);
        rta(i,2) = ret(2);
        rta(i,3) = ret(3);
    endfor

endfunction

