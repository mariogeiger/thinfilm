## usage: thinfilm theta(deg) lamda(nm) pol(deg0S) ninc kinc nexit kexit [ dlayer(nm) nlayer klayer ...]

function [reflectance,transmittance,absorptance] = thinfilm (varargin)

    ## write command
    cmd = "multilayer-thinfilm-octave";
    for i = 1:length (varargin)
        cmd = cstrcat(cmd, " ", num2str(varargin{i}));
    endfor
    
    ## run c++ program
    [status, ret] = unix(cmd);

    if (status != 0)
        error ("multilayer-thinfilm return error");
    endif
    
    ## parse output
    [reflectance,transmittance,absorptance] = sscanf (ret, "%f %f %f");

endfunction

