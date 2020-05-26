% This function is part of the code of 
%
% Sutour, C., Deledalle, C.A., Aujol, J.F. 
% Estimation of the noise level function based on a nonparametric detection
% of homogeneous image regions. 
%
% It is availible at  
% https://github.com/csutour/RNLF
%
function v = getoptions(options, name, v, mendatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end 