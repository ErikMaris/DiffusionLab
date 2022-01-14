function varargout = getLocVar(obj,type)
% getLocVar Gets diffusion constant for diffusion estimator
% 
% getLocVar extracts the localization variance for the Dest object and 
% converts it to the current units.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% varargout = getLocVar(obj,type)
%   type / char: 'tracks' outputs per track, 'mean' for the mean of the 
%    and population for the bootstrapped population.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% 
% 
% -------------------------------------------------------------
% Copyright (C) 2021 J.J. Erik Maris
% Inorganic Chemistry & Catalysis, Utrecht University
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.


if numel(obj) > 1
    error('tracks:getLocVar:multipleObjs',...
        '%d objects parsed, can only process one.',numel(obj))
end

if nargin < 2
    type = 'tracks';
end

[uF,uL] = obj.getUnitFactor('pixelsize.^2');

switch type
    case 'tracks'
        if isfield(obj.Dest.results,'locVar')
            out = obj.Dest.results.locVar.*uF;
        else
            error('tracks:getLocVar:noLocVar',...
                'Estimator does not compute localization variance.')
        end
    case 'mean'
        if isfield(obj.Dest.mResults,'locVar')
            out = obj.Dest.mResults.locVar*uF;
        else
            error('tracks:getLocVar:nomLocVar',...
                'Estimator does not compute mean localization variance.')
        end
    otherwise
        error('tracks:getLocVar:unrecognizedType',...
            'Type %s is not recognized by function %s', type,mfilename)
end

varargout{1} = out;
if nargout > 1
    varargout{2} = uL;
end
    

end

