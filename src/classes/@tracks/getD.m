function varargout = getD(obj,type)
% getD Gets diffusion constant for diffusion estimator
% 
% getD extracts the diffusion constant for the Dest object and converts it
% to the current units.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% varargout = getD(obj,type)
%   type / char: 'tracks' outputs per track, 'mean' for the mean of the 
%    population, and 'bootstrap' for the bootstrapped population.

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
    error('tracks:getD:multipleObjs',...
        '%d objects parsed, can only process one.',numel(obj))
end

if nargin < 2
    type = 'tracks';
end

[uF,uL] = obj.getUnitFactor('pixelsize.^2/dt');
switch type
    case 'tracks'
        if isfield(obj.Dest.derivedVals,'D')
            out = obj.Dest.derivedVals.D.*uF;
        else
            error('tracks:getD:noD',...
                'Estimator does not compute diffusion constant.')
        end
    case 'mean'
        if isfield(obj.Dest.derivedVals,'mD')
            out = obj.Dest.derivedVals.mD*uF;
        else
            error('tracks:getD:nomD',...
                'Estimator does not compute mean diffusion constant.')
        end
    case 'bootstrap'
        if isfield(obj.Dest.derivedVals,'bD')
            out = obj.Dest.derivedVals.bD*uF;
        else
            error('tracks:getD:nobD',...
                'Estimator does not compute bootstrapped diffusion constant.')
        end
    otherwise
        error('tracks:getD:unrecognizedType',...
            'Type %s is not recognized by function %s', type,mfilename)
end

varargout{1} = out;
if nargout > 1
    varargout{2} = uL;
end
    

end

