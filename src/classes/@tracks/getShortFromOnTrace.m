function onTraceShort = getShortFromOnTrace(obj,onTrace)
% onTraceShort = getShortFromOnTrace(onTrace)
%--------------------------------------------------------------------------
% Description
%
% generates an onTraceShort from an onTrace
%
% Literature:
%
% 
%--------------------------------------------------------------------------
% Necessary Inputs (name/data type)
% 
% onTrace/vector: vector with the number of frames, in which for each
%   frame is indicated whether the trace is off (0) or on (1). Must be a
%   column vector.
%
%--------------------------------------------------------------------------
% Outputs (name/data type)
%
% onTraceShort/vector: odd values are off, even are on. Use a 0 as first
%   odd value if trace starts with on.
%
%--------------------------------------------------------------------------
% Variable Inputs (flag/ data type /(default)):
% 
% 
%--------------------------------------------------------------------------
% Dependencies
%
%
%--------------------------------------------------------------------------

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
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

appendVal = [];
if onTrace(1) == 1 % if tracks starts on
    appendVal = 0; % add a zero
end
diffIdx = find(diff(onTrace));
onTraceShort = [appendVal; ([diffIdx; numel(onTrace)] - [0; diffIdx])]; % all shorthands: odd is off, even is on

end

