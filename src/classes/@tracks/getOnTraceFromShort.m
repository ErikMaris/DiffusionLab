function onTrace = getOnTraceFromShort(obj,onTraceShort)
% generates an onTrace from an onTraceShort
% input: 
%   onTraceShort/vector: odd values are off, even are on. Use a 0 as first
%   odd value if trace starts with on
% output:
%   onTrace/vector: vector with the number of frames, in which for each
%   frame is indicated whether the trace is off (0) or on (1).

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

numFrames = sum(onTraceShort);

if numFrames == 0
    error(['The variable onTraceShort does not contain any time steps empty in function ''' mfilename '''.'])
end

onTrace = zeros(numFrames,1);
if length(onTraceShort) == 1 % no on states (1); only zeros is OK
    return    
end

currPos = onTraceShort(1) + 1; % indicates the current index in the onTrace
for ii = 2:2:length(onTraceShort)-2
    onTrace(currPos:currPos+onTraceShort(ii)-1) = 1;
    currPos = currPos + onTraceShort(ii+1) + onTraceShort(ii); % adds the on (onTraceShort(ii)) and off (onTraceShort(ii-1)) values to the current position
end
if length(onTraceShort) == 3 % exception when length(onTraceShort) = 3
    ii = 0;
end
onTrace(currPos:currPos+onTraceShort(ii+2)-1) = 1; % ii+2 situation outside loop to prevent out of bounds of currPos + onTraceShort(ii+1) + onTraceShort(ii)



end

