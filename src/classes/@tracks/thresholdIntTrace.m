function obj = thresholdIntTrace(obj,threshold)
% thresholdIntTrace Thresholds the value for the intensity trace.
%
% obj = thresholdIntTrace(obj,threshold) thresholds the value for the 
%   intensity trace into an onTrace with x < threshold = 0 and 
%   x >= threshold = 1.


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

for ii = 1:numel(obj)
    obj(ii).onTrace = cell(obj(ii).nTracks,1);
    obj(ii).onTraceShort = cell(obj(ii).nTracks,1);
    for jj = 1:obj(ii).nTracks
        obj(ii).onTrace{jj} = obj(ii).intTrace{jj} >= threshold;
        obj(ii).onTraceShort{jj} = getShortFromOnTrace(obj,obj(ii).onTrace{jj});
    end
    obj(ii).onTrace_valid = true;
end