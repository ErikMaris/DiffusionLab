function obj = computeTraceBlinkingGap(obj,blinkingGap)
% ONTRACEBLINKINGGAP fills the blinking gap in the on-trace
%   as specified by blinkingGap

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

nObj = length(obj);

for ii = 1:nObj
    for jj = 1:size(obj(ii).onTraceShort,1)
        % --- get indices that must be merged
        idx = find(obj(ii).onTraceShort{jj} <= blinkingGap); % find indices of positions that must be merged
        if isempty(idx); continue; end % otherwise kk = 1, which gives out of bounds
        idx = idx((mod(idx,2) == 1),1); % get off positions
        if isempty(idx); continue; end % otherwise kk = 1, which gives out of bounds
        if idx(end) == size(obj(ii).onTraceShort{jj},1); idx = idx(1:end-1); end % remove last value if it is the last one of the onTraceShort
        if isempty(idx); continue; end % otherwise kk = 1, which gives out of bounds
        if idx(1) == 1; idx = idx(2:end); end % remove first value in onTrace if present
        if isempty(idx); continue; end % otherwise kk = 1, which gives out of bounds
        % --- perform merging
        count = 0; % counts the number of changed positions
        for kk = 1:size(idx,1)
            rp = idx(kk)-count;
            obj(ii).onTraceShort{jj} = [obj(ii).onTraceShort{jj}(1:rp-2);...
                sum(obj(ii).onTraceShort{jj}(rp-1:rp+1)); obj(ii).onTraceShort{jj}(rp+2:end)]; % merge position rp with neighbours
            count = count + 2; % new positions have shifted -2 with respect new positions
        end
        % --- calculate onTrace from onTraceShort
        obj(ii).onTrace{jj} = getOnTraceFromShort(obj(ii),obj(ii).onTraceShort{jj}); % store onTrace from onTraceShort
    end
end