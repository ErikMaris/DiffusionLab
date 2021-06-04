function obj = computeOnTrace(obj,allFramesFlag,commonTimes)
% computeOnTrace Computes the continuous on-trace based on the track's time tags
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj = computeOnTraceC(obj,allFramesFlag,commonTimes)
%   allFramesFlag/boolean: if true, the on-off trace is padded with zeros 
%   for the movie length.
%   commonTimes: must be an evenly spaced vector with nTimes + 1 entries

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

if nargin < 2
    allFramesFlag = true;
end

if nargin > 2
    if any(diff(commonTimes) <= 0)
        error('SC_tracks:computeOnTrace:BadTimeVector', ...
            'Common times vector is not strictly increasing.');
    end
    commonTimes = commonTimes(:);
end


for ii = 1:numel(obj)
    if nargin < 3
        commonTimes_s = obj(ii).getCommonTimes;
        dt = tracks.roundn(min(diff(commonTimes_s)), tracks.TOLERANCE); % estimate dt
        commonTimes = [commonTimes_s; commonTimes_s(end) + dt];
%         if mod(commonTimes(end) - commonTimes(1),dt) > 10^-tracks.TOLERANCE
%             error('SC_tracks:computeOnTraceC:BadCommonTimes', ...
%                 'Common times vector cannot be evenly spaced.');
%         end
%         commonTimes = dt:dt:commonTimes(end)+dt; % fill gaps in common times (assumes equally spaced times starting at frame=1)
    else
        commonTimes_s = commonTimes(1:end-1);
    end
    obj(ii).onTrace = cell(obj(ii).nTracks,1);
    for jj = 1:obj(ii).nTracks
        onTrace = ismembertol(commonTimes_s(:), obj(ii).time{jj}(:),10^-tracks.TOLERANCE);
        N = numel(onTrace);
        if sum(onTrace) == 0 % if is all 0
            obj(ii).onTrace{jj} = N;
            continue
        end
        if sum(onTrace) == N % if is all one
            obj(ii).onTrace{jj} = [0; N];
            continue
        end
%         onTrace = onTrace(2:end);
        if ~allFramesFlag
            thisFrame = find(onTrace);
            onTrace = onTrace(thisFrame(1):thisFrame(end));
        end
        appendVal = [];
        if onTrace(1) == 1 % if tracks starts on
            appendVal = 0; % add a zero
        end
        First_idx =  find(abs(onTrace-circshift(onTrace,1)) > 10^-tracks.TOLERANCE); % find(DT-circshift(DT,1) ~= 0); first idx is dummy to mark first frame, which is at Idx 2 (Idx = 1 is reserved fro t-1 (t = 0s))
        Last_idx = find(abs(onTrace-circshift(onTrace,-1)) > 10^-tracks.TOLERANCE); % find(DT-circshift(DT,-1) ~= 0)
        if First_idx(1) ~= 1 % ugly fix to include first index
            First_idx = [1; First_idx];
        end

        if Last_idx(end) ~= N % ugly fix to include last index
            Last_idx = [Last_idx; N];
        end
        obj(ii).onTrace{jj} = [appendVal; (commonTimes(Last_idx+1)-commonTimes(First_idx))]; % assumes equally spaced times starting at 0 (= firstIdx-1:lastIdx)
    end
    obj(ii).onTrace_valid = true;
end


end