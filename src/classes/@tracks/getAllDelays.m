function delays = getAllDelays(obj, indices)
% getAllDelays Computes all unique delays in all tracks.
% 
% Delays have been rounded to tolerance.
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% delays = getAllDelays(obj, indices)
%   indices: indices of tracks considered in the function

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from msdanalyzer https://tinevez.github.io/msdanalyzer/
% (C) 2007 Jean-Yves Tinevez
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

if nargin < 2 || isempty(indices)
    indices = 1 : obj.nTracks;
end

nTracks = numel(indices);
all_delays = cell(nTracks,1);
for ii = 1 : nTracks
    % First, find all possible delays in time vectors.
    % Time can be arbitrary spaced, with frames missings,
    % non-uniform sampling, etc...
    index = indices(ii);
    t = obj.time{index};
    [T1, T2] = meshgrid(t, t); % make raster with times along x and y
    dT = tracks.roundn(abs(T1(:)-T2(:)), tracks.TOLERANCE); % calculate all possible time differences
    all_delays{ii} = unique(dT); % remove duplicates
end
delays = unique( vertcat(all_delays{:}) ); % merge and remove duplicates
end