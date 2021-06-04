function overlapPairs = getTrackOverlap(obj,spatialThres,temporalThres)
% getTrackOverlap Finds the tracks within a spatial and temporal distance
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% overlapPairs = getTrackOverlap(obj,spatialThres,temporalThres)
%   spatialThres: minimal distance (before unit conversion) between tracks
%   to be overlap computed from the minimum bouding circle
%   temporalThres: minimal distance (before unit conversion) between tracks
%   to be overlap computed from the begin and end of the track

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

if ~obj.center_valid
    obj = obj.computeCenter;
end

overlapPairs = [];
R = obj.features.MinBoundCircleRadius;
CC = obj.features.MinBoundCircleCenter;
for ii = 1:obj.nTracks
    for jj = ii+1:obj.nTracks % loop over all pairs
        % --- check for overlap spatialThres + minboundcircle 
        r_tot = R(ii) + R(jj);
        r_CenterCenter = sqrt((CC(ii,1) - CC(jj,1))^2 + ...
            (CC(ii,2) - CC(jj,2))^2);
        if r_tot + spatialThres > r_CenterCenter
            % --- check for overlap in time + temporalThres
            if obj.time{ii}(end) + temporalThres > obj.time{jj}(1) || obj.time{jj}(end) + temporalThres > obj.time{ii}(1)
                overlapPairs = vertcat(overlapPairs,[ii jj]);
            end
        end
    end
end


end

