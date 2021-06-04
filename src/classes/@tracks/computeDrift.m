function obj = computeDrift(obj, method, extra, interpmethod)
% COMPUTEDRIFT Compute and store drift correction.
% 
% obj = obj.computeDrift(method) computes and stores the drift matrix using
% one of the following methods:
%
% 'clear' clears the stored drift matrix.
%
% 'movmean' derives the drift from the moving mean of a number of specified
% tracks. obj = obj.computeDrift('movmean',track_indices,k) with track_
% indices, the indices of the tracks taken for the calculation of drift,
% and k the window of the moving mean. The tracks used for this analysis
% should have a fixed position over time. In order to reduce the random
% error in the localization, a moving mean is taken. Since the error is
% uncorrelated, k of 10 usually sufficient, but can be taken larger.
%
% The drift matrix is stored in the field 'drift' in tracks. It is an
% matrix with dimensions [nTracks nDim+1] with columns [time coord1 coord2
% ...]. If present, it will be used by plotTracks and computeMSD.

% To do: apply drift correction in function computeFeatures.
%
% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Erik Maris 05/2019
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

    % First compute common time points
    time = obj(ii).getCommonTimes();
    n_times = numel(time);

    switch lower(method)

        case 'clear'

            obj(ii).drift = [];

        case 'movmean'
            
            track_indices = extra;
            k = interpmethod;
            
            if isempty(track_indices)
                
            end
            
            ldrift = zeros(n_times, obj(ii).nDim);
            n_drift = zeros(n_times, 1);
            for jj = track_indices

                % round track times to value set in getCommonTimes
                t = obj(ii).time{jj};
                t = tracks.roundn(t, tracks.TOLERANCE);

                % Determine target time index in bulk
                % returns indices(time) intections time and t
                [~, index_in_all_tracks_time, ~] = intersect(time, t);

                % calculate moving mean w.r.t. center of mass
                mm = movmean(obj(ii).coords{jj},k,1) - mean(obj(ii).coords{jj},1);
                
                % Add to mean accum for these indexes
                n_drift(index_in_all_tracks_time) = n_drift(index_in_all_tracks_time) + 1; % count additions per time bin
                ldrift(index_in_all_tracks_time, :) = ldrift(index_in_all_tracks_time, :) + mm ; % add coords to time bin drift

            end

            ldrift = ldrift ./ repmat(n_drift, [1 obj(ii).nDim]); % get mean by devision by number of particles in time bin
            obj(ii).drift = [time ldrift];


        otherwise
            error('Unknown drift correction method %s. Must be ''clear'' or ''movmean''.', method);
    end

    obj(ii).MSD_valid = false;
    
end


end