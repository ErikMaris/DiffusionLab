function obj = computeDrift(obj, method, extra, interpmethod)
% COMPUTEDRIFT Compute and store drift correction.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj = obj.computeDrift(method) computes and stores the drift
% using one of the 4 following methods:
%
% 'clear' does not compute drift and remove any prior drift
% computation results.
%
% 'manual' allow to specify manually the drift vector:
% obj = obj.computeDrift('manual', dv); where dv is a double array
% of size N x (nDim+1) (nDim being the problem dimensionality), and
% must be arranged as following: [ Ti Xi Yi ... ] etc...
% On top of this, the drift vector must cover all the possible time
% points specified in the tracks field of this object: It must
% start before the first point and end after the last one,
% otherwise an error is thrown.
%
% Missing values within these extremities are interpolated using a
% linear interpolation scheme by default. To specify another
% interpolation scheme, use the following syntax:
% obj.computeDrift('manual', dv, interpmethod), with interpmethod
% being any value accepted by interp1 ('linear', 'nearest',
% 'spline', 'pchip', 'cubic').
%
% 'centroid' derives the drift by computing the center of mass of
% all particles at each time point. This method best work for a
% large number of particle and when the same number of particles is
% found at every time points. It fails silently otherwise.
%
% 'velocity' derives drift by computing instantaneous velocities
% and averaging them all together, at each time point. If the
% particles are in sufficient number, and if their real movement is
% indeed uncorrelated, the uncorrelated part will cancel when
% averaging, leaving only the correlated part. We assume this part
% is due to the drift. This method is more robust than the
% 'centroid' method against particle disappearance and appearance.
%
% 'movmean' derives the drift from the moving mean of a number of specified
% tracks. obj = obj.computeDrift('movmean',track_indices,k) with track_
% indices, the indices of the tracks taken for the calculation of drift,
% and k the window of the moving mean. The tracks used for this analysis
% should have a fixed position over time. In order to reduce the random
% error in the localization, a moving mean is taken. Since the error is
% uncorrelated, k of 10 usually sufficient, but can be taken larger.
%
% Results are stored in the 'drift' field of the returned object.
% It is a double array of size N x (nDim+1) (nDim being the problem
% dimensionality), and must be arranged as following: [ Ti Xi Yi ... ]
% etc. If present, it will by used for any call to computeMSD andother
% methods.
% 

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Erik Maris 05/2019
% Added moving mean drift correction 
% 
% Adapted from msdanalyzer (Erik Maris 2019/03/29)
% https://tinevez.github.io/msdanalyzer/
% Copyright (C) 2013 - 2014 Jean-Yves Tinevez
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
    n_tracks = obj(ii).nTracks;

    switch lower(method)

        case 'manual'

            if nargin < 4
                interpmethod = 'linear';
            end

            drift_dim = size(extra, 2);
            if drift_dim ~= obj(ii).nDim + 1
                error('tracks:computeDrift:BadDimensionality', ...
                    'Drift must be of size N x (nDim+1) with [ T0 X0 Y0 ... ] etc...');
            end

            uninterpolated_time = extra(:, 1);
            uninterpolated_drift = extra(:, 2:end);
            if min(uninterpolated_time) > min(time) || max(uninterpolated_time) < max(time)
                error('msdanalyzer:computeDrift:BadTimeVector', ...
                    'For manual drift correction, time vector must cover all time vector from all tracks.');
            end

            ldrift = interp1(...
                uninterpolated_time, ...
                uninterpolated_drift, ...
                time, interpmethod); % interpolates using 'interpmethod' between
                                    % 'uninterpolated_time' and
                                    % 'uninterpolated_drift' for values 'time'

            obj(ii).drift = [time ldrift];

        case 'centroid'

            ldrift = zeros(n_times, obj(ii).nDim);
            n_drift = zeros(n_times, 1);
            for jj = 1 : n_tracks

                % round track times to value set in getCommonTimes
                t = obj(ii).time{jj};
                t = tracks.roundn(t, tracks.TOLERANCE);

                % Determine target time index in bulk
                % returns indices(time) intections time and t
                [~, index_in_all_tracks_time, ~] = intersect(time, t);

                % Add to mean accum for these indexes
                n_drift(index_in_all_tracks_time) = n_drift(index_in_all_tracks_time) + 1; % count additions per time bin
                ldrift(index_in_all_tracks_time, :) = ldrift(index_in_all_tracks_time, :) + obj(ii).coords{jj}; % add coords to time bin

            end

            ldrift = ldrift ./ repmat(n_drift, [1 obj(ii).nDim]); % get mean by devision by number of particles in time bin
            obj(ii).drift = [time ldrift];


        case 'velocity'

            sum_V = zeros(n_times, obj(ii).nDim);
            n_V = zeros(n_times, 1);

            for jj = 1 : n_tracks

                % round track times to value set in getCommonTimes
                t = obj(ii).time{jj};
                t = tracks.roundn(t, tracks.TOLERANCE);

                % Determine target time index in bulk
                % returns indices(time) intections time and t
                [~, index_in_all_tracks_time, ~] = intersect(time, t);

                % Remove first element
                index_in_all_tracks_time(1) = [];

                % Compute speed
                V = diff( obj(ii).coords{jj} ) ./ repmat(diff(t), [ 1 obj(ii).nDim]);

                % Add to mean accum for these indexes
                n_V(index_in_all_tracks_time) = n_V(index_in_all_tracks_time) + 1; % count additions per time bin
                sum_V(index_in_all_tracks_time, :) = sum_V(index_in_all_tracks_time, :) + V; % add velocity to time bin

            end

            % Build accumulated drift
            sum_V(1, :) = 0;
            n_V(1, :) = 1;
            % Integrate
            d_time = [0; diff(time) ];
            % average dx/dt per time * dt = dx; cumsum takes the sum of all deplacements as the drift vector 
            ldrift = cumsum( sum_V ./ repmat(n_V, [1 obj(ii).nDim]) .* repmat(d_time, [1 obj(ii).nDim]), 1); 
            obj(ii).drift = [time ldrift];

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
            error('msdanalyzer:computeDriftCorrection:UnknownCorrectionMethod', ...
                'Unknown correction method %s. Must be ''clear'', ''manual'', ''centroid'' or ''velocity''.', ...
                method);
    end

    obj(ii).MSD_valid = false;
%     obj(ii).vcorr_valid = false;
    
end


end