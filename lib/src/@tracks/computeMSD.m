function obj = computeMSD(obj)
% computeMSD Measures the mean squared displacement from a trajectories
% 
% Computation is vectorized!
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% 

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from Maxime Deforet, May 21 2013
% https://nl.mathworks.com/matlabcentral/fileexchange/
% 41858-kehl-a-fast-no-loop-method-to-compute-msd
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

%     delayIdx = obj(ii).dtSDidx;
    trackIndices = 1 : obj(ii).nTracks;

    nTracks = numel(trackIndices);
    fprintf('Computing MSD of %d tracks... \n', nTracks);

    obj(ii).MSD = cell(nTracks, 1); % preallocation
%     obj(ii).dtSD = cell(nTracks, 1); % preallocation
    if ~isempty(obj(ii).drift)
        tdrift = obj(ii).drift(:,1);
        xdrift = obj(ii).drift(:, 2:end);
    end
    
    dispVal = round(nTracks/10);
    
    % use getdR2 now
%     allDelays = obj(ii).getAllDelays;
%     
%     delayIdx = delayIdx + 1; % skip zero (= index 1)
%     if max(delayIdx+1) > max(allDelays)
%         delayIdx = delayIdx(delayIdx <= max(allDelays)); % only take valid indices
%     end
% 
%     takeDelay = allDelays(delayIdx); 
    
    for jj = 1 : nTracks
        
        if mod(jj,dispVal) == 0; fprintf('%d%% ',jj/dispVal*10); end % display

        index = trackIndices(jj);
        
        t = obj(ii).time{index};
        t = tracks.roundn(t, tracks.TOLERANCE);
        X = obj(ii).coords{index};

        % Determine drift correction
        if ~isempty(obj(ii).drift)
            % Determine target delay index in bulk
            % find to-be-corrected time indices in track
            [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
            % Keep only track times that can be corrected.
            X = X(index_in_track_time, :);
            t = t(index_in_track_time);
            % Subtract drift position to track position
            X = X - xdrift(index_in_drift_time, :);
        end

        n_detections = size(X, 1);
        
        % This code contains no loop. Each step is "vectorized", therefore pretty
        % fast.
        % The idea is to compute all the possible pairings in only one time.
        % Trajectory = list of T positions (x,y,t). Or (x,t) or (x,y,z,t) or even
        % higher dimension (x1,x2,x3,...t)
        % tau is the list of all the possible time intervals within the trajectory.
        % MSD is a list of mean squared displacements, for each value of tau.
        % MSD(tau) = sqrt( x(i)^2 - x(i+tau)^2 ), averaged over i.
        
        [I, j] = find(triu(ones(n_detections),1)); % list of indices of possible pairings: I row; j col
        D = zeros(n_detections, n_detections);
        D(I + n_detections*(j-1)) = (sum(abs(X(I,:) - X(j,:)).^2,2)); % Squared distance computation in one line !
        % Time intervals between the two points of each pairing :
        dt = zeros(n_detections, n_detections);
        dt(I + n_detections*(j-1)) = -(t(I) - t(j));
        % Then D is a list of squared distances. dt is a list of corresponding
        % time intervals. Now we have to average all D values over different dt
        % values
        % We first remove all 0 values from dt matrix, and from D as well.
        idx_0 = find(dt == 0);
        dt(idx_0) = [];
        D(idx_0) = [];
        % Then we sort dt in ascending order, and sort D in the same way.
        [DT,idx] = sort(dt(:));
        DT = tracks.roundn(DT, tracks.TOLERANCE); % this is important, otherwise identical delays will not be recognised
        DD = D(idx);
        % We now have DD and DT, respectively a list of squared distances, and
        % the corresponding time intervals.
        % Now we have to compute the mean DD for each possible value of DT.
        % Let's get the first and last indices of each set of DT
        if isempty(DD) % 1 track point
            continue
        elseif numel(DD) == 1 % 2 track points
            delays = DT;
            mean_msd = DD;
            n_msd = 1;
            std_msd = 0;
        else
            % tolerance in this step is not required to due rounding with tolerance of DT
            %    First_idx =  find(abs(DT-circshift(DT,1)) > 10^-tracks.TOLERANCE); 
            %    Last_idx = find(abs(DT-circshift(DT,-1)) > 10^-tracks.TOLERANCE);
            First_idx = find(DT-circshift(DT,1) ~= 0);
            Last_idx = find(DT-circshift(DT,-1) ~= 0);
            
            % For instance, DT=1 start at First_idx(1) and end at Last_idx(1)
            %               DT=2 start at First_idx(2) and end at Last_idx(2)...
            % To get the average, we first compute the cumulative (starting from 0), then we
            % get "the derivative of the cumulative".
            C = cumsum([0,DD]);
            % For each possible value of DT, the mean of DD is the difference between
            % the initial index and the last index of the cumulative, divided by the
            % number of elements between those two indices :
            mean_msd = (C(Last_idx+1)-C(First_idx))'./(Last_idx-First_idx+1); 
            delays = DT(First_idx); % list of intervals

            % based on the formula for the sample standard deviation
            % (Wikipedia) and tested against std() and msdanalyzer
            n_msd = (Last_idx-First_idx+1);
            CD = cumsum([0, (DD - mean_msd(repelem(1:numel(delays),n_msd))').^2]);
            std_msd = sqrt((CD(Last_idx+1)-CD(First_idx))'./(n_msd-1)); 
    %         std_msd = sqrt((CD(Last_idx+1)-CD(First_idx))'./(n_msd-1))./sqrt(n_msd);  % do standard error of the mean
            std_msd(isnan(std_msd)) = 0;

            % We replace points for which N=0 by Nan, to later treat
            % then as missing data. Indeed, for each msd cell, all the
            % delays are present. But some tracks might not have all
            % delays
            delay_not_present = n_msd == 0;
            mean_msd( delay_not_present ) = NaN;
        end
        
        obj(ii).MSD{index} = [0 0 0 n_detections; delays mean_msd std_msd n_msd ]; % add first line manually
        
        % use getdR2 for this functionality
%         takeSD = ismembertol(DT,takeDelay,10^-tracks.TOLERANCE);
%         obj(ii).dtSD{index} = [DT(takeSD) DD(takeSD)'];
        
    end
    fprintf('\n')

    obj(ii).MSD_valid = true;
    
end

end

