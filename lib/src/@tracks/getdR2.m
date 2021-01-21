function [dR2,delayTimes] = getdR2(obj,delayIdx,indices)
% getdR2 Computes the squared displacement from trajectories
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


ii = 1;
if numel(obj) > 1
    error('Can only process one object. Got %i.', numel(obj))
end

if nargin < 3 || isempty(indices)
    indices = 1 : obj(ii).nTracks;
end

% get number of tracks
nTracks = numel(indices);

% let's provide some user info
fprintf('Computing MSD of %d tracks... \n', nTracks);
dispVal = round(nTracks/10);

% get decays to output
delayIdx = delayIdx + 1; % skip zero (= index 1)
allDelays = obj(ii).getAllDelays;
if max(delayIdx+1) > numel(allDelays)
    delayIdx = delayIdx(delayIdx <= numel(allDelays)); % only take valid indices
end
delayTimes = allDelays(delayIdx);

% get drift vectors
if ~isempty(obj(ii).drift)
    tdrift = obj(ii).drift(:,1);
    xdrift = obj(ii).drift(:, 2:end);
end

% preallocation
dR2_track = cell(nTracks,1);

% core is the same as MSD computation, but then stripped down
for jj = 1:nTracks

    % display progress
    if mod(jj,dispVal) == 0; fprintf('%d%% ',jj/dispVal*10); end

    % get current index
    index = indices(jj);

    % get values for current track
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
    DD = D(idx);
    % We now have DD and DT, respectively a list of squared distances, and
    % the corresponding time intervals.
    
    % store 
    DT = tracks.roundn(DT, tracks.TOLERANCE); % this is important, otherwise identical delays will not be recognised
%     take_idx = ismembertol(DT,delayTimes,10^-tracks.TOLERANCE);
    take_idx = ismembertol(DT,delayTimes);
    dR2_track{index} = [DT(take_idx) DD(take_idx)'];
end
fprintf('\n')

% get a matrix with all [dt MSD]
dR2_all = vertcat(dR2_track{:});

% preallocation
dR2 = cell(numel(delayTimes),1);

% get only desired delay times
for ii = 1:numel(delayTimes)
    dR2{ii} = sort(dR2_all(dR2_all(:,1) == delayTimes(ii),2));
end


end

