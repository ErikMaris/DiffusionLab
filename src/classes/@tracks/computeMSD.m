function obj = computeMSD(obj)
% computeMSD Measures the mean squared displacement from a trajectories
%
% obj = obj.computeMSD returns an object with the mean squared displacement
% stored in the property 'MSD' in tracks.


% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
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

    trackIndices = 1 : obj(ii).nTracks;

    nTracks = numel(trackIndices);
    fprintf('Computing MSD of %d tracks... \n', nTracks);

    obj(ii).MSD = cell(nTracks, 1); % preallocation
    
    % get drift matrix
    if ~isempty(obj(ii).drift)
        tdrift = obj(ii).drift(:,1);
        xdrift = obj(ii).drift(:, 2:end);
    end
    
    dispVal = round(nTracks/10);
    
    for jj = 1 : nTracks
        
        % display
        if mod(jj,dispVal) == 0; fprintf('%d%% ',jj/dispVal*10); end 

        % get current index
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
        
        % do MSD computation
        MSD = tracks.static_getMSD(t,X); % does MSD computation
        
        % write to object
        if ~isempty(MSD)
            obj(ii).MSD{index} = MSD;
        end
        
    end
    fprintf('\n')

    obj(ii).MSD_valid = true;
    
end

end

