function varargout = plotTracks(obj, ha, indices, driftFlag)
% PLOTTRACKS Plot the tracks stored in this object.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj(ii).plotTracks plots the particle trajectories stored in the
% tracks object obj in the current axes. This method only
% works for 2D or 3D problems.
%
% obj(ii).plotTracks(ha) plots the trajectories in the axes with
% handle ha.
%
% obj(ii).plotTracks(ha, indices) where indices is a vector, allows
% to specify the track to be plotted, using their indices.
% Leave the vector empty to plot all trajectories.
%
% obj(ii).plotTracks(ha, indices, driftFlag) where driftFlag is a
% boolean flag, allows to specify whether the plot should
% display the trajectories corrected for drift (true) or
% uncorrected (false). A proper drift vector must be computed
% prior to setting this flag to true. See tracks.computeDrift.
%
% hps = obj(ii).plotTracks(...) returns the handles to the
% individual line objects created.
%
% [hps, ha] = obj(ii).plotTracks(...) returns also the handle to
% the axes handle the trajectories are plot in.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted fromn msdanalyzer https://tinevez.github.io/msdanalyzer/
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

if any([obj.nDim] < 2) || any([obj.nDim] > 3)
    error('tracks:plotTracks:UnsupportedDimensionality', ...
        'Can only plot tracks for 2D or 3D problems, got %dD.', obj.nDim);
end

if numel(unique([obj.nDim])) > 1
    error('tracks:plotTracks:NonequalDimensionality', ...
        'Can only plot tracks with an identical dimensionality.');
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 4 || isempty(driftFlag)
    driftFlag = false;
end



nObj = numel(obj);
hps  = [];
nTracks = zeros(nObj,1);

for ii = 1:nObj
    
    if nargin < 3 || isempty(indices)
        indices = 1 : obj(ii).nTracks;
    end
    
    [unitFactor,unitLabel] = obj(ii).getUnitFactor('pixelsize');
    
    nTracks(ii) = numel(indices);
    
    if nObj > 1
        colors = obj(ii).getColormap(nObj);
        colors = repmat(colors(ii,:),[nTracks(ii) 1]);
    else
        colors = obj(ii).getColormap(nTracks(ii));
    end
    
    hp = NaN(nTracks(ii), 1);
    
    if obj(ii).nDim == 2
        % 2D case
        for jj = 1 : nTracks(ii)

            index = indices(jj);
            track = obj(ii).coords{index};
            trackName = sprintf('Track %d', index );

            if isempty(track)
                continue
            end

            x = track(:,1).*unitFactor;
            y = track(:,2).*unitFactor;

            if driftFlag && ~isempty(obj(ii).drift)
                tdrift = obj(ii).drift(:,1);
                xdrift = obj(ii).drift(:,2);
                ydrift = obj(ii).drift(:,3);
                t = track(:,1);
                [~, index_in_drift_time, ~] = intersect(tdrift, t);
                % Subtract drift position to track position
                x = x - xdrift(index_in_drift_time);
                y = y - ydrift(index_in_drift_time);
            end

            hp(jj) =  plot(ha, x, y, ...
                'Color', colors(jj,:), ...
                'DisplayName', trackName,...
                'LineWidth',obj(ii).lineWidth);
            
            hold(ha, 'on');

        end

        xlabel(['x (' unitLabel ')'])
        ylabel(['y (' unitLabel ')'])
        
        axis image
    else
        % 3D case
        [unitFactor,unitLabel] = obj(ii).getUnitFactor('pixelsize');
        for jj = 1 : nTracks(ii)

            index = indices(jj);
            track = obj(ii).coords{index};
            trackName = sprintf('Track %d', index );

            if isempty(track)
                continue
            end

            x = track(:,1).*unitFactor;
            y = track(:,2).*unitFactor;
            z = track(:,3).*unitFactor;

            if driftFlag && ~isempty(obj(ii).drift)
                tdrift = obj(ii).drift(:,1);
                xdrift = obj(ii).drift(:, 2);
                ydrift = obj(ii).drift(:, 3);
                zdrift = obj(ii).drift(:, 4);
                t = track(:,1);
                [~, index_in_drift_time, ~] = intersect(tdrift, t);
                % Subtract drift position to track position
                x = x - xdrift(index_in_drift_time);
                y = y - ydrift(index_in_drift_time);
                z = z - zdrift(index_in_drift_time);
            end

            hp(jj) =  plot3(ha, x, y, z, ...
                'Color', colors(jj,:), ...
                'DisplayName', trackName, ...
                'LineWidth', obj(ii).lineWidth);
            
            hold(ha, 'on');
        end

        xlabel(['x (' unitLabel ')'])
        ylabel(['y (' unitLabel ')'])
        zlabel(['z (' unitLabel ')'])

    end
    hps = [hps;hp];
    axis equal
end

xlim auto
ylim auto

if nObj > 1
    % to get the legend for the populations right, we use the last plotted
    % object in the axes to make our legend.
    legendTxt = cell(nObj,1);
    legendN = cumsum(nTracks);
    legendHps = zeros(nObj,1);
    for ii = 1:nObj
        legendTxt{ii} = ['Population ' num2str(ii) ];
        legendHps(ii) = hps(legendN(ii));
    end
    legend(legendHps,legendTxt{:})
end

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end