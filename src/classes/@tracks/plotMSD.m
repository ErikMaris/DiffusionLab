function varargout = plotMSD(obj, ha, errorbarFlag, indices)
% PLOTMSD Plot the mean square displacement curves.
% 
% 
% obj(ii).plotMSD plots the MSD curves in the current axes.
%
% obj(ii).plotMSD(ha) plots the MSD curves in the axes
% specified by the handle ha.
%
% obj(ii).plotMSD(ha, errorbar), where errorbar is a
% boolean flag, allows to specify whether the curves should be
% plotted with error bars (equal to standard deviation). It is
% false by default.
%
% obj(ii).plotMSD(ha, errorbar, indices) plots the MSD curves for the
% particles with the specified indices only. Leave empty to
% plot MSD for all particles.
%
% hps =  obj(ii).plotMSD(...) returns the handle array for the
% lines generated.
%
% [hps, ha] =  obj(ii).plotMSD(...) also return the axes handle in
% which the lines were plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from msdanalyzer https://tinevez.github.io/msdanalyzer/
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


if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 3
    errorbarFlag = false;
end

all_indices = false;
if nargin < 4 || isempty(indices)
    all_indices = true;
end

nObj = numel(obj);

% preallocation hps
if all_indices
    nTracks = sum([obj.nTracks]);
else
    nTracks = nObj * numel(indices);
end
    
if errorbarFlag == 1
    hps  = gobjects(nTracks,4);
elseif strcmp(errorbarFlag,'errorbar')
    hps  = gobjects(nTracks,1);
else
    hps  = gobjects(nTracks,1);
end

trackcounter = 1;

for ii = 1:nObj
    
    if ~obj(ii).MSD_valid
        obj(ii) = obj(ii).computeMSD;
    end
    
    if all_indices
        indices = 1 : obj(ii).nTracks;
    end

    [xUnitFactor,xUnitLabel] = obj(ii).getUnitFactor('dt');
    [yUnitFactor,yUnitLabel] = obj(ii).getUnitFactor('pixelsize.^2');

    nTracks = numel(indices);
    
    if nObj > 1
        colors = obj(ii).getColormap(nObj);
        colors = repmat(colors(ii,:),[nTracks 1]);
    else
        colors = obj(ii).getColormap(nTracks);
    end
    
    for jj = 1:nTracks
        
        index = indices(jj);

        msd_spot = obj(ii).MSD{index};
        if isempty( msd_spot ) % do not plot empty spots in MSD
            continue
        end

        trackName = sprintf('Track %d', index );

        t = msd_spot(:,1).*xUnitFactor;
        m = msd_spot(:,2).*yUnitFactor;
        if errorbarFlag == 1
            s = msd_spot(:,3).*yUnitFactor; % standard deviation
            hp = tracks.errorShade(ha, t, m, s, '.-', colors(jj,:), true);
            hps(trackcounter,:) = [hp.mainLine hp.edge(1) hp.edge(2) hp.patch];
            set( hps(trackcounter,:), 'DisplayName', trackName );
        elseif strcmp(errorbarFlag,'errorbar') %
            s = msd_spot(:,3).*yUnitFactor; % standard deviation
            hps(trackcounter) = errorbar(ha, t, m, s, '.-',...
                'Color', colors(ii,:), ...
                'LineWidth', obj(ii).lineWidth, ...
                'DisplayName', trackName );
            set( hps(trackcounter), 'DisplayName', trackName );
        else
            hps(trackcounter) = plot(ha, t, m, '.-', ...
                'Color', colors(jj,:), ...
                'LineWidth', obj(ii).lineWidth, ...
                'DisplayName', trackName );
            set( hps(trackcounter), 'DisplayName', trackName );
        end
        trackcounter = trackcounter + 1;
        hold(ha, 'on');
    end
end

xlabel(['Delay time (' xUnitLabel ')'])
ylabel(['MSD (' yUnitLabel ')'])

if nObj > 1
    % to get the legend for the populations right, we use the last plotted
    % object in the axes to make our legend.
    legendTxt = cell(nObj,1);
    if all_indices
        legendN = cumsum([obj.nTracks]);
    else
        legendN = (1:nObj)'.*numel(indices);
    end
    legendHps = zeros(nObj,1);
    for ii = 1:nObj
        legendTxt{ii} = ['Population ' num2str(ii) ];
        legendHps(ii) = hps(legendN(ii));
    end
    legend(legendHps,legendTxt{:})
end


if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end
end