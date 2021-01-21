function varargout = plotMeanMSD(obj, ha, errorbarFlag, indices)
% PLOTMEANMSD Plot the weighted mean of the MSD curves.
% 
% obj,plotMeanMSD computes and plots the weighted of all MSD
% curves. See tracks.getMeanMSD.
%
% More on the errorbar
% https://tinevez.github.io/msdanalyzer/tutorial/MSDTuto_brownian.html
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj,plotMeanMSD(ha) plots the curve in the axes with the
% specified handle.
%
% obj,plotMeanMSD(ha, errorbar) where 'errorbar' is a boolean allow
% to specify whether to plot the curve with error bars indicating
% the weighted standard deviation. Default is false. Set 'errorbar' to
% 'errorbar' to use a simple error bar for plotting.
%
% obj,plotMeanMSD(ha, errorbar, indices) computes and plots the
% mean only fothe MSD curves whose indices are given on the
% 'indices' array.
%
% h = obj,plotMeanMSD(...) returns the handle to the line plotted.
%
% [h, ha] = obj,plotMeanMSD(...) also returns the handle of the
% axes in which the curve was plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from msdanalyzer https://tinevez.github.io/msdanalyzer/
% Copyright (C) 2013 - 2014 Jean-Yves Tinevez
%
% -------------------------------------------------------------
% Copyright (C) 2019 J.J. Erik Maris
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
if errorbarFlag == 1
    hps  = gobjects(nObj,4);
elseif strcmp(errorbarFlag,'errorbar')
    hps  = gobjects(nObj,1);
else
    hps  = gobjects(nObj,1);
end

for ii = 1:nObj
    
    if ~obj(ii).MSD_valid
        obj(ii) = obj(ii).computeMSD;
    end
    
    if all_indices
        indices = 1 : obj(ii).nTracks;
    end

    [xUnitFactor,xUnitLabel] = obj(ii).getUnitFactor('dt');
    [yUnitFactor,yUnitLabel] = obj(ii).getUnitFactor('pixelsize.^2');
    
    colors = obj(ii).getColormap(nObj);

    msmsd = obj(ii).getMeanMSD(indices);
    
    if isempty( msmsd ) % do not plot empty spots in MSD
        continue
    end

    popName = sprintf('Population %d', ii );

    t = msmsd(:,1).*xUnitFactor;
    m = msmsd(:,2).*yUnitFactor;
    if errorbarFlag == 1
        s = msmsd(:,3).*yUnitFactor; % standard error of the mean
        hp = tracks.errorShade(ha, t, m, s, colors(ii,:), true);
        hps(ii,:) = [hp.mainLine hp.edge(1) hp.edge(2) hp.patch];
        set( hps(ii,:), 'DisplayName', popName );
    elseif strcmp(errorbarFlag,'errorbar') %
        s = msmsd(:,3).*yUnitFactor; % standard error of the mean
        hps(ii) = errorbar(ha, t, m, s, ...
            'Color', colors(ii,:), ...
            'LineWidth', obj(ii).lineWidth, ...
            'DisplayName', popName );
        set( hps(ii), 'DisplayName', popName );
    else
        hps(ii) = plot(ha, t, m,  ...
            'Color', colors(ii,:), ...
            'LineWidth', obj(ii).lineWidth, ...
            'DisplayName', popName );
        set( hps(ii), 'DisplayName', popName );
    end
    hold(ha,'on') 
end

xlabel(['Delay time (' xUnitLabel ')'])
ylabel(['MSD (' yUnitLabel ')'])

if nObj > 1
    legend(hps(:,1))
end


if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end
end