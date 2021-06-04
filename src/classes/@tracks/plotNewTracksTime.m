function varargout = plotNewTracksTime(obj, ha, indices, timeRange)
% plotNewTracksTime Plot number of new tracks over time.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj.plotNewTracksTime plots number of new tracks over time in the current
% axes.
%
% obj.plotNewTracksTime(ha) plots the number of new tracks over time in the
% axes specified by the handle ha.
%
% obj.plotNewTracksTime(ha, indices) plots the number of new tracks over
% time for the particles with the specified indices only. Leave empty to
% plot for all particles.
%
% obj.plotNewTracksTime(ha, indices, timeRange) where timeRange is
% to-be-plotted time domain in the current units.
%
% hps =  obj.plotNewTracksTime(...) returns the handle array for the
% lines generated.
%
% [hps, ha] =  obj.plotNewTracksTime(...) also return the axes handle in
% which the lines were plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
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

xunit = 's';

if nargin < 2 || isempty(ha)
    ha = gca;
end

all_indices = false;
if nargin < 3 || isempty(indices)
    all_indices = true;
end

all_timeRange = false;
if nargin < 4 || isempty(timeRange)
    all_timeRange = true;
end

nObj = numel(obj);
hp  = gobjects(nObj,1);
pd = [];

for ii = 1:nObj
    if obj(ii).nTracks > 0
        
        if ~obj(ii).blinking_valid
            obj(ii) = obj(ii).computeBlinking;
        end
        
        if all_indices
            indices = 1 : obj(ii).nTracks;
        end
        
        nTracks = numel(indices);
        
        [uF,uL] = obj(ii).getUnitFactor(xunit);
        
       	if nObj > 1
            colors = obj(ii).getColormap(nObj);
            colors = repmat(colors(ii,:),[nTracks 1]);
        else
            colors = obj(ii).getColormap(nTracks);
        end
        
        if all_timeRange
            minval = min(obj(ii).firstFrame.*uF);
            maxval = max(obj(ii).firstFrame.*uF);
        else
            minval = timeRange(1);
            maxval = timeRange(2);
        end
        
        N = histcounts(obj(ii).firstFrame(indices).*uF,minval - 0.5*uF:uF:maxval + 0.5*uF);
        
        hp(ii) = plot(ha,minval:uF:maxval,N./uF, ...
        	'Color', colors(ii,:), ...
            'LineWidth', obj(ii).lineWidth, ...
            'DisplayName', ['Population ' num2str(ii)]);
        
        hold(ha, 'on');
        
    end
end

hold off

xlabel(['Time (' uL ')'])
ylabel(['Reactivity (counts/' uL ')'])

if nObj > 1
    legend on
end

if nargout > 0
    varargout{1} = hp;
    if nargout > 1
        varargout{2} = ha;
    end
end