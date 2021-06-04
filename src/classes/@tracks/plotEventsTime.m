function varargout = plotEventsTime(obj, ha, eventName, indices, binsize, timeRange)

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
if nargin < 4 || isempty(indices)
    all_indices = true;
end

if nargin < 5 || isempty(binsize)
    binsize = 1;
end

all_timeRange = false;
if nargin < 6 || isempty(timeRange)
    all_timeRange = true;
end

nObj = numel(obj);
hp  = gobjects(nObj,1);

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
        
        [x,N,yl] = getEventTime(obj(ii),indices,eventName,minval,maxval,xunit,binsize);
        
        hp(ii) = plot(ha,x,N, ...
        	'Color', colors(ii,:), ...
            'LineWidth', obj(ii).lineWidth, ...
            'DisplayName', ['Population ' num2str(ii)]);
        
        hold(ha, 'on');
        
    end
end

hold off

xlabel(['Time (' uL ')'])
ylabel(yl)

if nObj > 1
    legend
end

if nargout > 0
    varargout{1} = hp;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

function [X,out,l] = getEventTime(obj,indices,eventName,minval,maxval,xunit,binsize)
    
    [uF,uL] = obj.getUnitFactor(xunit);
    binsize = binsize*uF;
    Xedges = minval - 0.5*uF:binsize:maxval + 0.5*uF;
    X = Xedges(1:end-1) + 0.5.*diff(Xedges);
    
    switch eventName
        case 'Reaction rate'
            data = obj.firstFrame(indices);
            data = data(data > 1); % do not take first
            out = histcounts(data.*uF,Xedges);
            out = out./binsize;
            l  = ['Reaction rate (counts/' uL ')'];
        case 'Blinking rate ON->OFF'
            out = histcounts(vertcat(obj.ONOFF{indices}).*uF,Xedges);
            out = out./binsize;
            l  = ['Blinking rate ON->OFF (counts/' uL ')'];
        case 'Blinking rate OFF->ON'
            out = histcounts(vertcat(obj.OFFON{indices}).*uF,Xedges);
            out = out./binsize;
            l  = ['Blinking rate OFF->ON (counts/' uL ')'];
        case 'Number of points ON'
            out = histcounts(vertcat(obj.time{indices}).*uF,Xedges);
        case 'Number of points ON+OFF'
            temp = cellfun(@(x) x(1):x(end), obj.time(indices),'UniformOutput', false);
            out = histcounts(vertcat(temp{:}).*uF,Xedges);
        otherwise
            error('%s is not recognized',eventName)
    end
    
end