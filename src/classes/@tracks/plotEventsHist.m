function varargout = plotEventsHist(obj, ha, eventName, indices, timeRange, suppressFit)

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

all_timeRange = false;
if nargin < 5 || isempty(timeRange)
    all_timeRange = true;
end

if nargin < 6
    suppressFit = false;
end

nObj = numel(obj);
hp  = gobjects(nObj,2);
lambda = [];

for ii = 1:nObj
    if obj(ii).nTracks > 0
        
        if ~obj(ii).blinking_valid
            obj(ii) = obj(ii).computeBlinking;
        end
        
        if all_indices
            indices = 1 : obj(ii).nTracks;
        end
        
        nTracks = numel(indices);
        
        uF = obj(ii).getUnitFactor(xunit);
        
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
        
        [Nt,xl,lambda] = getEventTime(obj(ii),indices,eventName,minval,maxval,xunit);
        
        if ~isempty(lambda) && suppressFit == false
            hp(ii,:) = histfitpoisson(Nt(2:end),lambda);
%             hp(ii,:) = histfit(Nt(2:end),[],'Poisson');
            hp(ii,1).FaceColor =  colors(ii,:);
            hp(ii,1).LineWidth = obj(ii).lineWidth;
            hp(ii,1).DisplayName = ['Population ' num2str(ii)];
            hp(ii,2).Color =  colors(ii,:);
            hp(ii,2).LineWidth = obj(ii).lineWidth;
            hp(ii,2).LineStyle = '--';
            hp(ii,2).DisplayName = ['Lambda ' num2str(lambda)];
        else
            hp(ii,1) = histogram(Nt(2:end));
            hp(ii,1).FaceColor =  colors(ii,:);
            hp(ii,1).LineWidth = obj(ii).lineWidth;
            hp(ii,1).DisplayName = ['Population ' num2str(ii)];
        end
        
        hold(gca, 'on');
        
    end
end

hold off

xlabel(xl)
ylabel('Counts')

if nObj > 1 || suppressFit == false
    legend % to show lambda
end

if nargout > 0
    varargout{1} = hp(:);
    if nargout > 1
        varargout{2} = ha;
        if nargout > 1
            varargout{3} = lambda;
        end
    end
end

end

function [out,l,lambda] = getEventTime(obj,indices,eventName,minval,maxval,xunit)
    
    [uF,uL] = obj.getUnitFactor(xunit);
    Xedges = minval - 0.5*uF:uF:maxval + 0.5*uF;
    lambda = [];
    
    switch eventName
        case 'Reaction rate'
            data = obj.firstFrame(indices);
            data = data(data > 1); % do not take first
            out = histcounts(data.*uF,Xedges);
            lambda = poissfit(out);
            l  = ['Reaction rate (counts/' num2str(uF) ' ' uL ')'];
        case 'Blinking rate ON->OFF'
            out = histcounts(vertcat(obj.ONOFF{indices}).*uF,Xedges);
            lambda = poissfit(out);
            l  = ['Blinking rate ON->OFF (counts/' num2str(uF) ' ' uL ')'];
        case 'Blinking rate OFF->ON'
            out = histcounts(vertcat(obj.OFFON{indices}).*uF,Xedges);
            lambda = poissfit(out);
            l  = ['Blinking rate OFF->ON (counts/' num2str(uF) ' ' uL ')'];
        case 'Number of points ON'
            out = histcounts(vertcat(obj.time{indices}).*uF,Xedges);
        case 'Number of points ON+OFF'
            temp = cellfun(@(x) x(1):x(end), obj.time(indices),'UniformOutput', false);
            out = histcounts(vertcat(temp{:}).*uF,Xedges);
        otherwise
            error('%s is not recognized',eventName)
    end
    
end

function hp = histfitpoisson(data,lambda)

hp = gobjects(1,2);

minval = floor(min(data))-0.5;
maxval = ceil(max(data))+0.5;
Xedges = minval:maxval;

N = numel(data);
h = 1; % binwidth 

% Histogram
hp(1) = histogram(data,Xedges);
% Overlay the distribution
hold on;
f = @(z) N.*h.*poisspdf(z,lambda);
x = minval+0.5:maxval-0.5;
hp(2) = plot(x,f(x));
hold off;
end