function varargout = plotdR2(obj, ha, indices, delayIdx, logbinFlag, normalization, Xedges)
% plotdR2 Plots a histogram of the squared displacement per population
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotSqDisp(obj, ha, indices, delayIdx, logbinFlag, normalization)
%   obj: allows >1 objs to be parsed
%   ha: handle to axes. If empty the current axis is taken
%   indices: indeces of the tracks to be displayed
%   delayIdx: the index/indices of the delay time(s) to be plotted
%   units: display units of the property. Default is empty
%   logbinFlag: bin logarithmically instead of linear
%   normalization / string: default 'probaility', other options are 
%   'count', 'countdensity', 'pdf', 'cumcount', or 'cdf'.
%
% h = plotSqDisp(...) returns the handle to the line plotted.
%
% [h, ha] = plotSqDisp(...) also returns the handle of the axes in which 
% the curve was plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% 
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

if nargin < 4 || isempty(delayIdx)
    delayIdx = 1;
end

if nargin < 5
    logbinFlag = false;
end

if nargin < 6
    normalization = 'probability';
end

if nargin < 7 || isempty(Xedges)
    Xedges = false;
end

nObj = numel(obj);

if nObj > 1 && numel(delayIdx) > 1
    error('tracks:plotSqDisp:multipleObjDelayIdx',...
        'Cannot plot multiple delay indices for multiple objects.')
end

all_indices = false;
if nargin < 3 || isempty(indices)
    all_indices = true;
end

if nObj > 1 % plot multiple populations; get cell with this
    data = cell(nObj,1);
    for ii = 1:nObj
        if all_indices
            indices = 1:obj(ii).nTracks;
        end
        data{ii} = obj(ii).getdR2(delayIdx,indices);
        data{ii} = data{ii}{1}; % only one delayIdx
        [uF,uL] = obj(ii).getUnitFactor('pixelsize.^2');
        data{ii} = data{ii}.*uF;
    end
    nPlot = numel(obj);
else % plot multiple delay times; get cell with this
    if all_indices
        indices = 1:obj.nTracks;
    end
    [data, dt] = obj.getdR2(delayIdx,indices);
    [uF,uL] = obj.getUnitFactor('pixelsize.^2');
    data = cellfun(@(x) x.*uF, data,'UniformOutput',false);
    nPlot = numel(delayIdx);
end

if ~Xedges % if edges are not specified by user
    Xedges = getXedges(vertcat(data{:}),logbinFlag); % get edges based on all data points
end

for ii = 1:nPlot
    
    if isempty(data{ii})
        warning('tracks:plotSqDisp:EmptySquaredDisplacement',...
        'Squared displacement of object %i has not been calculated. Skipping object.',ii)
        continue
    end

    if nObj > 1
        jj = ii;
        plotName = sprintf('Population %d', ii );
    else
        jj = 1;
        [uFtime, uLtime] = obj(jj).getUnitFactor('dt');
        plotName = sprintf('Delay time %f%s', dt(ii).*uFtime, uLtime);
    end
    
    colors = obj(jj).getColormap(nPlot);
    LW = obj(jj).lineWidth;
    
    histogram(ha, ...
        data{ii}, Xedges, ...
        'FaceColor', colors(ii,:), ...
        'LineWidth', LW, ...
        'DisplayName', plotName, ...
        'Normalization', normalization);
    hold(ha, 'on');
end

hold(ha, 'off');

if logbinFlag
    set(ha,'XScale','log');
end

if nPlot > 1
    legend
end

xlabel(['Squared displacement (' uL ')'])
ylabel(normalization)


if nargout > 0
    varargout{1} = ha.Children;
    if nargout > 1
        varargout{2} = ha;
    end
end
end

function Xedges = getXedges(data,logbinFlag)
    if logbinFlag
        [~,Xedges] = histcounts(log10(data)); % in log
        Xedges = [Xedges  (1:5) .* diff(Xedges(1:2)) + Xedges(end)]; % add some trailing X values
        Xedges = 10.^Xedges; % in lin scale
    else
        [~,Xedges] = histcounts(data);
    end
end

