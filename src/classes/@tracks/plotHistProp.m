function varargout = plotHistProp(obj, ha, property, indices, logbinFlag, normalization)
% plotHistProp Plots a histogram of a property per population
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotHistProp(obj, ha, property, indices, units, logbinFlag, normalization, propertyLabel)
%   obj: allows >1 objs to be parsed
%   ha: handle to axes. If empty the current axis is taken
%   property / string: name of the property to be plotted
%   indices: indeces of the tracks to be displayed
%   units: display units of the property. Default is empty
%   logbinFlag: bin logarithmically instead of linear
%   normalization / string: default 'probaility', other options are 
%   'count', 'countdensity', 'pdf', 'cumcount', or 'cdf'.
%   propertyLabel / string: display label of property
%
% h = plotHistProp(...) returns the handle to the line plotted.
%
% [h, ha] = plotHistProp(...) also returns the handle of the axes in which 
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

% try obj.getProp(property);
% catch
%     error('tracks:plotHistProp:PropertyDoesNotExist',...
%         'Property %s does not exist in all objects.',property)
% end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 4
    indices = [];
end

if nargin < 5
    logbinFlag = false;
end

if nargin < 6
    normalization = 'probability';
end

nObj = numel(obj);

[data,units,propertyLabel] = obj.getFeature(property);

Xedges = getXedges(unpackCells(obj,units,data,indices),logbinFlag); % get edges based on all data points // deals with empty indices

for ii = 1:nObj

    
    [~,uL] = obj(ii).getUnitFactor(units); % empty units accepted
    
    colors = obj(ii).getColormap(nObj);

    popName = sprintf('Population %d', ii );
    
    data = unpackCells(obj,units,obj(ii).getFeature(property),indices); % // deals with empty indices
    
    if isempty(data)
        continue
    end
    
    histogram(ha, ...
        data, Xedges, ...
        'FaceColor', colors(ii,:), ...
        'LineWidth', obj(ii).lineWidth, ...
        'DisplayName', popName, ...
        'Normalization', normalization);
    hold(ha, 'on');
end

hold(ha, 'off');

if logbinFlag
    set(ha,'XScale','log');
end

if isempty(uL)
    xlabel(propertyLabel)
else
    xlabel([propertyLabel ' (' uL ')'])
end
ylabel(normalization)


% legendTxt = cell(nObj,1);
% for ii = 1:nObj
%     legendTxt{ii} = ['Population ' num2str(ii) ];
% end
% legend(legendTxt{:})
if nObj > 1
    legend
end


if nargout > 0
    varargout{1} = ha.Children;
    if nargout > 1
        varargout{2} = ha;
    end
end
end

function Xedges = getXedges(data,logbinFlag)
    if logbinFlag
        dataLTzero = data < 0;
        s = sum(dataLTzero);
        if sum(dataLTzero) > 0
            data = data(~dataLTzero);
            warning('tracks:plotHistProp:negValsLogBinning',...
                'Data contains negative values, which cannot be binned logarithmically. Throwing %i data points.',s)
        end
        [~,Xedges] = histcounts(log10(data)); % in log
        Xedges = [Xedges  (1:5) .* diff(Xedges(1:2)) + Xedges(end)]; % add some trailing X values
        Xedges = 10.^Xedges; % in lin scale
    else
        [~,Xedges] = histcounts(data);
    end
end

function data = unpackCells(obj,unit,data,indices)
for ii = 1:numel(data) % unpack first layer and convert with unit factor
    while iscell(data{ii}) % unpack fully
        data{ii} = vertcat(data{ii}{:});
    end
    if ~isnumeric(data{ii})
        data{ii} = []; % sometimes the user refers to an object with a property which is empty. This gives non-numeric output and cannot be plotted
        continue
    end
    if ~isempty(indices)
        data{ii} = data{ii}(indices,:); % do data selection
    end
    if ~isempty(unit)
        data{ii} = data{ii}.*obj(ii).getUnitFactor(unit); % convert to desired units
    end
end
data = vertcat(data{:});
end

