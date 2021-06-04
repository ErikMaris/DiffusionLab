function varargout = plotBoxProp(obj, ha, property, indices)
% plotBoxProp Plots the boxplot of a property per population
% 
% The central mark indicates the mean, the box top and bottom edges the 
% 25th and 75th percentiles, respectively, and the whiskers the most 
% extreme data points not considered outliers. Outliers are indicated with 
% plus symbols. 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotBoxProp(obj, ha, property, indices, units, propertyLabel)
%   obj: allows >1 objs to be parsed
%   ha: handle to axes. If empty the current axis is taken
%   property / string: name of the property to be plotted
%   indices: indices of the tracks to be displayed
%   units: display units of the property
%   propertyLabel / string: display label of property
%
% h = plotBoxProp(...) returns the handle to the line plotted.
%
% [h, ha] = plotBoxProp(...) also returns the handle of the axes in which 
% the curve was plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
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

nObj = numel(obj);
data = cell(nObj,1);
popName = cell(nObj,1);

for ii = 1:nObj
    [data,units,propertyLabel] = obj(ii).getFeature(property);
    data{ii} = unpackCells(obj,units,data,indices); % // deals with empty indices
end

Ndata = cellfun(@(x) numel(x), data);

if sum(Ndata) == 0
    return
end

gr = cell(nObj,1);
for ii = 1:nObj
    gr{ii} = repmat({sprintf('Population %d', ii )},Ndata(ii),1);
end
data = vertcat(data{:});
gr = vertcat(gr{:});
boxplot(data,gr)

[~,uL] = obj(1).getUnitFactor(units); % empty units accepted

if isempty(uL)
    ylabel(propertyLabel)
else
    ylabel([propertyLabel ' (' uL ')'])
end

if nargout > 0
    varargout{1} = ha.Children;
    if nargout > 1
        varargout{2} = ha;
    end
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

