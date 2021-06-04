function varargout = plotPropCorrelation(obj, ha, propertyX, propertyY, indices)
% plotPropCorrelation Plots the correlation of two properties
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotPropCorrelation(obj, ha, propertyX, propertyY, indices, unitX, unitY, propertyLabelX, propertyLabelY)
%   obj: allows >1 objs to be parsed
%   ha: handle to axes. If empty the current axis is taken
%   propertyX / string: name of the property to be plotted on the x-axis
%   propertyY / string: name of the property to be plotted on the y-axis
%   indices: indeces of the tracks to be displayed
%   unitsX: display symunits of propertyX. Default is empty
%   unitsY: display symunits of propertyY. Default is empty
%   propertyLabelX / string: display label of property on the x-axis
%   propertyLabelY / string: display label of property on the y-axis
%
% h = plotPropCorrelation(...) returns the handle to the line plotted.
%
% [h, ha] = plotPropCorrelation(...) also returns the handle of the axes in 
% which the curve was plotted.

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

try obj.getFeature(propertyX);
catch
    error('tracks:plotPropCorrelation:PropertyDoesNotExist',...
        'Property %s does not exist in all objects.',propertyX)
end

try obj.getFeature(propertyY);
catch
    error('tracks:plotPropCorrelation:PropertyDoesNotExist',...
        'Property %s does not exist in all objects.',propertyY)
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 5
    indices = [];
end

nObj = numel(obj);
hps  = nan(nObj,1);

for ii = 1:nObj


    
    popName = sprintf('Population %d', ii );
    colors = obj(ii).getColormap(nObj);
    
    [dataX,unitX,propertyLabelX] = obj(ii).getFeature(propertyX,true);
    [dataY,unitY,propertyLabelY] = obj(ii).getFeature(propertyY,true);
    
    [uFx,uLx] = obj(ii).getUnitFactor(unitX);
    [uFy,uLy] = obj(ii).getUnitFactor(unitY);
    
    if ~isempty(indices)
        dataX = dataX(indices,:);
        dataY = dataY(indices,:);
    end
    
    if ~isempty(uFx)
        dataX = dataX.*uFx;
    end
    if ~isempty(uFy)
        dataY = dataY.*uFy;
    end
    if isempty(dataX) || isempty(dataY)
        continue
    end

    hps(ii) = scatter(ha, dataX, dataY, 'x', ...
        'MarkerEdgeColor', colors(ii,:), ...
        'LineWidth', obj(ii).lineWidth, ...
        'DisplayName', popName );
    hold(ha, 'on');
end
hold(ha, 'off');

if isempty(uLx)
    xlabel(propertyLabelX)
else
    xlabel([propertyLabelX ' (' uLx ')'])
end
if isempty(uLy)
    ylabel(propertyLabelY)
else
    ylabel([propertyLabelY ' (' uLy ')'])
end

if nObj > 1
%     % to get the legend for the populations right, we use the last plotted
%     % object in the axes to make our legend.
%     legendTxt = cell(nObj,1);
%     legendHps = zeros(nObj,1);
%     for ii = 1:nObj
%         legendTxt{ii} = ['Population ' num2str(ii) ];
%         legendHps(ii) = hps(ii);
%     end
%     legend(legendHps,legendTxt{:})
    legend
end


if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end
end

