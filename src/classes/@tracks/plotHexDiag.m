function varargout = plotHexDiag(obj,ha,binSize,property,binType,logFlag,noEdgeFlag)
% plotHexDiag Plots a hex grid colour coded by the mean value of a property
%
% Function is work in progress...
%
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj.plotHexDiag plots the MSD curves in the current axes.
%
% obj.plotHexDiag(ha) plots the MSD curves in the axes
% specified by the handle ha.
%
% obj.plotHexDiag(ha,binSize) where binSize is the size of the hex in the
% current units.
%
% obj.plotHexDiag(ha,binSize,property) where property is the property used
% for the colour coding.
%
% obj.plotHexDiag(ha,binSize,property,binType) deprecated
%
% obj.plotHexDiag(ha,binSize,property,binType,logFlag) where logFlag
% indicates whether the property is plotted in a log scale (true/false).
%
% obj.plotHexDiag(ha,binSize,property,binType,logFlag,noEdgeFlag) where 
% noEdgeFlag indicates whether the edges are plotted (true/false).
%
% hps =  obj.plotHexDiag(...) returns the handle array for the
% lines generated.
%
% [hps, ha] =  obj.plotHexDiag(...) also return the axes handle in
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

% Florian Meirer
% no unit conversion

try obj.getFeature(property);
catch
    error('tracks:plotHexDiag:PropertyDoesNotExist',...
        'Property %s does not exist in all objects.',property)
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 3 || isempty(binSize)
    binSize = 1000e-9; % 100 nm
end

if nargin < 6 || isempty(logFlag)
    logFlag = false;
end

if nargin < 7 || isempty(noEdgeFlag)
    noEdgeFlag = false;
end

nObj = numel(obj);
hp = [];

% binType is CoM, make sure it is computed
for ii = 1:nObj
    if ~ismember('CenterOfMass',obj(ii).features.Properties.VariableNames)
        obj(ii) = obj(ii).computeFeatures(true,'CenterOfMass');
    end
end

% make grid

% switch binType
%     case 'localizations'
%         XY = cell(numel(obj),1);
%         for ii = 1:numel(obj)
%             XY{ii} = vertcat(obj(ii).coords{:}).*obj(ii).getUnitFactor('pixelsize');
%         end
%     case 'CoM'
%         
%     otherwise
%         error('binType: ''%s'' is not recognised.',binType)
% end

% % sometimes output is not cell, and vertcat will give an error
% if ~iscell(XY)
%     XY = {XY};
% end

XY = cell(numel(obj),1);
for ii = 1:numel(obj)
    XY{ii} = obj(ii).features.CenterOfMass.*obj(ii).getUnitFactor('pixelsize');
end


% get properties
prop = cell(numel(obj),1);
for ii = 1:numel(obj)
    [dummy,dummyUnit] = obj(ii).getFeature(property,true);
    if isempty(dummy)
        XY{ii} = [];
        continue
    end
    prop{ii} = dummy.*obj(ii).getUnitFactor(dummyUnit);
end
prop = vertcat(prop{:});

XY = vertcat(XY{:});

[v,c,idxhex] = tesselate(XY,'hexagonal',binSize);

ne_idx = unique(idxhex); % loop over unique indices
propHex = nan(numel(c),1);
for ii = 1:numel(ne_idx)
    propHex(ne_idx(ii)) = mean(prop(idxhex == ne_idx(ii)));
end

if logFlag
    if sum(propHex <= 0) > 0
        warning('Deleting hexagons with negative values')
        propHex(propHex <= 0) = []; % clear negative values
    end
    propHex = log10(propHex);
    minPropHex = min(propHex);
end

maxPropHex = max(propHex); % for setting color scale

nColour = 128;
propHex = ceil(normalize(propHex,'range').*(nColour-1)+1); % prevent 0 index

cmap = obj(1).getColormap(nColour);

% Plot the hexagonal mesh, including cell borders
% [XVhex, YVhex] = voronoi(Xhex(:),Yhex(:));


for ii = 1:length(c)
    if all(c{ii}~=1)  % If at least one of the indices is 1,
                      % then it is an open region and we can't
                      % patch that.
        if isnan(propHex(ii))
            continue
%             p.FaceAlpha = 0;
        end
        if noEdgeFlag
            p = patch(v(c{ii},1),v(c{ii},2),cmap(propHex(ii),:),'EdgeColor','none');
        else
            p = patch(v(c{ii},1),v(c{ii},2),cmap(propHex(ii),:));
        end
        hp = [hp p];
    end
end

[~,uL] = obj(1).getUnitFactor('pixelsize'); % assume same unit system in all objs
xlabel(['x (' uL ')'])
ylabel(['y (' uL ')'])

% xBorder = 0;
% yBorder = 0;
% xlim([XYmin(1)-xBorder XYmax(1)+xBorder].*);
% ylim([XYmin(2)-yBorder XYmax(2)+yBorder]);

colormap(cmap)
h = colorbar;
if logFlag
    if minPropHex == maxPropHex
        ha.CLim = [minPropHex maxPropHex + 0.01.*maxPropHex];
    else
        ha.CLim = [minPropHex maxPropHex];
    end
    ylabel(h,['log10(' property ')'])
else
    ha.CLim = [0 maxPropHex];
    ylabel(h,property)
end


axis image;

if nargout > 0
    varargout{1} = hp;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

function [v,c,idxhex] = tesselate(XY,method,binSize)
% XY in desired units

% get min and max coordinate values
XYmin = min(XY);
XYmax = max(XY);

% set bin limits
XYmin = floor(XYmin./binSize)-2; % 2 buffer bin, because border hexagons are not drawn
XYmax = ceil(XYmax./binSize)+2;

switch method
    case 'hexagonal'
        % Generate hexagonal grid with interval 1
        Rad3Over2 = sqrt(3) / 2;
        [Xhex, Yhex] = meshgrid(XYmin(1):1:XYmax(1),XYmin(2):1:XYmax(2));
        n = size(Xhex);
        Xhex = Rad3Over2 * Xhex;
        Yoff = repmat([0 0.5],[n(1),ceil(n(2)/2)]);
        % Yoff only generates even values in y-direction
        Yhex = Yhex + Yoff(:,1:size(Yhex,2));
        
        % get all the points that make up the hexagonal grid and scale size
        XYhex = [Xhex(:),Yhex(:)].*binSize;
        
        % let's assign all XY to a point XYhex that makes up the voronoi using KNN
        idxhex = knnsearch(XYhex,XY);
        
        % store in cell array sorted per cell
%         hexcont = cellfun(@(x) idxhex(idxhex == x),num2cell(1:numel(XYhex)),'UniformOutput',false);
        
        % get vertices and cells voronoi diagram
        [v,c] =  voronoin(XYhex);
        
%         T = table(v,c,hexcont,'VariableNames',{'Vertices','Cells','Contents'});
    otherwise
        error('Method ''%s'' is not recognised.',method)
end

end