function varargout = plotVoronoiDiag(obj,ha,property,noEdgeFlag,noCoMFlag,colourBool)
% plotVoronoiDiag Plots a Voronoi diagram of a property per population
% 
% Voronoi diagram of the tracks center of mass coordinates. The tracks at 
% the edges are enclosed by a convex hull. 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj(ii).plotVoronoiDiag plots the particle Voronoi diagram stored in the
% tracks object obj in the current axes. It shows the center of mass and 
% the edges. This method only works for 2D problems. If multiple objects
% are parsed, the 
%
% obj(ii).plotVoronoiDiag(ha) plots the diagram in the axes with
% handle ha.
%
% obj(ii).plotVoronoiDiag(ha, property) where the property is a string
% with the name of the to-be-plotted property.
%
% obj(ii).plotVoronoiDiag(ha, property, noEdgeFlag) where the noEdgeFlag is
% a boolean, allows to specify whether the edges of the Voronoi diagram
% should be plotted.
%
% obj(ii).plotVoronoiDiag(ha, property, noEdgeFlag, noCoMFlag) where the 
% noCoMFlag is a boolean, allows to specify whether the center of mass of 
% the tracks, which were used to construct the Voronoi diagram should be 
% plotted.
%
% obj(ii).plotVoronoiDiag(ha, property, noEdgeFlag, noCoMFlag, colourBool)
% where the colourBool is a logical vector indicating whether the Voronoi
% cell of track should be coloured. The track indices continue counting
% through mutilple objects.
%
% h = plotVoronoiDiag(...) returns the handle to the line plotted.
%
% [h, ha] = plotVoronoiDiag(...) also returns the handle of the axes in which 
% the curve was plotted.

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from Florian Meirer
% (C) 2016
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

nObj = numel(obj);

if nObj == 1
    try obj.getFeature(property);
    catch
        error('tracks:plotVoronoiDiag:PropertyDoesNotExist',...
            'Property %s does not exist in all objects.',property)
    end
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 4
    noEdgeFlag = false;
end

if nargin < 5
    noCoMFlag = false;
end

nTracks = [obj.nTracks];

if nargin < 6
    colourBool = ones(sum(nTracks),1); % take all
end

xyBorders = [10 10];

hp = [];
% hp = NaN(sum(nTracks),1);

if nObj > 1
    colors = cell(nObj,1);
    for ii = 1:nObj
        tcolors = obj(ii).getColormap(nObj);
        colors{ii} = repmat(tcolors(ii,:),[nTracks(ii) 1]);
    end
    colors = vertcat(colors{:});
else
    colors = cell(nObj,1);
    for ii = 1:nObj
        [colors{ii},units] = obj(ii).getFeature(property,true);
        colors{ii} = colors{ii} .* obj(ii).getUnitFactor(units);
    end
    colors = vertcat(colors{:});
%     colors = colors(:,1); % only take first column
end

CoM = cell(nObj,1);
for ii = 1:nObj
    if ~ismember('CenterOfMass',obj(ii).features.Properties.VariableNames)
        obj(ii) = obj(ii).computeFeatures(false,'CenterOfMass');
    end
    CoM{ii} = obj(ii).features.CenterOfMass;
end

CoM = vertcat(CoM{:});

% get the convex hull and extend it by 5 pixels; then add the hull points
% to the data to cut of the voronoi diagram at the extended hull
hullPs = obj(1).getHullPoints(CoM,10);

startOfhPs = size(CoM,1)+1;
MSDxyOrig = CoM;
CoM = vertcat(CoM,hullPs);
[v,c] =  voronoin(CoM);

for ii = 1:length(c)
    if all(c{ii}~=1)  % If at least one of the indices is 1,
                      % then it is an open region and we can't
                      % patch that.
        if ii < startOfhPs % only patch real data points
            if noEdgeFlag
                p = patch(v(c{ii},1),v(c{ii},2),colors(ii,:),'EdgeColor','none');
            else
                p = patch(v(c{ii},1),v(c{ii},2),colors(ii,:));
            end
            if colourBool(ii) == 0
                p.FaceAlpha = 0;
            end
            hp = [hp p];
        end
    end
end  

xBorder = xyBorders(1);
yBorder = xyBorders(2);
xlim([min(CoM(:,1))-xBorder max(CoM(:,1))+xBorder]);
ylim([min(CoM(:,2))-yBorder max(CoM(:,2))+yBorder]);
CoM = MSDxyOrig;

if ~noCoMFlag
    hold on;
    plot(CoM(:,1),CoM(:,2),'k.','LineWidth',obj(1).lineWidth);
end
    
axis image;

if nObj > 1
    % to get the legend for the populations right, we use the last plotted
    % object in the axes to make our legend.
    legendTxt = cell(nObj,1);
    legendHps = zeros(nObj,1);
    idx = cumsum(nTracks);
    for ii = 1:nObj
        legendTxt{ii} = ['Population ' num2str(ii) ];
        legendHps(ii) = hp(idx(ii));
    end
    legend(legendHps,legendTxt{:})
else
    colorbar
end


if nargout > 0
    varargout{1} = hp;
    if nargout > 1
        varargout{2} = ha;
    end
end


end