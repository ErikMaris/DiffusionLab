function varargout = plotPCAScorePlot(obj, ha, properties)
% plotPCAScorePlot Plots a score plot of the computed properties
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% plotPCAScorePlot(obj) plots the PCA score plot for all computed
% properties.
%
% plotPCAScorePlot(obj, ha) plots in the axis handle ha.
%
% plotPCAScorePlot(obj, ha, properties) only uses the properties
% compute and plot the score plot. Properties is a cell with the property
% name (VariableDescriptions) or property variable name (VariableNames).
% Only Properties with a single number are taken. If empty, all propertues
% are used.
%
% h = plotPCAScorePlot(...) returns the handle to the line plotted.
%
% [h, ha] = plotPCAScorePlot(...) also returns the handle of the axes in which 
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

if nargin < 3
    properties = [];
end

features = obj.exportPropertyTable(true,true,false);

if isempty(features)
    error('No properties computed')
end

% clean all variables with two columns
clearProperties = zeros(size(features.Properties.VariableNames)); % all properties to zero
for ii = 1:size(features,2)
    if numel(features{1,ii}) > 1
        clearProperties(ii) = 1;
    end
end
features(:,logical(clearProperties)) = [];



% select properties given by user
if ~isempty(properties)
    
    % undocumented GUI to select features
    if strcmp(properties,'GUI')
        features(:,ismember(features.Properties.VariableNames,'Classifier')) = []; % clear the classifier as property
        idx = listdlg('ListString',features.Properties.VariableDescriptions,'SelectionMode','multiple','PromptString','Please select the features used for segmentation:','InitialValue',1:numel(features.Properties.VariableDescriptions));
        if isempty(idx)
            return % user pressed cancel
        end
        properties = features.Properties.VariableDescriptions(idx);
    end
    
    % select properties
    if ischar(properties) || isstring(properties)
        properties{1} = properties; % put in cell
    end
    
    takeProperties = zeros(size(features.Properties.VariableNames)); % all properties to zero
    
    for ii = 1:numel(properties) % loop over cell (array)
        takeProperties(ismember(features.Properties.VariableNames,properties)) = 1;
        takeProperties(ismember(features.Properties.VariableDescriptions,properties)) = 1;
    end
    
    features = features(:,logical(takeProperties)); % take columns in properties (takeProperties)
else
    features(:,ismember(features.Properties.VariableNames,'Classifier')) = []; % clear the classifier as property
end

X = table2array(features);
X = zscore(X);
Xnames = features.Properties.VariableDescriptions(:);
[coeff,score,~,~,explained] = pca(X);

% make a table to display the variance explained in the command window
T = array2table([(1:numel(explained))' cumsum(explained)],'VariableNames',{'Principal component','Cummulative variance explained'}); 
disp(T)

nFeatures = size(X,2);
if nFeatures < 3
    error('Please select at least 3 features')
else
    axes(ha); %set the current axes to ha
    h = biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',Xnames);
end

if nargout > 0
    varargout{1} = h;
    if nargout > 1
        varargout{2} = ha;
    end
end


end