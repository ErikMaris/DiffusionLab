function features = exportPropertyTable(obj,SIflag,excludeIncomplete,addUnitsToHeader)
% exportPropertyTable Exports the features to a table.
%
%
% features = exportPropertyTable(obj,SIflag) where SIflag resets the units 
% to SI (true) or keeps the current units (false). Features is a table
% containing all features.
%
% features = exportPropertyTable(obj,model,SIflag,excludeIncomplete) where
% excludeIncomplete exludes the features that are incomplete for one or 
% more track objects and uses these for the segmentation (true) or uses the
% incomplete set.
%
% features = exportPropertyTable(obj,model,SIflag,excludeIncomplete,
% addUnitsToHeader) where addUnitsToHeader markers wheter the unit names 
% are included in the header or feature name of the model (true/false). 
% This depends on the feature names used in the Classification Trainer.
%

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

if nargin < 2
    SIflag = false;
end

if nargin < 3
    excludeIncomplete = false;
end

if nargin < 4
    addUnitsToHeader = true;
end

if SIflag
    % reset usit system to SI
    obj = obj.resetUnitSystem;
end

% check if all are in the same unit system
if numel(obj) > 1 
    if ~isequal(obj.unitSystem)
        error('Unit systems differ between objects. Cancel export. Consider setting the SIflag to ''true'' for export.')
    end
end

nObj = numel(obj);

% collect all property names
propertyNames = cell(nObj,1);
features = cell(nObj,1);
emptymarker = zeros(nObj,1);
for ii = 1:nObj
    features{ii} = obj(ii).features;
    if isempty(features{ii})
        emptymarker(ii) = 1;
        continue
    end
    propertyNames{ii} = obj(ii).features.Properties.VariableNames;
end

if sum(emptymarker) == nObj % if no features are available, return
    features = [];
    return
end

propertyNames = unique(horzcat(propertyNames{:}));

% give tables same variables
if ~excludeIncomplete
    for ii = 1:nObj
        add = ~ismember(propertyNames,features{ii}.Properties.VariableNames);
        addNames = {propertyNames{add}};
        for jj = 1:numel(addNames)
            features{ii}.(addNames{jj}) = nan(size(features{ii},1),1);
        end
    end
else
    keep = ones(1,numel(propertyNames));
    for ii = 1:nObj
        keep = keep.*ismember(propertyNames,features{ii}.Properties.VariableNames);
    end
    keepNames = propertyNames(~keep);
    for ii = 1:nObj
        clear = ismember(features{ii}.Properties.VariableNames,keepNames);
        features{ii}(:,clear) = [];
    end
end

tableSize = nan(1,nObj);
objNames = cell(nObj,1);

% convert units
for ii = 1:nObj
    for jj = 1:size(features{ii},2)
        if isempty(features{ii}.Properties.VariableUnits{jj})
            continue
        end
        [uF,uL] = obj(ii).getUnitFactor(features{ii}.Properties.VariableUnits{jj});
        features{ii}.(features{ii}.Properties.VariableNames{jj}) = ...
            features{ii}.(features{ii}.Properties.VariableNames{jj}) .* uF;
        features{ii}.Properties.VariableUnits{jj} = uL;
    end
    % store info for classifiers
    tableSize(ii) = size(features{ii},1);
    objNames{ii} = obj(ii).name;
end

% merge table
features = vertcat(features{:});

% clear trackID
try
    features.trackID = [];
catch
    % do nothing
end

% generate classifiers
% if one object does not have a name, just number them
if ~any(~cellfun(@isempty,objNames)) 
    classifier = repelem(1:nObj,tableSize)';
else % if they all have a name, add as cell array
    classifier = cell(nObj,1);
    for jj = 1:nObj
        classifier{jj} = cell(tableSize(jj),1);
        % fill in name with 'No name #'
        if isempty(objNames{jj})
            objNames{jj} = ['No name ' num2str(jj)];
        end
        classifier{jj}(:) = objNames(jj);
    end
    classifier = string(vertcat(classifier{:}));
end
% add classifier
features = addvars(features,classifier,'Before',features.Properties.VariableNames{1},'NewVariableNames','Classifier'); % add before first property to make sure it taken as classifier by the Classification Learner
features.Properties.DimensionNames = {'Tracks','Track properties'};

if addUnitsToHeader
    % add units to header
    for ii = 1:size(features,2)
        if isempty(features.Properties.VariableUnits{ii})
            continue
        end
        features.Properties.VariableNames{ii} = [features.Properties.VariableNames{ii} '_' features.Properties.VariableUnits{ii}];
    end
end

end

% 
% 
% % proallocate cell with (properties,object)
% propCell = cell(size(TrackProperties,1)+1,nObj); % one extra for classification
% headers = cell(size(TrackProperties,1),1);
% 
% % get all properties; ouch so ugly..
% for jj = 1:size(TrackProperties,1)
%     for jj = 1:nObj
%         propCell(jj,jj) = obj(jj).getFeature(TrackProperties{jj,2},false);
%         % check if column does exist
%         if size(propCell{jj,jj},2) < TrackProperties{jj,4}
%             % if not, make empty
%             propCell{jj,jj} = [];
%             continue % for unit label
%         else
%             % if does exist, get column and do unit conversion
%             propCell{jj,jj} = propCell{jj,jj}(:,TrackProperties{jj,4}).*obj(jj).getUnitFactor(TrackProperties{jj,3});
%         end
%         % construct headers
%         [~,unitLabel] =  obj(1).getUnitFactor(TrackProperties{jj,3});
%         if ischar(unitLabel)
%             headers{jj} = [TrackProperties{jj,1} '_' unitLabel];
%         else
%             headers{jj} = TrackProperties{jj,1};
%         end
%     end
% end
% 
% % take only ones which are available for all populations
% % check for empty
% empty = cellfun(@isempty,propCell);
% empty = sum(empty,2);
% % only take which are non-empty for all track objects
% take = (empty == 0);
% propCell = propCell(take,:); % do not overwrite last line for classification 
% headers = headers(take,:);
% 
% % check equal length of all
% numelProp = cellfun(@numel,propCell);
% objNames = cell(nObj,1);
% for jj = 1:nObj
%     if ~unique(numelProp(:,jj))
%         error('Number of property elements are not identical for object %i.',jj)
%     end
%     objNames{jj} = obj(jj).name;
% end
% numelProp = numelProp(1,:);
% 
% % construct table
% T = array2table(cell2mat(propCell'),'VariableNames',headers);