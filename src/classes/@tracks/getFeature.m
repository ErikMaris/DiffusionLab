function [feature,varargout] = getFeature(obj,featureName,noCellFlag)
% Deals with referencing fields within a structure within an
% object. Output is cell(nObj,1) if noCellFlag = false,
% otherwise if noCellFlag = true only the first object is taken
% and the output is given as is, i.e. without cell.

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

if nargin < 3
    noCellFlag = false;
end

feature = cell(numel(obj),1);
featureUnit = cell(numel(obj),1);
featureDescription = cell(numel(obj),1);
for ii = 1:numel(obj) % I don't see a way to avoid a loop here :-(
    try
        feature{ii} = obj(ii).features.(featureName);
        idx = ismember(obj(ii).features.Properties.VariableNames,featureName);
        % don't retrieve others
        try
            featureUnit{ii} = obj(ii).features.Properties.VariableUnits{idx};
        catch
            featureUnit{ii} = [];
        end
        try
            featureDescription{ii} = obj(ii).features.Properties.VariableDescriptions{idx};
        catch
            featureDescription{ii} = [];
        end
    catch
        feature{ii} = [];
        featureUnit{ii} = [];
        featureDescription = [];
    end
end

if noCellFlag % could be faster MEM-wise if cat of all fieldvals is prevented
    feature = feature{1};

end

if nargout > 1
    if numel(featureUnit) > 1 && numel(unique(featureUnit)) > 1
        warning('Different unit definitions, picking first.')
    end
    varargout{1} = featureUnit{1};
    if nargout > 2
        if isempty(featureDescription)       
            varargout{2} = [];
        else
            if numel(featureDescription) > 1 && numel(unique(featureDescription)) > 1
                warning('Different descriptions, picking first.')
            end
            varargout{2} = featureDescription{1};
        end
    end
end
end

