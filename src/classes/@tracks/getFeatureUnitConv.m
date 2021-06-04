function feature = getFeatureUnitConv(obj,featureName,noCellFlag)
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
            warning('No units available.')
            featureUnit{ii} = [];
        end
        if ~isempty(featureUnit{ii})
            feature{ii} = feature{ii} .* obj(ii).getUnitFactor(featureUnit{ii}); % unit conversion
        end
    catch
        feature{ii} = [];
    end
end

if noCellFlag % could be faster MEM-wise if cat of all fieldvals is prevented
    feature = feature{1};

end

end

