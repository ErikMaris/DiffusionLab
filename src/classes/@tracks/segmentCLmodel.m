function tseg = segmentCLmodel(obj,model,SIflag,excludeIncomplete,addUnitsToHeader)
% segmentCLmodel Segments tracks based on a Classification Trainer model.
%
% tseg = segmentCLmodel(obj,model) segments the tracks in obj based on the 
% Classification Trainer model.
%
% tseg = segmentCLmodel(obj,model,SIflag) where SIflag resets the units to
% SI (true) or keeps the current units (false).
%
% tseg = segmentCLmodel(obj,model,SIflag,excludeIncomplete) where
% excludeIncomplete exludes the features that are incomplete for one or 
% more track objects and uses these for the segmentation (true) or uses the
% incomplete set.
%
% tseg = segmentCLmodel(obj,model,SIflag,excludeIncomplete,addUnitsToHeader)
% where addUnitsToHeader markers wheter the unit names are included in the
% header or feature name of the model (true/false). This depends on the
% feature names used in the Classification Trainer.
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

if numel(obj) > 1
    error('Can only segment a single object.')
end

% use same defaults as 'exportPropertyTable' to keep compatibility
if nargin < 3
    SIflag = false;
end

if nargin < 4
    excludeIncomplete = false;
end

if nargin < 5
    addUnitsToHeader = true;
end

T = obj.exportPropertyTable(SIflag,excludeIncomplete,addUnitsToHeader);

% Classification Learner generates valid variable names. We also have
% to do this with our headers
T.Properties.VariableNames = strrep(T.Properties.VariableNames,' ','');
T.Properties.VariableNames = strrep(T.Properties.VariableNames,'^','');
T.Properties.VariableNames = strrep(T.Properties.VariableNames,'/','');

% check whether they are all available
if any(~ismember(lower(model.RequiredVariables),lower(T.Properties.VariableNames))) % case insensitive
    error('Missing variables: %s.',...
        strjoin(model.RequiredVariables(~ismember(lower(model.RequiredVariables),lower(T.Properties.VariableNames))),', '))
end


% Classification Learner seems to be quite incosistent with capitalization.
% We make the headers the same brute force-style. Wish 
tableVarNames = T.Properties.VariableNames;
for ii = 1:numel(model.RequiredVariables)
    b = strcmpi(tableVarNames,model.RequiredVariables{ii});
    tableVarNames{b} = model.RequiredVariables{ii};
end

T.Properties.VariableNames = tableVarNames;

Y = model.predictFcn(T);

u = unique(Y);

nSeg = numel(u);

if numel(u) == 1
    tseg = obj; % all of same category
    tseg.name = u;
else
    tseg = tracks.empty(nSeg,0);
    for ii = 1:nSeg
        tseg(ii) = obj.segment(ismember(Y,u(ii))); % works for both numeric and string (= categorical)
        tseg(ii).name = u(ii);
    end
end

end

