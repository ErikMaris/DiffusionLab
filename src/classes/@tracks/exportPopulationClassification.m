function exportPopulationClassification(obj,filepath)
% exportPopulationClassification(obj,filepath) exports the datasetID,
% importID (which is the track number in the import file) and the
% population.
%
% Supported file extensions are: .txt, .dat, .csv, .xls, .xlsm, .xlsx

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
%
% -------------------------------------------------------------
% Copyright (C) 2022 J.J. Erik Maris
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

[~,~,ext] = fileparts(filepath);
if ~ismember(ext,{'.txt','.dat','.csv','.xls','.xlsm','.xlsx'}) % check for correct file format
    error('File format should be either .xls, .xlsm, or .xlsx.')
end

if isfile(filepath) % since the function does not clear the excel sheet before writing, it is good practise to make sure a new file is made
    error(['File ''' filepath ''' does already exist.'])
end


nObj = numel(obj);

% preallocate output tables
T = cell(nObj,1);

for ii = 1:nObj
    % make output table
    T{ii} = obj(ii).fitProps(:,1:2); % datasetID, importID
    Population = repelem(ii,obj(ii).nTracks,1); % population number
    T{ii}.Population = Population;
end

T = vertcat(T{:});
T = sortrows(T,[1 2 3]); % sort rows first datasetID, importID, then population

writetable(T,filepath)







end