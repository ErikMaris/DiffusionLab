function exportMSD(obj,filepath)

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

if any(~[obj.MSD_valid])
    error('Please compute the MSD first.')
end

[~,~,ext] = fileparts(filepath);
if ~ismember(ext,{'.xls','.xlsm','.xlsx'}) % check for correct file format
    error('File format should be either .xls, .xlsm, or .xlsx.')
end

if isfile(filepath) % since the function does not clear the excel sheet before writing, it is good practise to make sure a new file is made
    error(['File ''' filepath ''' does already exist.'])
end

nObj = numel(obj);

warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ; % if multiple sheets have to be written
for ii = 1:nObj
    allDelays = obj(ii).getAllDelays;
    msd_out = nan(numel(allDelays)+1,obj(ii).nTracks+1); % preallocate with nan; nan gives an empty cell in excel sheet
    msd_out(2:end,1) = allDelays .* obj(ii).getUnitFactor('dt'); % add delay times
    msd_out(1,2:end) = 1:obj(ii).nTracks; % add track index
    std_out = msd_out; % same headers
    
    [~,time_unit] = obj(ii).getUnitFactor('dt'); % get time unit name
    [uF_MSD,msd_unit] = obj(ii).getUnitFactor('pixelsize.^2');
    for jj = 1:obj(ii).nTracks
        [~, idx_allDelays, ~] = intersect(allDelays, obj(ii).MSD{jj}(:,1));
        msd_out(idx_allDelays+1,jj+1) = obj(ii).MSD{jj}(:,2).*uF_MSD; % +1 because we skip the header line
        std_out(idx_allDelays+1,jj+1) = obj(ii).MSD{jj}(:,3).*uF_MSD; % +1 because we skip the header line
    end
    writematrix(msd_out,filepath,'Sheet',['Population ' num2str(ii)]); % make file
    writecell({['Delay time (' time_unit ') \ MSD (' msd_unit ')']},filepath,'Sheet',['Population ' num2str(ii)],'Range','A1')
    writematrix(std_out,filepath,'Sheet',['Population ' num2str(ii)],'Range',['A' num2str(size(msd_out,1)+2)]); % make file
    writecell({['Delay time (' time_unit ') \ standard deviation MSD (' msd_unit ')']},filepath,'Sheet',['Population ' num2str(ii)],'Range',['A' num2str(size(msd_out,1)+2)])
end
warning( 'on', 'MATLAB:xlswrite:AddSheet' ) ;

end