function exportMeanMSD(obj,filepath,indices)
% exportMeanMSD(obj,filepath, indices)

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

if nargin < 3
    indices = [];
end

nObj = numel(obj);

warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ; % if multiple sheets have to be written
for ii = 1:nObj

    [xUnitFactor,time_unit] = obj(ii).getUnitFactor('dt');
    [yUnitFactor,msd_unit] = obj(ii).getUnitFactor('pixelsize.^2');

    msmsd = obj(ii).getMeanMSD(indices);
    
    if isempty( msmsd ) % do not plot empty spots in MSD
        continue
    end

    t = msmsd(:,1).*xUnitFactor;
    m = msmsd(:,2).*yUnitFactor;
    s = msmsd(:,3).*yUnitFactor; % standard error of the mean
    msd_out = [NaN NaN NaN; t m s];
    
    writematrix(msd_out,filepath,'Sheet',['Population ' num2str(ii)]); % make file
    writecell({['Delay time (' time_unit ')'],['MSD (' msd_unit ')'],['Standard error (' msd_unit ')']},filepath,'Sheet',['Population ' num2str(ii)],'Range','A1:C1')
    
end
warning( 'on', 'MATLAB:xlswrite:AddSheet' ) ;

end