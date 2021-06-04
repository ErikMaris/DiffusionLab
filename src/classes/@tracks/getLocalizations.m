function obj_loc = getLocalizations(obj)
% getLocalizations Gets a localizations object from tracks
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj_loc = getLocalizations(obj)

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

% if nargin < 2 || isempty(indices)
%     indices = 1 : obj.nTracks;
% end

nii = numel(obj);

% --- run linking algorithms
obj_loc = localizations.empty(nii,0);
for ii = 1:nii
    if obj(ii).verbose == 1; disp(['Caluclating track ' num2str(ii) ' of ' num2str(nii) '.']); end
    
    allcoords = [vertcat(obj(ii).time{:}) vertcat(obj(ii).coords{:})];
    allcoords = sortrows(allcoords,1);
    allcoords(:,1) = tracks.roundn(allcoords(:,1), tracks.TOLERANCE); % round to allow grouping
    
%     First_idx =  find(abs(time-circshift(time,1)) > 10^-tracks.TOLERANCE); % find(DT-circshift(DT,1) ~= 0);
%     Last_idx = find(abs(time-circshift(time,-1)) > 10^-tracks.TOLERANCE); % find(DT-circshift(DT,-1) ~= 0)
    
    % --- store into general format
    dx = diff(allcoords(:,1))~=0;
    startIdx = [1 ;find(dx ~= 0)+1; size(dx,1) + 2]; % + 2 because startIdx should be padded 1 at begin and end of vector
    dev = diff(startIdx);
    
    % --- convert output into tracks structure
    coords = mat2cell(allcoords(:,2:end),dev,obj(ii).nDim);
    
    time = obj(ii).getCommonTimes;
    
%     for jj = 1:numel(indices)
%         
%         kk = indices(jj);
%         
%         [~, index_in_all_delays, ~] = intersect(time, obj(ii).time{kk});
%         
%         
%         coords{index_in_all_delays} = vertcat(coords{index_in_all_delays},obj(ii).coords{kk});
%         
%     end
        
    
    % --- assign to localization object

    obj_loc(ii)  = localizations(time,coords,obj(ii).getUnitFactor('m'),obj(ii).getUnitFactor('s'));
    
end

end

