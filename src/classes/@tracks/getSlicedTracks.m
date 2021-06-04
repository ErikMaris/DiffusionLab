function [obj,idx] = getSlicedTracks(obj,method,k)
% GETSLICEDTRACKS slices the tracks in shorter segments
%
% [obj,idx] = obj.getSlicedTracks(method,k) returns an object "obj" with 
% sliced tracks and original track index "idx" according to the following 
% methods:
%
%   'uniquepoints': takes subsequent points with a window k, without
%   overlap. Localizations that come in a set smaller than k are discarded.
%
%   'uniquetimes': takes subsequent time with a window k, without
%   overlap. Blinking gaps result in sliced tracks with less than k 
%   localizations. Localizations that come in a set smaller than k are 
%   discarded.

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

% loop over all objects
for ii = 1:numel(obj)
    idx = cell(numel(obj),1);
    switch lower(method)
        case 'uniquepoints'

            nSlicedTracks = floor(obj(ii).nTrackPoints./k); % the number of times k fits in a track
            idx{ii} = repelem(1:obj(ii).nTracks,nSlicedTracks);
            
            % preallocation fitprops
            FP = obj(ii).fitProps;
            FP = sortrows(FP,'trackID'); % sort based on track IDs
            FPheaders = FP.Properties.VariableNames;
            FPtypes = cell(1,numel(FPheaders));
            FPtypes(1) = {'double'};
            FPtypes(2:end) = {'cell'};
            FPsliced = table('Size',[sum(nSlicedTracks) numel(FPheaders)]',...
                'VariableTypes',FPtypes,'VariableNames',FPheaders);
            
            % take only segments in window
            coords = cellfun(@(x,y) x(1:y*k,:),obj(ii).coords,num2cell(nSlicedTracks),'UniformOutput',false); % creates empty if y*k = 0
            time = cellfun(@(x,y) x(1:y*k,:),obj(ii).time,num2cell(nSlicedTracks),'UniformOutput',false);
            
            % merge segments in vector
            coords = vertcat(coords{:});
            time = vertcat(time{:});
            
            % split back into tracks of k size
            coords = mat2cell(coords,repelem(k,sum(nSlicedTracks)),obj(ii).nDim);
            time = mat2cell(time,repelem(k,sum(nSlicedTracks)));
            
            % write
            obj(ii).coords = coords;
            obj(ii).time = time;
            
            % repeat for all columns in table
            FP = table2cell(FP);
            
            for jj = 2:numel(FPheaders) % trackID is first
                temp = cellfun(@(x,y) x(1:y*k,:),FP(:,jj),num2cell(nSlicedTracks),'UniformOutput',false); % creates empty if y*k = 0
                temp = vertcat(temp{:});
                FPsliced(:,jj) = mat2cell(temp,repelem(k,sum(nSlicedTracks)));
            end
            FPsliced(:,1) = num2cell((1:sum(nSlicedTracks))');
                
            obj(ii).fitProps = FPsliced;
            
        case 'uniquetimes'
            warning('Function ''%s'' with method ''%s'' assumes that the time is defined in frames.',mfilename,method)
            
            endtimes = cellfun(@(x) x(end), obj(ii).time);
            
            nSlicedTracks = floor(endtimes./k); % the number of times k fits in a track
            idx{ii} = repelem(1:obj(ii).nTracks,nSlicedTracks);
            
            % preallocation
            coords = cell(sum(nSlicedTracks),1);
            time = cell(sum(nSlicedTracks),1);
            
            % preallocation fitprops cell
            FP = obj(ii).fitProps;
            FP = sortrows(FP,'trackID'); % sort based on track IDs
            FPheaders = FP.Properties.VariableNames;
            nFPheaders = numel(FPheaders);
            % Direct assignement of properties to table leads to conversion
            % from 'double' to 'cell' error
            FPsliced_cell = cell(sum(nSlicedTracks),numel(FPheaders)-1); % do preallocate trackID
            
            % repeat for all columns in table
            FP = table2cell(FP);
            
            cc = 1; % counter
            % So many loops, apologies to my future self... I suppose
            % performance doesn't really matter here.
            for jj = 1:obj(ii).nTracks
                % loop over all segments
                for kk = 1:nSlicedTracks(jj)
                    % make logical with values to take, for k = 5 e.g. [1 2
                    % 3 4 5], [6 7 8 9 10], etc.
                    take = ismembertol(obj(ii).time{jj},((kk-1)*k+1:kk*k)',10^-tracks.TOLERANCE);
                    time{cc} = obj(ii).time{jj}(take);
                    coords{cc} = obj(ii).coords{jj}(take,:);
                    for ll = 2:nFPheaders % trackID is first, skip this one
                        FPsliced_cell{cc,ll-1} = FP{jj,ll}(take,:); 
                    end
                    cc = cc + 1; % update counter
                end
            end
            
            % write
            obj(ii).coords = coords;
            obj(ii).time = time;
            
            % preallocation fitprops table
            FPtypes = cell(1,numel(FPheaders));
            FPtypes(1) = {'double'};
            FPtypes(2:end) = {'cell'};
            FPsliced = table('Size',[sum(nSlicedTracks) numel(FPheaders)]',...
                'VariableTypes',FPtypes,'VariableNames',FPheaders);
                        
            % assign values to table
            FPsliced(:,1) = num2cell((1:sum(nSlicedTracks))');
            for jj = 2:numel(FPheaders) % skip trackID
                FPsliced(:,jj) = FPsliced_cell(:,jj-1);
            end
                
            obj(ii).fitProps = FPsliced;
        otherwise
            error('The method ''%s'' is not recognised by %s.',method,mfilename)
    end
    
    obj(ii).features_valid = false;
    obj(ii).MSD_valid = false;
    obj(ii).Dest_valid = false;
    obj(ii).onTrace_valid = false;
    obj(ii).distOnOff_valid = false;
    
end


end




