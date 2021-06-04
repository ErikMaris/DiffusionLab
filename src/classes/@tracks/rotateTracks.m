function obj = rotateTracks(obj,angle,centerFlag)
% ROTATETRACKS Rotates the coordinates by a certain angle.
%
% obj = rotateTracks(obj,angle) returns a tracks object with the
% coordinates rotated by an angle specified by 'angle'. This means:
%   2D: angle is [theta]
% The angle can  either be a scalar or row vector, and the same rotation is
% applied to all objects, or a column vector or matrix [number_objects 
% angle] to apply a specific rotation to each object.
%
% obj = rotateTracks(obj,angle,centerFlag) centers the rotation around the
% minimum bounding circle center

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

if nargin < 3
    centerFlag = true;
end

% format angle
if size(angle,1) == 1
    % repeat angle for all objects
    angle = repmat(angle,[numel(obj),1]);
elseif size(angle,1) == numel(obj)
    % do nothing
else
    error('Angle does not have the expected dimensions. Please read the help.')
end
    
    

for ii = 1:numel(obj)
    
    if centerFlag
        % get coordinates center minimum bounding circle
        allCoords = vertcat(obj(ii).coords{:});
        [center,~] = tracks.minboundcircle(allCoords(:,1),allCoords(:,2),true);
        center = center(:);
        clear allCoords
    else
        center = [0;0];
    end

    if obj(ii).nDim == 2
        theta = angle(ii,1);
        % define rotation matrix
        % https://nl.mathworks.com/matlabcentral/answers/93554-how-can-i-
        % rotate-a-set-of-points-in-a-plane-by-a-certain-angle-about-an-arbitrary-point
        R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
        % loop over all tracks
        for jj = 1:obj(ii).nTracks
            obj(ii).coords{jj} = (R*(obj(ii).coords{jj}' - center) + center)';
        end
    else
        error('Only 2 dimensional data is supported.')
    end
    
    obj(ii).features_valid = false;
    obj(ii).Dest_valid = false;
    
end

end

