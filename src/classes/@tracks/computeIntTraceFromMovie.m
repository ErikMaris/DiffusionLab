function obj = computeIntTraceFromMovie(obj,moviepath,ROI_px)
% computes the mean intensity trace from a movie based on the center of 
% mass of the tracks and the edge length given by ROI_px.

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

% --- check whether blinking superclass exists
if ~strcmp(superclasses(obj),'blinking')
    return
end

% I = bfopen(moviepath);
% I = I{1}(:,1);
I = Localizer('readccdimages', 1, -1, moviepath);
% I = permute(I,[2 1 3]);

ROI_px = floor((ROI_px-1)/2); % get ROI size in -ROI_px to + ROI_px

for ii = 1:numel(obj)
    if ~obj.center_valid
        % --- recompute if necessary
        obj(ii) = obj(ii).computeCenter;
    end
    
    obj(ii).intTrace = cell(obj(ii).nTracks,1);
    CoM = round(obj(ii).CenterOfMass);
    for jj = 1:obj(ii).nTracks
        obj(ii).intTrace{jj} = squeeze(mean(mean(I(CoM(jj,1)-ROI_px:CoM(jj,1)+ROI_px,CoM(jj,2)-ROI_px:CoM(jj,2)+ROI_px,:),1),2));
    end
    
end

end

