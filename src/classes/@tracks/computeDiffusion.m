function obj = computeDiffusion(obj)
% computeDiffusion Wrapper to compute the Dest class from the parent
% 
% obj = computeDiffusion(obj)

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


for ii = 1:numel(obj)
    if isempty(obj(ii).Dest)
        warning('tracks:computeDiffusion:noEstimator',...
            'No estimator loaded for population %i.',ii)
        continue
    end
    obj(ii).Dest = obj(ii).Dest.compute(obj(ii));
    % load all results into features table
    try
        fn = obj(ii).Dest.results.Properties.VariableNames;
        if isempty(obj(ii).features)
            obj(ii) = obj(ii).computeFeatures(false,'nTrackPoints');
        end
        clear = ismember(obj(ii).features.Properties.VariableNames,fn);
        obj(ii).features(:,clear) = []; % clear duplicates and prevent error
        obj(ii).features = [obj(ii).features obj(ii).Dest.results];
        disp('Diffusion features loaded.')
    catch
        % do nothing
    end
    % set to valid
    obj(ii).Dest_valid = true;
end

end