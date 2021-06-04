function [dR2,iCPD,dt] = getR2CPD(obj,indices,delayIdx)

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

if nargin < 3 || isempty(delayIdx)
    delayIdx = 1;
end

nObj = numel(obj);

if nObj > 1 && numel(delayIdx) > 1
    error('tracks:getR2CPD:multipleObjDelayIdx',...
        'Cannot get multiple delay indices for multiple objects.')
end

all_indices = false;
if nargin < 3 || isempty(indices)
    all_indices = true;
end

if nObj > 1 % plot multiple populations; get cell with this
    dR2 = cell(nObj,1);
    for ii = 1:nObj
        if all_indices
            indices = 1:obj(ii).nTracks;
        end
        dR2{ii} = obj(ii).getdR2(delayIdx,indices);
        dR2{ii} = dR2{ii}{1}; % only one delayIdx
    end
else % plot multiple delay times; get cell with this
    if all_indices
        indices = 1:obj.nTracks;
    end
    [dR2, dt] = obj.getdR2(delayIdx,indices);
end

% data has already been sorted in getdR2 from low to high dR2
% compute the cumsum...
iCPD = cellfun(@(x) linspace(1,0,numel(x))', dR2,'UniformOutput',false);
% ...and normalize to get the cpd. Compute the inverse cpd for plotting
% iCPD = cellfun(@(x) 1-x./x(end), iCPD,'UniformOutput',false);

end

