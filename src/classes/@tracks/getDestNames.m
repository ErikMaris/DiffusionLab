function names = getDestNames(folderpath)
% getDestNames Gets the names of the diffusion estimators in folder
% 
% Escapes files DestModel.m, obsolete, and CDFfuncfree.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% 
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

if nargin < 2
    folderpath = mfilename('fullpath');
    [folderpath,~,~] = fileparts(folderpath);
    folderpath = [folderpath filesep '..' filesep '+Dest']; % default
end

listing = dir(folderpath);
names = {listing.name};

filter = [zeros(1,2) ones(1,numel(names)-2)];
filter = filter .* ~strcmp(names,{'DestModel.m'});
filter = filter .* ~strcmp(names,{'obsolete'});
filter = filter .* ~strcmp(names,{'CDFfuncfree'});
names = names(find(filter == 1));

if isempty(names)
    return
end

for ii = 1:numel(names)
    fn = names{ii};
    N = numel(fn);
    names{ii} = fn(1:N-2);
end

end

