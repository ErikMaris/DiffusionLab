function reportDest(obj)
% reportDest Reports the results of the diffusion estimator (Dest).
%
% reportDest(obj) reports the results of the diffusion estimator contained 
%   in tracks.Dest.

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

% fprintf('--------------------------------------------------------\n\n')
for ii = 1:numel(obj)
    for jj = 1:numel(obj(ii).Dest)
        fprintf('Population %d, diffusion estimator %d.\n',ii,jj)
        obj(ii).Dest(jj).DModel.report(obj(ii).Dest(jj),obj(ii))
    end
    fprintf('\n')
end

% fprintf('--------------------------------------------------------\n')


end