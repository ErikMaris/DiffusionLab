function time = getCommonTimes(obj)
% getCommonTimes finds the existing time tags in data set
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% time = getCommonTimes(obj)

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from msdanalyzer https://tinevez.github.io/msdanalyzer/ with
% permission.
% (C) 2007 Jean-Yves Tinevez
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

time = vertcat(obj.time{:});
time = unique(tracks.roundn(time, tracks.TOLERANCE));