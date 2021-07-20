function msmsd = getMeanMSD(obj, indices)
% GETMEANMSD Compute the weighted mean of all MSD curves.
%
% Linear fit is the MSD of all tracks
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% msd = obj.getMeanMSD computes and return the weighted mean of all
% MSD curves stored in this object. All possible delays are first
% derived, and for each delay, a weighted mean is computed from all
% the MSD curves stored in this object. Weights are set to be the
% number of points averaged to generate the mean square
% displacement value at the given delay. Thus, we give more weight
% to MSD curves with greater certainty (larger number of elements
% averaged).
%
% Results are returned as a N x 5 double array, and ordered as
% following: [ dT M STD N W] with:
% - dT the delay vector
% - M the weighted mean of MSD for each delay
% - STD the weighted standard deviation
% - N the number of degrees of freedom in the weighted mean
% - W the number of data points contributing to mean
% (see http://en.wikipedia.org/wiki/Weighted_mean)
%
% msd = obj.getMeanMSD(indices) only takes into account the MSD
% curves with the specified indices.

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

% Adapted and checked from msdanalyzer - Erik Maris 2019/04/02

if numel(obj) > 1
    error('Can only process one tracks object at the time.')
end

if ~obj.MSD_valid
    obj = obj.computeMSD(indices);
end

if nargin < 2 || isempty(indices)
    indices = 1 : obj.nTracks;
end

MSD = obj.MSD(indices); % get desired MSD's

msmsd = tracks.static_getMeanMSD(MSD);

end