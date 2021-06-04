function dR2 = static_getdR2(t,X,delayTimes)
% static_getdR2 Computes the squared displacements
%
% dR2 = static_getdR2(time, coords) returns a matrix 'dR2' with [delays 
% squared_displacement] for a vector 'time' and 'coords' [X Y ..], and a
% delayTime.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% 

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Adapted from Maxime Deforet, May 21 2013
% https://nl.mathworks.com/matlabcentral/fileexchange/
% 41858-kehl-a-fast-no-loop-method-to-compute-msd
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

n_detections = size(X, 1);

% This code contains no loop. Each step is "vectorized", therefore pretty
% fast.
% The idea is to compute all the possible pairings in only one time.
% Trajectory = list of T positions (x,y,t). Or (x,t) or (x,y,z,t) or even
% higher dimension (x1,x2,x3,...t)
% tau is the list of all the possible time intervals within the trajectory.
% MSD is a list of mean squared displacements, for each value of tau.
% MSD(tau) = sqrt( x(i)^2 - x(i+tau)^2 ), averaged over i.
[I, j] = find(triu(ones(n_detections),1)); % list of indices of possible pairings: I row; j col
D = zeros(n_detections, n_detections);
D(I + n_detections*(j-1)) = (sum(abs(X(I,:) - X(j,:)).^2,2)); % Squared distance computation in one line !
% Time intervals between the two points of each pairing :
dt = zeros(n_detections, n_detections);
dt(I + n_detections*(j-1)) = -(t(I) - t(j));
% Then D is a list of squared distances. dt is a list of corresponding
% time intervals. Now we have to average all D values over different dt
% values
% We first remove all 0 values from dt matrix, and from D as well.
idx_0 = find(dt == 0);
dt(idx_0) = [];
D(idx_0) = [];
% Then we sort dt in ascending order, and sort D in the same way.
[DT,idx] = sort(dt(:));
DD = D(idx);
% We now have DD and DT, respectively a list of squared distances, and
% the corresponding time intervals.

% store 
DT = tracks.roundn(DT, tracks.TOLERANCE); % this is important, otherwise identical delays will not be recognised
%     take_idx = ismembertol(DT,delayTimes,10^-tracks.TOLERANCE);
take_idx = ismembertol(DT,delayTimes);
dR2 = [DT(take_idx) DD(take_idx)'];


end

