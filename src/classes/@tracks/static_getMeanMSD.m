function msmsd = static_getMeanMSD(MSD)
% STATIC_GETMEANMSD Compute the weighted mean of all MSD curves
%
% msmsd = static_getMeanMSD(MSD) returns the mean MSD in 'msmsd'
% as [ T mean std Nfreedom N] for a cell with MSDs grouped per track 
% [delays mean_msd std_msd n_msd] as output by 'static_getMSD',
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
% Reproduced from msdanalyzer https://tinevez.github.io/msdanalyzer/
% (C) 2014 Jean-Yves Tinevez
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

% MATLAB only makes a physical copy MSD if we alter it, so as long as we
% read-only, using this function with a wrapper should not cost extra
% memory.

n_tracks = numel(MSD);

% First, collect all possible delays; faster to get them from tracks than
% to calculate them via getAllDelays (?)
all_delays = cell(n_tracks, 1);
for ii = 1 : n_tracks
    if isempty( MSD{ii} )
        continue
    end
    all_delays{ii} = MSD{ii}(:,1); % get delays
end
delays = unique( vertcat( all_delays{:} ) );
n_delays = numel(delays);

% Collect
sum_weight          = zeros(n_delays, 1); % sum number of displacements contributing
sum_weighted_mean   = zeros(n_delays, 1); % sum MSD * number of displacements

% 1st pass
for ii = 1 : n_tracks
    
    if isempty( MSD{ii} )
        continue
    end
    
    t = MSD{ii}(:,1);
    m = MSD{ii}(:,2);
    n = MSD{ii}(:,4);
    
    % Do not take NaNs
    valid = ~isnan(m);
    t = t(valid);
    m = m(valid);
    n = n(valid);
    
    % Find common indices
    [~, index_in_all_delays, ~] = intersect(delays, t);
    
    % Accumulate
    sum_weight(index_in_all_delays)           = sum_weight(index_in_all_delays)         + n; % V1
    sum_weighted_mean(index_in_all_delays)    = sum_weighted_mean(index_in_all_delays)  + m .* n; % MSD * n
end

% Compute weighted mean
mmean = sum_weighted_mean ./ sum_weight; % ensemble mean

% 2nd pass: unbiased variance estimator
sum_weighted_variance = zeros(n_delays, 1);
sum_square_weight     = zeros(n_delays, 1);

for ii = 1 : n_tracks
    
    if isempty( MSD{ii} )
        continue
    end
    
    t = MSD{ii}(:,1);
    m = MSD{ii}(:,2);
    n = MSD{ii}(:,4);
    
    % Do not take NaNs
    valid = ~isnan(m);
    t = t(valid);
    m = m(valid);
    n = n(valid);
    
    % Find common indices
    [~, index_in_all_delays, ~] = intersect(delays, t);
    
    % Accumulate "Weighted sample variance" https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
    sum_weighted_variance(index_in_all_delays)    = sum_weighted_variance(index_in_all_delays)  + n .* (m - mmean(index_in_all_delays)).^2 ; % w*(msd - mean_msd)^2 = nominator
    sum_square_weight(index_in_all_delays)        = sum_square_weight(index_in_all_delays)      + n.^2; % V2
end

% Standard deviation
mstd = sqrt( sum_weight ./ (sum_weight.^2 - sum_square_weight) .* sum_weighted_variance ); 
% sqrt( sum_[n] / ((sum_[n]^2 - sum[n^2]) * (msd - mean_msd)^2 ) )
% V1 * nominator / (V1^2 - V2) = nominator/ (V1 - (V2/V1)) = unbiased
% estimate of sample variance for "Reliability weights" according to
% Wikipedia

% If in the rare case, the "sum_weighted_variance" is zero, this yields a
% mstd of NaN. We set the mstd to Inf for the weights calculation in the 
% mean MSD fit.
mstd(isnan(mstd)) = Inf;

% Output [ T mean std Nfreedom N]
msmsd = [ delays mmean mstd (sum_weight.^2 ./ sum_square_weight) sum_weight];

end