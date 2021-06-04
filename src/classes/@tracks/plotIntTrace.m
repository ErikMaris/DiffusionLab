function varargout = plotIntTrace(obj, ha, indices)
% plotIntTrace Plots the intensity trace.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj.plotIntTrace plots the intensity trace in the current axes.
%
% obj.plotIntTrace(ha) plots the intensity trace in the axes
% specified by the handle ha.
%
% obj.plotIntTrace(ha, indices) plots the intensity trace for the
% particles with the specified indices only. Leave empty to
% plot for all particles.
%
% hps =  obj.plotIntTrace(...) returns the handle array for the
% lines generated.
%
% [hps, ha] =  obj.plotIntTrace(...) also return the axes handle in
% which the lines were plotted.

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

if nargin < 2
    ha = gca;
end
if nargin < 3 || isempty(indices)
    indices = 1 : obj.nTracks;
end

xunit = 's';
[xUnitFactor,xUnitLabel] = obj.getUnitFactor(xunit);

n_spots = numel(indices);
colors = jet(n_spots);

hps = NaN(n_spots, 1);

hold(ha, 'on');

for ii = 1 : n_spots
    
    index = indices(ii);
    
    traceName = sprintf('Trace %d', index );
    
    t = 1:length(obj.intTrace{index});
    t = t'.* xUnitFactor;
    int = obj.intTrace{index};

    hps(ii) = plot(ha, t, int, ...
        'Color', colors(ii,:), ...
        'DisplayName', traceName );
    
end

xlabel(['time (' xUnitLabel ')'])
ylabel('Mean counts')

hold(ha, 'off');

if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

