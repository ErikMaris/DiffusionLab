function varargout = plotOnTrace(obj, ha, index)
% plotOnTrace Plot the On Trace of one track.
%
% obj.plotOnTrace plots the On Trace in the current axes.
% 
% obj.plotOnTrace(ha, index) plots the On Trace of track with index in the 
% axes specified by the handle ha.
%
% hps =  obj.plotOnTrace(...) returns the handle array for the
% lines generated.
%
% [hps, ha] =  obj.plotOnTrace(...) also return the axes handle in
% which the lines were plotted.


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

if numel(index) ~= 1 || ~isnumeric(index)
    error('plotOnTrace:badInput:ID',...
        'Track ID %f not recognized',index)
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

nObj = numel(obj);
hp = gobjects(nObj,2);

for ii = 1:nObj
    
    colors = obj(ii).getColormap(nObj);
    
    trackName = sprintf('Track %d, population %d', index, ii );
    
    onTraceFull = getOnTraceFromShort(obj(ii),obj(ii).onTrace{index});
    
    X = 1:obj(ii).nFrames;
    
    hp(ii,1) = stairs(ha, X, onTraceFull, ...
        'Color', colors(ii,:), ...
        'LineWidth', obj(ii).lineWidth, ...
        'DisplayName', trackName );
    
    hold(ha, 'on');
    
    X2 = [X;X];
    Y = [onTraceFull';onTraceFull'];
    hp(ii,2) = area(X2([2:end end]),Y(1:end),'FaceColor',colors(ii,:));
    
end

hp = hp(:); 

xlabel('Frame')
ylabel('ON?')

if nargout > 0
    varargout{1} = hp;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

