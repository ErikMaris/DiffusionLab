function varargout = plotDistOnOff(obj, plotOnFlag, varargin)

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

% --- parse input

defaultNormalization = 'count';
expectedNormalization = {'count','countdensity','cumcount','probability','pdf','cdf'};
defaultAxes = 'Semilogy';
expectedTypeAxes = {'Linear','Semilogx','Semilogy','Loglog'};
defaultHa = gca;
defaultIndices = 1:obj.nTracks;
defaultEdges = [];

p = inputParser;

addRequired(p,'obj',@(x) isobject(x));
addRequired(p,'plotOnFlag',@(x) islogical(x) || x == 0 || x == 1);
addOptional(p,'ha',defaultHa,@isobject);
addOptional(p,'indices',defaultIndices,@(x) isnumeric(x) && isvector(x));
addParameter(p,'Normalization',defaultNormalization,@(x) any(validatestring(x,expectedNormalization)));
addParameter(p,'Axes',defaultAxes,@(x) any(validatestring(x,expectedTypeAxes)));
addParameter(p,'Edges',defaultEdges,@(x) isnumeric(x) && isvector(x));
parse(p,obj,plotOnFlag, varargin{:});

ha = p.Results.ha;
indices = p.Results.indices;
Normalization = p.Results.Normalization;
Axes = p.Results.Axes;
Edges = p.Results.Edges;

clear p

if numel(obj) > 1
    error('SC_blinking:multipleObjParsed', ...
        'Function does not accept object arrays.')
end

% --- get on and off values
if plotOnFlag
    [histvals,~] = obj.getOnOffFromOnTraceShort(true,indices);
else
    [~,histvals] = obj.getOnOffFromOnTraceShort(true,indices);
end

tunit = 's';
[unitFactor,unitLabel] = obj.getUnitFactor(tunit);
histvals = histvals.*unitFactor;

colors = cmapMCEC(1);

if ~isempty(Edges)
    hps = histogram(ha,histvals,Edges,'Normalization',Normalization,'FaceColor',colors);
else
    if strcmp(Axes,{'Semilogx','Loglog'})
        % devide axes
        [~,Xedges] = histcounts(log10(histvals)); % in log
        hps = histogram(ha,histvals,10.^Xedges,'Normalization',Normalization,'FaceColor',colors);
    else
        hps = histogram(ha,histvals,'Normalization',Normalization,'FaceColor',colors);
    end
end



xlabel(['Length (' unitLabel ')'])
ylabel(Normalization)

switch(Axes)
    case 'Linear'
        set(ha, 'XScale', 'linear')
        set(ha, 'YScale', 'linear')
    case 'Semilogx'
        set(ha, 'XScale', 'log')
        set(ha, 'YScale', 'linear')
    case 'Semilogy'
        set(ha, 'XScale', 'linear')
        set(ha, 'YScale', 'log')
    case 'Loglog'
    	set(ha, 'XScale', 'log')
        set(ha, 'YScale', 'log')
end


if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

