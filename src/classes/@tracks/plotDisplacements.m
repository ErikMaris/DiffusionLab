function plotDisplacements(obj,noCovarFlag)
% noCovarFlag defines wether localizations may occur in multiple
% displacements

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

colours = cmapMCEC;

X = cell(obj.nTracks,1);

if noCovarFlag
    for ii = 1:obj.nTracks
        X{ii} = obj.coords{ii}(2:end,:) - obj.coords{ii}(1:end-1,:);
    end
else
    for ii = 1:obj.nTracks
        X{ii} = obj.coords{ii}(2:2:end,:) - obj.coords{ii}(1:2:end-1,:);
    end
end

X = vertcat(X{:});

[N,Xedges] = histcounts(X(:,1),200);
Xvals = Xedges(1:end-1) + 0.5.*diff(Xedges);
f = fit(Xvals(:),N(:),'gauss1');


subplot(1,2,1)
histogram(X(:,1),Xedges,'FaceColor',colours(1,:))
hold on
histogram(X(:,1),Xedges,'FaceColor',colours(2,:))
plot(f)
hold off
subplot(1,2,2)
scatter(X(:,1),X(:,2))
axis square


disp(f)


end

