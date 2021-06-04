function [pd, varargout] = getLocErrFromImmobile(obj,movmean_window,indices)

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

pd = cell(numel(obj),1);
ha = cell(numel(obj),1);

for ii = 1:numel(obj)

    if nargin < 3
        indices = 1:obj(ii).nTracks;
    end
    
    for jj = 1:numel(indices)

        ha{ii} = gobjects(numel(indices),2);

        index = indices(jj);
        
        t = obj(ii).time{index};
        t = tracks.roundn(t, tracks.TOLERANCE);
        X = obj(ii).coords{index};  
        
        if movmean_window ~= false

            obj = obj.computeDrift('movmean',index,movmean_window);

            tdrift = obj(ii).drift(:,1);
            xdrift = obj(ii).drift(:,2:end);

            % Determine target delay index in bulk
            % find to-be-corrected time indices in track
            [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
            % Keep only track times that can be corrected.
            X = X(index_in_track_time, :);
            % Subtract drift position to track position
            X = X - xdrift(index_in_drift_time, :);
            
        end


%                     colors = obj(ii).getColormap(numel(obj));
%                     figure
%                     subplot(1,2,1)
%                     scatter(obj(ii).coords{trackID}(:,1),obj(ii).coords{trackID}(:,2),'.','MarkerEdgeColor',colors(1,:))
%                     axis image
%                     box on
%                     subplot(1,2,2)
%                     scatter(tcorr(ii).coords{trackID}(:,1),tcorr(ii).coords{trackID}(:,2),'.','MarkerEdgeColor',colors(1,:))
%                     axis image
%                     box on

        X = X - mean(X,1); % error coordinates, i.e. true position at X = (0,0)

%                     ft = fittype('a1*exp(-0.5*((x-b1)/c1)^2)','coefficients',{'a1','b1','c1'}); % normal distribution

        uf = obj(ii).getUnitFactor('nm');

        X = X.*uf; % convert to nm

%                     [Nx,edgesx] = histcounts(X(:,1));
%                     [Ny,edgesy] = histcounts(X(:,2));
%                     Xv = edgesx(1:end-1)'+0.5*diff(edgesx)';
%                     Yv = edgesy(1:end-1)'+0.5*diff(edgesy)';
%                     fox = fit(Xv,Nx',ft,'StartPoint',[max(Nx),0,10]);
%                     foy = fit(Yv,Ny',ft,'StartPoint',[max(Ny),0,10]);
        pdx = fitdist(X(:,1),'Normal');
        pdy = fitdist(X(:,2),'Normal');
%         figure
%         ha{ii}(jj,:) = histfit(X(:,1));
%         xlabel('Displacement error (nm)')
%         ylabel('Probability')

%                     figure
%                     hold on
%                     for ii = 1:numel(tcorr)
%                         histogram(errory{ii}(:,1).*uf,edgesy,'Normalization','pdf',"FaceColor",colors(ii,:))
%                         hold on
%                         plot(Xv,pdf(pdy{ii},Xv),'LineWidth',2,'Color','r')
%                     end
%                     xlabel('x - x_{mean} (nm)')
%                     ylabel('Probability')

        pd{ii}(jj,:) = [pdx.sigma pdx.mu pdy.sigma pdy.mu];
    end

end

if nargin > 1
    varargout{1} = ha;
end

end