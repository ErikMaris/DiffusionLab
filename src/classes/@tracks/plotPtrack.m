function varargout = plotPtrack(obj,ha,Nstep,D_mat,pixeljump_m,loc_err,dt)
% PLOTPTRACK Plot the probability map for full track detection.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 

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

% take the blinking fraction of the first population
% This function need to be documented!

if nargin < 2 || isempty(ha)
    ha = gca;
end

if isempty(obj.coords)
    blinking_fraction = {1};
else
    blinking_fraction = getBlinkingFraction(obj);
end

P_mat = getWorkingRange(D_mat,Nstep,blinking_fraction{1},pixeljump_m,loc_err,dt);

hps = surf(ha,Nstep+1.5,D_mat,P_mat,'EdgeColor','none');
set(gca,'YScale','log');
caxis([0 1])
view(2)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'P_{track}';
set(colorTitleHandle ,'String',titleString);
xlabel('Number of frames')
ylabel('D (m^2/s)')

if nargin > 3
    hold on
    s = {'ko','kx','k>','kd','k+','k*','ks','kv','kp','k<','kh'}; % black symbols
    d = hps;
    hps = gobjects(numel(obj)+1,1);
    hps(1) = d;
    for ii = 1:numel(obj)
        try % D has not been computed
            hps(ii+1) = scatter3(ha,cellfun(@(x) x(end)-x(1) +1,obj(ii).time),...
                obj(ii).Dest.results.D.*obj(ii).getUnitFactor('pixelsize.^2./dt'),...
                repelem(2,obj(ii).nTracks),s{mod(ii,11)+1});
        catch
            warning(['Diffusion constant of population ' num2str(ii) ' cannot be retrieved.'])
        end
    end
    if nargin > 4
        legend_names = cell(numel(obj),1);
        for ii = 1:numel(obj)
            legend_names{ii} = ['Population ' num2str(ii) ];
        end
        legend(hps(2:end),legend_names{:})
    end
end


% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end


function blinking_fraction = getBlinkingFraction(t)
% the blinking fraction is the fraction of all steps (also frames in which
% the molecule is OFF) that belong to a certain blinking gap.
% The algorithm:
% (1) N_step = N_trackpoints - 1
% (2) N_tot_step = N_step + sum over i[N_BG_i * BG_i]
% (3) blinking_fraction = (BG_i + 1) * N_BG_i / N_tot_step


t = t.computeOnTrace;
off = cell(numel(t),1);
blinking_fraction = cell(numel(t),1);
for ii = 1:numel(t)
    [~, off{ii}] = t(ii).getOnOffFromOnTraceShort(true); % get distribution off times
    N_step = sum(t(ii).nTrackPoints - 1); % (1)
    u = unique(off{ii}); % fund unique blinking gaps
    N_tot_step = N_step;
    N_BG = nan(numel(u),1);
    for jj = 1:numel(u) % loop over unique blinking gaps
        N_BG(jj) = sum(off{ii} == u(jj)); % count occurence of each unique blinking gap
        N_tot_step = N_tot_step + N_BG(jj) * u(jj); % (2), add missed steps to compute total number of steps including blinking gaps
    end
    blinking_fraction{ii} = nan(numel(u)+1,1);
    blinking_fraction{ii}(1) = (N_step - sum(N_BG))/N_tot_step; % blinking fraction for blinking gap is 0
    for jj = 1:numel(u)
        blinking_fraction{ii}(jj+1) = ((u(jj) + 1) * N_BG(jj)) / N_tot_step; % (3), blinking fraction for other blinking gaps
    end
end

end


function P_mat = getWorkingRange(D_mat,Nstep,blinking_fraction,pixeljump_m,loc_err_m,dt)
% D_mat: vector with diffusion constants (y-axis)
% N_step: vector with number of displacements (=number of localizations +1
% x-axis)
% blinking_fraction: vector with fraction wrt number of steps corresponding 
% to blinking gap = 0,1,2,3,4,etc
% pixeljump: pixel jump linking
% loc_err: localization error trajectory
% dt: frame time
%
% Idea for optimization: compute large set of displacements and scale the
% pixel jump instead of the diffusion constant

% define analytical solution
Brownian_cdf = @(r,D,dt,n,loc_err) 1-exp(-r.^2./(4.*D.*n.*dt + 4.*loc_err.^2)); % rowland 2017, 2D squared displacement CDF. CDF = linking probability for a certain pixel jump

% normalize fractions
blinking_fraction = blinking_fraction/sum(blinking_fraction); 
nBF = numel(blinking_fraction);

% preallocate
P_mat = ones(numel(D_mat),numel(Nstep)); % preallocate with ones
P_D_plot = nan(numel(D_mat),nBF);

% matrix with n (number of steps; blinking gap + 1)
n_mat = ones(size(P_D_plot)) .* (1:nBF);

% compute linking probability
P_D_plot = Brownian_cdf(pixeljump_m,repmat(D_mat',1,nBF),dt,n_mat,loc_err_m);

% compute track detection probability
for ii = 1:numel(blinking_fraction) % loop over blinking gap + 1
    for jj = 1:numel(Nstep)
        % every loop over ii multiplies the probability with the
        % probability corresponding to that blinking gap. In this way the
        % full track detection probability is computed.
        P_mat(:,jj) = P_D_plot(:,ii).^(blinking_fraction(ii).*Nstep(jj)./ii).*P_mat(:,jj); % P_mat is preallocated with ones
    end
end

end

