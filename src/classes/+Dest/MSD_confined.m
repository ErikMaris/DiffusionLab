classdef MSD_confined < Dest.DModel
    %NORMAL Model for confined diffusion
    
    % ---------------------------------------------------------------------
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
    
    % general
    properties (Constant)
        name = 'MSD_confined';
    end

    % MSD
    properties (Constant)
        fitfunction = 'p1*(1-exp(-x/p2))';
        lowerBounds = [0 0];
        upperBounds = [Inf Inf];
        nCoeff = 2;
    end
    
    methods (Static)
        function Dest = computeResults(Dest,~)
            Dest.results = table(Dest.coeffvals(:,1),Dest.coeffvals(:,2));
            Dest.results.Properties.VariableNames = {'r2','tau'};
            Dest.results.Properties.VariableDescriptions = {'r2','tau'};
            Dest.results.Properties.VariableUnits = {'pixelsize.^2','dt'};
        end
        
        function Dest = computeMeanResults(Dest,~)
            Dest.mResults = struct('r2',[],'tau',[],'r2Sigma',[],'tauSigma',[]);
            Dest.mResults.R2 = Dest.mCoeffvals(:,1);
            Dest.mResults.tau = Dest.mCoeffvals(:,2);
            
            % compute errors
            try % if enough data points for error estimation are present
                CI = confint(Dest.mFit_obj,0.6827); % https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
                CI = (CI(2,:) - CI(1,:))./2;
                Dest.mResults.R2Sigma = CI(1);
                Dest.mResults.tauSigma = CI(2);
            catch e
                warning(e.message) % retrow error as warning
                Dest.mResults.R2Sigma = nan; % prevent error in reportDest
                Dest.mResults.tauSigma = nan;
            end
        end
        
        function StartPoints = getStartPoints(x,y)
            % get start points for fit y=p1*(1-exp(x/p2))
            
            n = numel(y);
            n = ceil(n/2);
            StartPoints = [mean(y(n:end)) x(n)]; % r^2 = mean second part MSD; tau = 0
            
        end
        
        function report(Dest,parent)
            
            [uF_R2,uL_R2] = parent.getUnitFactor('pixelsize.^2');
            [uF_tau,uL_tau] = parent.getUnitFactor('dt');
            fprintf('Mean squared displacement analysis for confined diffusion.\n')
            fprintf('The mean confinement parameter R^2 is %e %s %e %s.\n',Dest.mResults.R2.*uF_R2,char(177),Dest.mResults.R2Sigma.*uF_R2,uL_R2)
            fprintf('The mean tau is %e %s %e %s.\n',Dest.mResults.tau.*uF_tau,char(177),Dest.mResults.tauSigma.*uF_tau,uL_tau)
            
        end
        
    end
end