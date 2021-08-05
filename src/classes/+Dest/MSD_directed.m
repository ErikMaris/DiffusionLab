classdef MSD_directed < Dest.DModel
    %NORMAL Model for directed diffusion
    
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
        name = 'MSD_directed';
    end
    
    % MSD
    properties (Constant)
        fitfunction = 'p1*x+p2 + p3*x^2';
        lowerBounds = [-Inf -Inf 0];
        upperBounds = [Inf Inf Inf];
        nCoeff = 3;
    end
    
    methods (Static)
        function Dest = computeResults(Dest,parent)
            % computes the derived values and stores in estimator object
            
            D = Dest.coeffvals(:,1) ./ 2 ./ parent.nDim; % Eq. B1
            locVar = Dest.coeffvals(:,2) ./ 2 ./ parent.nDim + 2 .* parent.getR .*D .* parent.dt; % Eq. B1
            diffusionSNR = real(sqrt(D.*parent.dt./locVar)); % definition Vestergaard
            diffusionSNR(D < 0 | locVar < 0) = nan;
            velocity = sqrt(Dest.coeffvals(:,3));
            
            Dest.results = table(D,locVar,diffusionSNR,velocity);
            Dest.results.Properties.VariableDescriptions = {'Diffusion constant','Localization variance','Diffusion SNR','Velocity'};
            Dest.results.Properties.VariableUnits = {'pixelsize.^2/dt','pixelsize.^2','','pixelsize./dt'};

        end
        
        function Dest = computeMeanResults(Dest,parent)
            % computes the derived values and stores in estimator object
            
            Dest.mResults = struct('D',[],'locVar',[],'velocity',[],'Dsigma',[],'locVarSigma',[],'velocitySigma',[]);
            
            Dest.mResults.D = Dest.mCoeffvals(1) ./ 2 ./ parent.nDim; % Eq. B1
            Dest.mResults.locVar = Dest.mCoeffvals(2) ./ 2 ./ parent.nDim + 2 .* parent.getR .*Dest.mResults.D .* parent.dt; % Eq. B1
            Dest.mResults.velocity = sqrt(Dest.mCoeffvals(:,3));
            
            % compute errors
            try % if enough data points for error estimation are present
                CI = confint(Dest.mFit_obj,0.6827); % https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
                CI = (CI(2,:) - CI(1,:))./2;
                Dest.mResults.Dsigma = CI(1) ./ 2 ./ parent.nDim; % Eq. B1
                Dest.mResults.locVarSigma = CI(2) ./ 2 ./ parent.nDim + 2 .* parent.getR .*Dest.mResults.D .* parent.dt; % Eq. B1
                Dest.mResults.velocitySigma = sqrt(CI(3));
            catch e
                warning(e.message) % retrow error as warning
                Dest.mResults.Dsigma = nan; % prevent error in reportDest
                Dest.mResults.locVarSigma = nan;
                Dest.mResults.velocitySigma = nan;
            end
        end
        
        function StartPoints = getStartPoints(x,y)
            % get start points for fit y=p1*x+p2 + p3*x^2
            
            StartPoints = [mean(y)/mean(x) y(1) 0];
            
        end
        
        function report(Dest,parent)
            
            [uF_D,uL_D] = parent.getUnitFactor('pixelsize.^2./dt');
            [uF_v,uL_v] = parent.getUnitFactor('pixelsize./dt');
            fprintf('Mean squared displacement analysis for directed diffusion.\n')
            fprintf('The mean diffusion constant is %e %s %e %s.\n',Dest.mResults.D.*uF_D,char(177),Dest.mResults.Dsigma.*uF_D,uL_D)
            if Dest.mResults.locVar < 0
                [uF_LV,uL_LV] = parent.getUnitFactor('pixelsize.^2');
                fprintf('The mean localization error is < 0, giving the localization variance: %e %s %e %s.\n',Dest.mResults.locVar.*uF_LV,char(177),Dest.mResults.locVarSigma.*uF_LV,uL_LV)
            else
                [uF_LV,uL_LV] = parent.getUnitFactor('pixelsize');
                fprintf('The mean localization error is %e %s %e %s.\n',sqrt(Dest.mResults.locVar.*uF_LV),char(177),sqrt(Dest.mResults.locVarSigma.*uF_LV),uL_LV)
            end
            fprintf('The mean velocity is %e %s %e %s.\n',Dest.mResults.velocity.*uF_v,char(177),Dest.mResults.velocitySigma.*uF_v,uL_v)
            
        end
        
    end
end