classdef MSD_normal < Dest.DModel
    %NORMAL Model for normal diffusion
    
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
        name = 'MSD_normal';
    end
    
    % MSD
    properties (Constant)
        fitfunction = 'p1*x + p2';
        lowerBounds = [];
        upperBounds = [];
        nCoeff = 2;
    end
    
    methods (Static)
        function Dest = computeResults(Dest,parent)
            % computes the derived values and stores in estimator object
            
            D = Dest.coeffvals(:,1) ./ 2 ./ parent.nDim; % Eq. B1
            locVar = Dest.coeffvals(:,2) ./ 2 ./ parent.nDim + 2 .* parent.getR .* D .* parent.dt; % Eq. B1
            Dest.results = table(D,locVar);
            Dest.results.Properties.VariableDescriptions = {'Diffusion constant','Localization variance'};
            Dest.results.Properties.VariableUnits = {'pixelsize.^2/dt','pixelsize.^2'};
            if parent.dt ~= 0
                diffusionSNR = real(sqrt(D.*parent.dt./locVar)); % definition Vestergaard
                diffusionSNR(D < 0 | locVar < 0) = nan;
                Dest.results.diffusionSNR = diffusionSNR;
                Dest.results.Properties.VariableDescriptions{end} = 'Diffusion SNR';
                Dest.results.Properties.VariableUnits{end} = '';
            end
        end
        function Dest = computeMeanResults(Dest,parent)
            % computes the derived values and stores in estimator object
            
            Dest.mResults = struct('D',[],'LocVar',[],'Dsigma',[],'LocVarSigma',[]);

            if parent.dt ~= 0
                Dest.mResults.D = Dest.mCoeffvals(1) ./ 2 ./ parent.nDim; % Eq. B1
                Dest.mResults.locVar = Dest.mCoeffvals(2) ./ 2 ./ parent.nDim + 2 .* parent.getR .*Dest.mResults.D .* parent.dt; % Eq. B1
            end

            % compute errors
            try % if enough data points for error estimation are present
                CI = confint(Dest.mFit_obj,1-2*tcdf(-1,Dest.mGof_obj.dfe)); % https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule, level = 2*tcdf(-1,Dest.mGof_obj.dfe); is same as level: https://nl.mathworks.com/matlabcentral/answers/34234-how-to-obtain-std-of-coefficients-from-curve-fitting
                CI = (CI(2,:) - CI(1,:))./2;
                Dest.mResults.Dsigma = CI(1) ./ 2 ./ parent.nDim; % Eq. B1
                Dest.mResults.locVarSigma = CI(2) ./ 2 ./ parent.nDim + 2 .* parent.getR .*Dest.mResults.D .* parent.dt; % Eq. B1
            catch e
                warning(e.message) % retrow error as warning
            end

        end

        function StartPoints = getStartPoints(x,y)
            % get start points for fit y=p1*x+p2
            
            StartPoints = [mean(y)/mean(x) y(1)]; % function computes start points for poly1 automatically
            
        end
        
        function report(Dest,parent)
            
            [uF_D,uL_D] = parent.getUnitFactor('pixelsize.^2./dt');
            [uF_LV,uL_LV] = parent.getUnitFactor('pixelsize.^2');
            fprintf('Mean squared displacement analysis for normal diffusion.\n')
            fprintf('The mean diffusion constant is %e %s %e %s.\n',Dest.mResults.D.*uF_D,char(177),Dest.mResults.Dsigma.*uF_D,uL_D)
            if Dest.mResults.locVar < 0
                fprintf('The mean localization error is < 0, giving the localization variance: %e %s %e %s.\n',Dest.mResults.locVar.*uF_LV,char(177),Dest.mResults.locVarSigma.*uF_LV,uL_LV)
            else
                fprintf('The mean localization error is %e %s %e %s.\n',sqrt(Dest.mResults.locVar.*uF_LV),char(177),(Dest.mResults.locVarSigma.*uF_LV)/(2*sqrt(Dest.mResults.locVar.*uF_LV)),uL_LV) % error propgagation https://chemistry.stackexchange.com/questions/91212/how-to-calculate-absolute-error-under-square-root-sign
            end
            
        end
        
    end

end