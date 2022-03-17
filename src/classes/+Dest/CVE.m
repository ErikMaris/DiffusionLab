classdef CVE < Dest.DEstimator
    % CVE estimator
    
    % ---------------------------------------------------------------------
    % Copyright (C) 2022 J.J. Erik Maris
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
    
    properties (Constant)
        estName = 'CVE'; % name of the diffusion estimator
        supportedModels = {}; % model available for this diffusion estimator
    end
    
    % parameters for GUI
    properties (Constant, Hidden)
        inputdlgPrompt = {'Compute population MSD fit only:'}; % inputdlg prompt
        inputdlgVariables = {'meanOnlyFlag'}; % names of variables corresponding to prompt
        plotdlgNames = {}; % plot function names
        plotdlgDescription = {}; % description for gui
        plotdlgType = {}; % track, population, bootstrap
    end
    
    % output values
    properties
        results % structure // coeffvals 
        mResults % structure // coeffvals 
    end
    
    % input parameters
    properties
        DModel = []; % object // model properties, dummy, is not implemented
        meanOnlyFlag = false; % logical // compute only mean fit?
    end
    
    methods
        % --- constructor
        function obj = CVE(varargin)
            %MSD Construct an instance of this class
            if mod(nargin,2) ~= 0
                error('Variables can only be parsed as Name-Value pairs')
            end
            for ii = 1:nargin/2
                try
                    obj.(varargin{2*ii-1}) = varargin{2*ii};
                catch
                    error(['The property ''' varargin{2*ii-1} ''' is not recognised.'])
                end
            end
        end
        
        function obj = set.DModel(obj,value)
            if isa(value,'Dest.DModel')
                validatestring(value.name,obj.supportedModels,mfilename,'DModel');
                obj.DModel = value;
            else
                validatestring(value,obj.supportedModels,mfilename,'DModel');
                obj.DModel = Dest.(value);
            end
        end
        
        function obj = compute(obj,parent)
            % COMPUTE Computes this diffusion estimator
            % 
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % obj = compute(obj,parent)
            % more info, see: tracks.computeDiffusion
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2022 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            % --- input checking
            if numel(parent) > 1
                error('MSD:compute:multipleObjs',...
                    '%d objects parsed, can only process one.',numel(parent))
            end
            
            % --- compute dR2, dR and dT
            
            dR2 = cellfun(@(x) sum((x(1:end-1,:) - x(2:end,:)).^2,2),parent.coords,'UniformOutput',false); % Squared displacement 
            dR2mp1 = cellfun(@(dR2) sqrt(dR2(1:end-1)).*sqrt(dR2(2:end)),dR2,'UniformOutput',false); % displacement
            dt = cellfun(@(x) x(2:end,:) - x(1:end-1,:),parent.time,'UniformOutput',false); % delta t
            R = parent.getR;
            
            % only compute all tracks if meanOnlyFlag is false
            if ~obj.meanOnlyFlag
                
                % --- compute per track
                
                D = cellfun(@(dR2,dR2mp1,dt) mean(dR2)/(2*parent.nDim*mean(dt)) + mean(dR2mp1)/(parent.nDim*mean(dt)),dR2,dR2mp1,dt);
                locVar = cellfun(@(dR2,dR2mp1) R*mean(dR2)/parent.nDim + (2*R-1)*mean(dR2mp1)/parent.nDim,dR2,dR2mp1);
                diffusionSNR = real(sqrt(D.*parent.dt./locVar)); % definition Vestergaard
                diffusionSNR(D < 0 | locVar < 0) = nan;
                
                obj.results = table(D,locVar,diffusionSNR);
                obj.results.Properties.VariableDescriptions = {'Diffusion constant','Localization variance','Diffusion SNR'};
                obj.results.Properties.VariableUnits = {'pixelsize.^2/dt','pixelsize.^2',''};
            end

            D = mean(vertcat(dR2{:}))/(2*parent.nDim*mean(vertcat(dt{:}))) + mean(vertcat(dR2mp1{:}))/(parent.nDim*mean(vertcat(dt{:})));
            locVar = R*mean(vertcat(dR2{:}))/parent.nDim + (2*R-1)*mean(vertcat(dR2mp1{:})/parent.nDim);
            diffusionSNR = real(sqrt(D.*parent.dt./locVar)); % definition Vestergaard
            diffusionSNR(D < 0 | locVar < 0) = nan;

            obj.mResults = table(D,locVar,diffusionSNR);
            obj.mResults.Properties.VariableDescriptions = {'Diffusion constant','Localization variance','Diffusion SNR'};
            obj.mResults.Properties.VariableUnits = {'pixelsize.^2/dt','pixelsize.^2',''};

            fprintf('Done\n')
            
        end
        
    end
end

