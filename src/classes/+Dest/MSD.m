classdef MSD < Dest.DEstimator
    % MSD Mean squared displacement diffusion class
    %
    % Linear fit of the mean squared displacement. Incorporates the
    % camera-based diffusion model of Berglund et al.:
    % PHYSICAL REVIEW E 85, 061916 (2012) Optimal diffusion 
    % coefficient estimation in single-particle tracking
    
    % ---------------------------------------------------------------------
    % Copyright (C) 2019 J.J. Erik Maris
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
        estName = 'MSD'; % name of the diffusion estimator
        supportedModels = {'MSD_normal','MSD_confined','MSD_directed'}; % model available for this diffusion estimator
    end
    
    % parameters for GUI
    properties (Constant, Hidden)
        inputdlgPrompt = {'Clipping factor (fraction < 1; N points > 1):','Minimum points taken for fit:','Compute population MSD fit only:'}; % inputdlg prompt
        inputdlgVariables = {'clip_factor','minFitRange','meanOnlyFlag'}; % names of variables corresponding to prompt
        plotdlgNames = {'plotMSDfit','plotMeanMSDfit'}; % plot function names
        plotdlgDescription = {'Plot MSD fit of track','Plot MSD fit of all tracks'}; % description for gui
        plotdlgType = {'track','population'}; % track, population, bootstrap
    end
    
    % output values
    properties
        r2fit % vector // tracks R-squared goodness of fit
        fit_obj % vector // tracks fit object
        gof_obj % vector // tracks goodness of fit object
        coeffvals % vector // raw data used for plotting
        Nfit % vector // number of data points contributing to fit
        results % structure // coeffvals 
        mR2fit % vector // population R-squared goodness of fit
        mFit_obj % vector // population fit object
        mGof_obj % vector // population goodness of fit object
        mCoeffvals % vector // raw data used for plotting mean
        mNfit % vector // number of data points contributing to fit
        mResults % structure // coeffvals 
    end
    
    % input parameters
    properties
        clip_factor = 0.25; % scalar // specifies the part of the MSD curve which is fitted; clip_factor < 1 is interpreted as a fraction >= as # data points
        minFitRange = 3; % scalar // minimum number of data points for fit (overrules the clip_factor)
        DModel = Dest.MSD_normal; % object // model properties
        meanOnlyFlag = false; % logical // compute only mean fit?
    end
    
    methods
        % --- constructor
        function obj = MSD(varargin)
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
        
        function obj = set.clip_factor(obj,value)
            validateattributes(value,{'numeric'},{'scalar','nonnegative'},mfilename,'clip_factor')
            obj.clip_factor = value;
        end
        
        function obj = set.minFitRange(obj,value)
            validateattributes(value,{'numeric'},{'scalar','nonnegative','integer'},mfilename,'minFitRange')
            obj.minFitRange = value;
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
            % Florian Meirer 2016: added clipping factor and minimum fit
            % points
            %
            % Based on msdanalyzer https://tinevez.github.io/msdanalyzer/
            % Copyright (C) 2013 - 2014 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            % --- input checking
            if numel(parent) > 1
                error('MSD:compute:multipleObjs',...
                    '%d objects parsed, can only process one.',numel(parent))
            end
            
            if ~parent.MSD_valid
                parent = parent.computeMSD;
            end
            
            % only compute all tracks if menaOnlyFlag is false
            if ~obj.meanOnlyFlag
                
                nTracks = parent.nTracks;

                if obj.clip_factor < 1
                    fprintf('Fitting %d curves of MSD = f(x), taking only the first %d%% of each curve...\nf(x) = %s\n',...
                        nTracks, ceil(100 * obj.clip_factor), obj.DModel.fitfunction )
                else
                    fprintf('Fitting %d curves of MSD = f(x), taking only the first %d points of each curve...\nf(x) = %s\n',...
                        nTracks, round(obj.clip_factor), obj.DModel.fitfunction )
                end

                minMSDpoints = min(parent.nTrackPoints)-1;

                if obj.minFitRange > minMSDpoints
                    error('MSD:compute:largeMinFitRange',...
                        'The minimum points taken for fit (%d) is larger than the number of points in the MSD curve of the shortest track (%d). The latter is the number of track points minus one. Please decrease the minimum points taken for fit.',obj.minFitRange,minMSDpoints)
                end

                % --- preallocation

                r2fit_ = NaN(nTracks, 1);
                Nfit_ = NaN(nTracks, 1);
                fit_obj_ = cell(nTracks,1); % addition F. Meirer 2016-01-23
                gof_obj_ = cell(nTracks,1); % addition F. Meirer 2016-01-23
                coeffvals_ = zeros(nTracks,obj.DModel.nCoeff); % addition E. Maris 2019-05-11

                % --- core
                dispVal = round(nTracks/10);
                excluded_curves = 0;
                n_minFitRange_used = 0;

                for ii = 1:nTracks

                    if mod(ii,dispVal) == 0; fprintf('%d%% ',ii/dispVal*10); end % display

                    msd_spot = parent.MSD{ii};

                    t = msd_spot(:,1);
                    y = msd_spot(:,2);
                    w = msd_spot(:,4); % Number of points

                    % as the code was it used the percentage of the longest msd curve to calculate the
                    % number of points to fit; this of course makes no sense so this was added:
                    % correction F.Meirer 13.01.2016
                    % correct to the actual length of the track first:
                    nonnan = ~isnan(y);
                    t = t(nonnan);
                    y = y(nonnan);
                    w = w(nonnan);

                    % Clip data, never take the first one dt = 0
                    if obj.clip_factor < 1
                        t_limit = 2 : round(numel(t) * obj.clip_factor);
                    else
                        t_limit = 2 : min(1+round(obj.clip_factor), numel(t));
                    end
                    % addition F.Meirer to have both a percentage AND a minimum number of
                    % points to fit; 2016-01-25
                    % check whether percentage is too short to provide at least minFitRange
                    % range:
                    if numel(t_limit) < obj.minFitRange
                        % and correct if necessary:
                        t_limit = 2:obj.minFitRange+1;
                        n_minFitRange_used = n_minFitRange_used+1;
                    end

                    Nfit_(ii) = numel(t_limit); % Erik Maris 30-04-2019

                    t = t(t_limit);
                    y = y(t_limit);
                    w = w(t_limit);

                    if numel(y) < 2
                        excluded_curves = excluded_curves + 1;
                        continue % need at least 2 points to fit => skip this  msd curve and move to next
                    end

                    [fo, gof] = fit(t, y, obj.DModel.fitfunction, 'Weights', w,...
                        'StartPoint', obj.DModel.getStartPoints(t,y),...
                        'Lower', obj.DModel.lowerBounds, 'Upper', obj.DModel.upperBounds);
    %                 [fo, gof] = fit(t, y, 'poly1', 'Weights', w); % linear relation for free diffusion

                    % adjusted r-square is a bit tricky to use here...
                    %r2fit(i_spot) = gof.adjrsquare;
                    % use real r-square instead: correction F.Meirer 24.10.2015
                    r2fit_(ii) = gof.rsquare;
                    % also store fit results: addition F. Meirer 2016-01-23
                    fit_obj_{ii} = fo;  
                    gof_obj_{ii} = gof;
                    coeffvals_(ii,:) = coeffvalues(fo);

                end
                fprintf('\n')
                fprintf('Excluded %i MSD curves from fitting based on user-defined fitting range.\n',excluded_curves);
                fprintf('For %i MSD curves the minimum fitting range was used instead of the provided percentage.\n',n_minFitRange_used);


                % => corrected by F. Meirer 2016-01-23 / E.Maris 2019-05-11:
                obj.r2fit = r2fit_;
                obj.fit_obj = fit_obj_;
                obj.gof_obj = gof_obj_;
                obj.coeffvals = coeffvals_;
                obj.Nfit = Nfit_;

                % computes the values we're interested in, such as diffusion
                % constant and localization variance
                obj = obj.DModel.computeResults(obj,parent); % computes per track 
            end
            
            
            % --- mean MSD analysis
            
            fprintf('Fitting mean curve.\n')

            msmsd = getMeanMSD(parent);

            t = msmsd(:,1);
            y = msmsd(:,2);
            w = 1./msmsd(:,3).^2; % 1/var https://forum.image.sc/t/question-regarding-ensemble-msd-fit-in-msdanalyzer/44066/3

            % Clip data, never take the first one dt = 0
            if obj.clip_factor < 1
                t_limit = 2 : round(numel(t) * obj.clip_factor);
            else
                t_limit = 2 : min(1+round(obj.clip_factor), numel(t));
            end
            Nfit_ = numel(t_limit); % Erik Maris 16-9-2019

            % addition F.Meirer to have both a percentage AND a minimum number of
            % points to fit; 2016-01-25
            % check whether percentage is too short to provide at least minFitRange
            % range:
            if numel(t_limit) < obj.minFitRange
                % and correct if necessary:
                t_limit = 2:obj.minFitRange+1;
            end

            t = t(t_limit);
            y = y(t_limit);
            w = w(t_limit);
            
            [fo, gof] = fit(t, y, obj.DModel.fitfunction, 'Weights', w, 'StartPoint', obj.DModel.getStartPoints(t,y)); 
            
            % following the same protocol as above
            obj.mR2fit = gof.rsquare;
            obj.mFit_obj = fo;  
            obj.mGof_obj = gof;
            obj.mCoeffvals = coeffvalues(fo);
            obj.mNfit = Nfit_;
            
            % computes the values we're interested in, such as diffusion
            % constant and localization variance
            obj = obj.DModel.computeMeanResults(obj,parent); % computes per mean fit

            fprintf('Done\n')
            
        end
        
        function varargout = plotMSDfit(obj, parent, ha, index, errorbar)
            % PLOTMSDFIT Plots the fit of the MSD for one track
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % plotMSDfit(obj, parent, ha, index, errorbar)
            %   ha: handle to axes. If empty the current axis is taken
            %   index: index of the track to be displayed
            %   error: if true, an errorbar with the standard devation is
            %   shown, otherwise no errorbar is displayed
            %
            % h = plotMSDfit(...) returns the handle to the line plotted.
            %
            % [h, ha] = plotMSDfit(...) also returns the handle of the
            % axes in which the curve was plotted.
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % Based on msdanalyzer https://tinevez.github.io/msdanalyzer/
            % Copyright (C) 2013 - 2014 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University

            % --- input checking
            if numel(parent) > 1
                error('MSD:plotMSDfit:multipleObjs',...
                    '%d objects parsed, can only process one.',numel(parent))
            end
            
            if nargin < 3 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 4
                index = 1;
            end
            
            if ~isnumeric(index) && ~isscalar(index)
                error('plotMSDfit:badInput', ...
                    'Invalid index.')
            end
            
            if nargin < 5
                errorbar = false;
            end

            % --- core
            msd_spot = parent.MSD{index};
            if isempty( msd_spot ) % do not plot empty spots in MSD
                return
            end
            
            [xUnitFactor,xUnitLabel] = parent.getUnitFactor('dt');
            [yUnitFactor,yUnitLabel] = parent.getUnitFactor('pixelsize^2');

            t = msd_spot(:,1).*xUnitFactor;
            m = msd_spot(:,2).*yUnitFactor;

            X = linspace(0,max(msd_spot(:,1)),100);
            Y = feval(obj.fit_obj{index},X).*yUnitFactor;
            X = X.*xUnitFactor;
            % MSDfun(X,obj.coeffvals(index,1),obj.coeffvals(index,2))
            
            % --- plotting
            hps = NaN(2,1);
            if errorbar
                s = msd_spot(:,3).*yUnitFactor; % standard error of the mean
                H = tracks.errorShade(ha,t,m,s, [0 0 0], false); % non-transparent; black
                H.mainLine.LineWidth = 1.5;
                H.mainLine.Marker = 'o';
                hps(1) = H.mainLine;
            else
                hps(1) = plot(ha,t,m,'-o','LineWidth', parent.lineWidth*(2/3));
            end
            hold on;
            hps(2) = plot(ha,X,Y,'r','LineWidth', parent.lineWidth);
            title(sprintf('Fit to MSD of track %i',index));
            xlabel(['Delay time (' xUnitLabel ')'])
            ylabel(['MSD (' yUnitLabel ')'])
            hold off
            
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
        function varargout = plotMeanMSDfit(obj, parent, ha, errorbar, indices)
            % plotMeanMSDfit Plots the fit of the mean MSD for one track
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % plotMeanMSDfit(obj, parent, ha, errorbar, indices)
            %   ha: handle to axes. If empty the current axis is taken
            %   errorbar: if true, an errorbar with the standard devation is
            %   shown, otherwise no errorbar is displayed
            %   indices: indices of the tracks taken into the mean MSD
            %
            % h = plotMeanMSDfit(...) returns the handle to the line plotted.
            %
            % [h, ha] = plotMeanMSDfit(...) also returns the handle of the
            % axes in which the curve was plotted.
            % 
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % Based on msdanalyzer https://tinevez.github.io/msdanalyzer/
            % Copyright (C) 2013 - 2014 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University

            % --- input checking
            if numel(parent) > 1
                error('MSD:plotMeanMSDfit:multipleObjs',...
                    '%d objects parsed, can only process one.',numel(parent))
            end
            
            if nargin < 5
                indices = [];
            end

            if nargin < 4
                errorbar = false;
                if nargin < 3
                    ha = gca;
                end
            end

            % --- core

            [xUnitFactor,xUnitLabel] = parent.getUnitFactor('dt');
            [yUnitFactor,yUnitLabel] = parent.getUnitFactor('pixelsize^2');
            
            % --- plotting
            [hps,ha] = parent.plotMeanMSD(ha, errorbar, indices);
            xlims = ha.XLim;
            
%             MSDfun = @(X,a,b) X.*a.*(yUnitFactor./xUnitFactor) +b.*yUnitFactor;
            X = linspace(0,xlims(2)/xUnitFactor,100); % do feval in reduced units
            Y = feval(obj.mFit_obj,X).*yUnitFactor;
            X = X.*xUnitFactor;
            
            hold on
            ht = plot(ha,X,Y,'r','LineWidth',parent.lineWidth);
            xlabel(['Delay time (' xUnitLabel ')'])
            ylabel(['MSD (' yUnitLabel ')'])
            hold off
            
%             fprintf('Estimation of the diffusion coefficient from linear fit of the mean MSD curve:\n')
%             fprintf('D = %.3g %s/%s and localization variance = %.3g %s\n', parent.getD('mean'), yUnitLabel,xUnitLabel,parent.getLocVar('mean'),yUnitLabel);

            if nargout > 0
                varargout{1} = [hps(:); ht];
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
    end
end

