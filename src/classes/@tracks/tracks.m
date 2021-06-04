classdef tracks < classHandling & unitProps & plotProps
    % TRACKS class to deal with coordinates/time organised data. 
    %
    % Contains methods to segment tracks based on their properties
    %
    % obj = tracks()
    %   Creates an empty instance of this class.
    %
    % obj = tracks(time,coords,pixel_m,time_s)
    %   time: cell array column vector with a vector with the
    %   numeric time tags of the track
    %   coords: cell array column vector with a matrix with the
    %   numeric coordinates [N,#dimensions]
    %   pixel_m: the length unit of the coordinates as string (e.g.
    %   'nm') or as value in meters (e.g. 64e-9 for 64 nm pixels)
    %   time_s: the time unit of the time tags as string (e.g.
    %   'ms') or as value in seconds (e.g. 0.03 for 30 ms)
    %   
    % SEE ALSO: importDoM and ImportLocalizer
    
    
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
    
    % properties values
    properties %(SetAccess = protected)
        features % table with track features
        doFeatures = {...
            'nTrackPoints',...
            'nTrackPointsAll',...
            'firstFrame',...
            'onFraction',...
            'blinkingRate',...
            'OFFON',...
            'ONOFF',...
            'MinBoundCircleRadius',... % radius of the minimum enclosing circle
            'MinBoundCircleCenter',... % coordinates of the center of mass
            'CenterOfMass',... % coordinates of the center of mass
            'MBCCminusCoM',... % minimum enclosing circle minus center of mass distance
            'entropy',... %  Shannon's entropy of the track path
            'trackLength',... % length of track
            'tortuosity',... % length of track divided by the Euclidean end to end distance
            'EVec1',... % vector first principal component
            'EVec2',... % vector second principal component
            'CVE1',... % weight first principal component
            'CVE2',... % weight second principal component
            'wEV1',... % normalized weight first principal component
            'EV1angle',... % angle between first and second principal component [degree]
            'voronoiSA'... % area of the Voronoi surface
            };
        % cell with property names that have to be computed
        resolution = 10; % resolution of the entropy
    end
    
    % output properties
    properties
        segTree = segTree(); % segmentation tree
        name; % name of the tracks
    end
    
    % input values
    properties(SetAccess = protected)
        filepath % file path of DoM results table
        fitProps % fit properties
        coordsZero % zero of the cartesian coordinate system
        coords % track coordinates
        time % frame numbers of coordinates
        dt = 1; % integration time (before unit conversion!, e.g. in frames) // dt = 0 means that the tracks have not been recorded as a time lapse
        allTimes = []; % all times
%         locVar =[]; % localization error
    end
    
    % drift correction
    properties(SetAccess = protected)
        drift % drift vectors [t coords/nDim] (computeDrift)
    end
    
    % MSD
    properties(SetAccess = protected)
        MSD % mean squared displacements (computeMSD)
    end
    
    % blinking
    properties(SetAccess = protected)
        intTrace % intensity trace
        onTrace % shorthand notation of the on-off trace
        onTraceBinEdges % onTrace bin edges
        ONOFF
        OFFON
    end

    % diffusion constant properties
    properties
        Dest % object for diffusion constant estimation (computeDiffusion)
        dte = 0; % exposure time (before unit conversion!, e.g. in frames) // dte = 0 no motion blur
    end
    
    properties (Dependent)
        nDim % number of dimensions
        nFrames % number of frames (time step)
        nTracks % number of tracks
        nTrackPoints % number of points per track it is on
%         timeLapseMode; % data has been recorded with equidistant time intervals
    end
    
    properties (Hidden = true)
        features_valid = false;
        MSD_valid = false;
        Dest_valid = false;
        onTrace_valid = false; % If false, onTrace have to be recomputed
        distOnOff_valid = false; % If false, distOnOff needs to be recomputed
    end
    
    properties (Constant)
        TOLERANCE = 12;  % Tolerance for binning delays together. Two delays will be binned together if they differ in absolute value by less than 10^-TOLERANCE.
    end
    
    % constructor function
    methods
        function obj = tracks(time,coords,pixel_m,time_s,varargin)
            % TRACKS Constructs an instance of this class
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % obj = tracks()
            %   Creates an empty instance of this class.
            %
            % obj = tracks(time,coords,pixel_m,time_s)
            %   time: cell array column vector with a vector with the
            %   numeric time tags of the track
            %   coords: cell array column vector with a matrix with the
            %   numeric coordinates [N,#dimensions]
            %   pixel_m: the length unit of the coordinates as string (e.g.
            %   'nm') or as value in meters (e.g. 64e-9 for 64 nm pixels)
            %   time_s: the time unit of the time tags as string (e.g.
            %   'ms') or as value in seconds (e.g. 0.03 for 30 ms)
            %   
            % SEE ALSO: importDoM and ImportLocalizer
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % 
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            
            % construct superclasses
            
            obj = obj@unitProps();
            
            if nargin <= 2 % make sure that pixelsize and dt always do exist
                u = symunit;
                obj = obj.setUnit('pixelsize',1*u.m,'dt',1*u.s);
%                 obj.unitSystem = [obj.unitSystem newUnit('pixelsize',1*u.m) newUnit('dt',1*u.s)];
            end
            
            if nargin > 0
                obj.coords = coords;
                obj.time = time;
                u = symunit;
            end
            
            if nargin > 2
                if isnumeric(pixel_m)
                    pixel_m = pixel_m*u.m; % only symunits and chars accepted
                end
                if isnumeric(time_s)
                    time_s = time_s*u.s;
                end
                obj = obj.setUnit('pixelsize',pixel_m,'dt',time_s);
            end
            
            if nargin > 4
                nIn = numel(varargin);
                if mod(nIn,2) ~= 0
                    error('Variables can only be parsed as Name-Value pairs')
                end
                for ii = 1:nIn/2
                    try
                        obj.(varargin{2*ii-1}) = varargin{2*ii};
                    catch
                        error(['The property ''' varargin{2*ii-1} ''' is not recognised.'])
                    end
                end
            end
            
        end
        
        % dependent properties
        
        function ndim = get.nDim(obj)
            if isempty(obj.coords)
                ndim = 0;
            else
                ndim = size(obj.coords{1},2);
            end
        end
        
        function nFrames = get.nFrames(obj)
            nFrames = max(vertcat(obj.time{:}));
        end
        
        function nTracks = get.nTracks(obj)
            nTracks = size(obj.coords,1);
        end
        
        function nTrackPoints = get.nTrackPoints(obj)
            nTrackPoints = zeros(obj.nTracks,1);
            for ii = 1:obj.nTracks
                nTrackPoints(ii) = size(obj.coords{ii},1);
            end
        end
        
        % set properties
        
        function obj = set.dt(obj,value)
            validateattributes(value,{'numeric'},{'nonnegative'})
            obj.dt = value;
        end
        
        function obj = set.dte(obj,value)
            validateattributes(value,{'numeric'},{'nonnegative'})
            obj.dte = value;
        end
        
        function obj = set.Dest(obj,value)
%             validateattributes(value,{''},{'nonempty'})
            obj.Dest = value;
        end
        
        function obj = set.resolution(obj,value)
            validateattributes(value,{'numeric'},{'nonnegative','integer'})
            obj.resolution = value;
        end
        
        % standard functions
        
        function R = getR(obj)
            % computeR Computes the motion blur coefficient
            % 
            % Computes the motion blur coefficient from the frame time
            % (obj.dt) and exposure time (obj.dte) based on Eq. S9 in 
            % M. Linden and J. Elf. Biophys. J. 115.2 (2018): 276-282
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % obj = computeR(obj)
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % 
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            R = nan(numel(obj),1);
            for ii = 1:numel(obj)
                R(ii) = (1/6)*obj(ii).dte/obj(ii).dt; % no unit conversion required
            end
        end
         
     end
    
    methods (Static,Hidden)
        
        function Y = roundn(X,N)
            % ROUNDN Round towards nearest number with Nth decimals.
            % 
            % ROUNDN(X,N) rounds the elements of X to the nearest numbers 
            % with the precision given by N.
            %
            % Examples:     roundn(8.73,0) = 9
            %               roundn(8.73,1) = 8.7
            %               roundn(8.73,2) = 8.73
            %               roundn(8.73,-1) = 10
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % Reproduced from msdanalyzer https://tinevez.github.io/msdanalyzer/
            % (C) 2007 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2021 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University

            Y = 10^(-N) .* round(X.*10^(N));
        end
        
        % computeCenter
        [center,radius] = minboundcircle(x,y,hullflag)
        
        %thresholdTracks
        keepBoolean = thresholdTracksBoolean(tracks,thresholdField,minThreshold,maxThreshold)
        
        %plotVoronoiDiag
        hPs = getHullPoints(XY,n)
        
        function wm = weightedmean(x, w)
            % weightedmean Computes the weighted mean
            % 
            % wm = weightedmean(x, w) computes the weighted mean wm of the
            % elements x with a weight w.
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % Reproduced from msdanalyzer https://tinevez.github.io/msdanalyzer/
            % (C) 2007 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            wm = sum( x .* w) ./ sum(w);
        end
        
        function sewm = standarderrorweightedmean(x, w)
            % standarderrorweightedmean Computes the standard error weighted mean
            % 
            % sewm = standarderrorweightedmean(x, w) computes the standard
            % error weighted mean sewm of the elements x with a weight w.
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % Reproduced from msdanalyzer https://tinevez.github.io/msdanalyzer/
            % (C) 2007 Jean-Yves Tinevez
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2021 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            n = numel(w);
            wbar = mean(w);
            xbar = sum( x .* w) ./ sum(w);
            sewm = n /((n-1) * sum(w)^2) * (sum( (w.*x - wbar*xbar).^2) ...
                - 2 * xbar * sum( (w-wbar).*(w.*x - wbar*xbar)) ...
                + xbar^2 * sum((w-wbar).^2));
        end
        
        H = errorShade(ha, x, y, errBar, col, transparent)
        
    end
    
    
    methods (Static)
        
        obj = importDoM(csvpath,falsePositivesFlag,pixel_m,time_s)
        obj = importLocalizer(csvpath,pixel_m,time_s)
        obj = importTrackpy(csvpath,pixel_m,time_s)
        obj = importCOMSOL(csvpath,nDim,pixel_m,time_s)
        names = getDestNames(folderpath)
        
    end
    
    methods (Static, Hidden)
        
        MSD = static_getMSD(t,X);
        msmsd = static_getMeanMSD(MSD);
        dR2 = static_getdR2(t,X,delayTimes);
        
    end
end

