function obj = computeFeatures(obj,preserveFlag,featureNames)
%COMPUTEFEATURES Computes the track features and stores them in a table
%       'firstFrame': first frame in which track appears
%       'firstFrameCoords': coordinates of the first frame of a track
%       'offFrames': frame numbers when no particle is detected during the
%           track. Is a result of a blinking > 0
%       'nTrackPointsOff': number of tracks points at which no particle is
%           detected
%       'nTrackPointsAll': number of tracks points between first and last
%           detection of a track.
%       'MinBoundCircleRadius'/vector: radius of minimum enclosing circle
%           of the track (XY)
%       'MinBoundCircleCenter'/matrix: center coordinates of minimum 
%           enclosing circle of the track (XY)
%       'CenterOfMass'/matrix: center of mass of the track
%       'MBCCminusCoM'/matrix: % get difference between CoM and center of 
%           minimum bounding circle in percent. Idea: if this difference is
%           small the track is a 'ball of wool', if the difference is large
%           it is an 'uncoiled ball of wool': in percent of 
%           'MinBoundCircleRadius' (i.e. between 0 to 1).
%
% computeEntropy Computes the Shannon entropy of a track
% 
% Track coordinates are set on a binary grid with the size of the minimum
% bounding cirle's diameter and divided in resolution^2 equally sized
% blocks. The entropy is Shannon's entropy of a binary image.

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

nObj = numel(obj);

if nargin < 2
    preserveFlag = false;
end

if nargin < 3
    featureNamesFlag = false;
else
    featureNamesFlag = true;
end

wb = waitbar(0,'Computing Features','Name','Computing all track features');

for ii = 1:nObj
    waitstr = ['Computing Features (' num2str(ii) '/' num2str(nObj) ')'];
    waitbar(0.2/nObj+(ii-1)/nObj,wb,waitstr)
    
    %% Dimensionality check
    
    nDim = obj(ii).nDim;
    if nDim == 0 % is empty
        continue
    end
    
    
    %% create new table if does not exist
    if isempty(obj(ii).features)
        T = table();
        T.trackID = (1:obj(ii).nTracks)'; % only create track IDs
        T.Properties.VariableUnits = {''};
        T.Properties.VariableDescriptions = {'Track ID'};
    else
        if preserveFlag
            % copy existing table
            T = obj(ii).features;
            T = sortrows(T,'trackID'); % code assumes sorted columns
            N = numel(T.Properties.VariableNames);
            if isempty(T.Properties.VariableUnits) % make sure it exists
                T.Properties.VariableUnits = cell(1,N);
            end
            if isempty(T.Properties.VariableDescriptions) % make sure it exists
                T.Properties.VariableDescriptions = cell(1,N);
            end
        else
            % make new
            T = table();
            T.trackID = (1:obj(ii).nTracks)'; % only create track IDs
            T.Properties.VariableUnits = {''};
            T.Properties.VariableDescriptions = {'Track ID'};
        end
    end
    
    if featureNamesFlag
        do = featureNames;
    else
        do = obj(ii).doFeatures;
    end
    
    % clear to-be-computed features
    bool = ismember(T.Properties.VariableNames,do);
    if any(bool)
        T = removevars(T,T.Properties.VariableNames{bool});
    end
    
    waitbar(0.3/nObj+(ii-1)/nObj,wb,waitstr)
    
    %% nTrackPoints
    
    thisFeatures = {'nTrackPoints'};
    
    if any(ismember(thisFeatures,do))
        if ismember('nTrackPoints',do); T.nTrackPoints = obj(ii).nTrackPoints;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Number of Points'; end
    end

    
    %% nTrackPointsAll, onFraction
    
    thisFeatures = {'nTrackPointsAll','onFraction'};
    
    if any(ismember(thisFeatures,do))
        nTrackPointsAll  = cellfun(@(x) x(end) - x(1) + 1,obj(ii).time);
        onFraction = obj(ii).nTrackPoints./nTrackPointsAll;
    
        if ismember('nTrackPointsAll',do); T.nTrackPointsAll = nTrackPointsAll;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Number of all points'; end
        if ismember('onFraction',do); T.onFraction = onFraction;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'ON-fraction'; end
    end
    
    %% firstFrame
    
    thisFeatures = {'firstFrame'};
    
    if any(ismember(thisFeatures,do))
        firstFrame = cellfun(@(x) x(1),obj(ii).time); % take first frame of appearance from frames
    
        if ismember('firstFrame',do); T.firstFrame = firstFrame;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'First frame'; end
    end
    
    %% Blinking
    
    thisFeatures = {'blinkingRate','OFFON','ONOFF'};
    
    if any(ismember(thisFeatures,do))
        if ~obj(ii).onTrace_valid
            obj(ii) = obj(ii).computeOnTrace; % is slow
        end
        
        % preallocation
        blinkingRate = zeros(obj(ii).nTracks,1);
        OFFON = cell(obj(ii).nTracks,1);
        ONOFF = cell(obj(ii).nTracks,1);

        % compute
        for jj = 1:obj(ii).nTracks
            nShort = numel(obj(ii).onTrace{jj});
            % weirdly coded - check
            if obj(ii).onTrace{jj}(1) == 0
                nShortZ = nShort - 1; % compensate for zero to indicate ON
            else
                nShortZ = nShort;
            end
            % blinking is ON-OFF-ON
            nTrackPointsAll  = cellfun(@(x) x(end) - x(1) + 1,obj(ii).time);
            if mod(nShort,2 == 0) % is even, ends on ON
                blinkingRate(jj) = floor(nShortZ-1/2)/nTrackPointsAll(jj);
            else % is odd, ends on OFF
                blinkingRate(jj) = floor(nShortZ/2)/nTrackPointsAll(jj);
            end

            % testing
            oT = obj(ii).onTrace{jj};
            coT = cumsum(oT);
            OFFON{jj} = coT(2:2:end); % the last time with OFF marks the transition
            ONOFF{jj} = coT(1:2:end);
            OFFON{jj} = OFFON{jj}(OFFON{jj} < obj(ii).nFrames);
            ONOFF{jj} = ONOFF{jj}(ONOFF{jj} > 0);

        end

        % write to table as table
        if ismember('blinkingRate',do); T.blinkingRate = blinkingRate;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Blinking rate'; end
        if ismember('OFFON',do); obj(ii).OFFON = OFFON; end
        if ismember('ONOFF',do); obj(ii).ONOFF = ONOFF; end
    end
    
    waitbar(0.4/nObj+(ii-1)/nObj,wb,waitstr)
    %% Center
    
    thisFeatures = {'MinBoundCircleRadius','MinBoundCircleCenter',...
        'CenterOfMass','MBCCminusCoM'};

    if any(ismember(thisFeatures,do))

        % preallocation

        MinBoundCircleRadius = zeros(obj(ii).nTracks,1);
        MinBoundCircleCenter = zeros(obj(ii).nTracks,2);
        CenterOfMass = zeros(obj(ii).nTracks,2);
        MBCCminusCoM = zeros(obj(ii).nTracks,1);

        % get the minimum bounding circle of the point cloud forming the
        % track (see: smallest-circle problem)
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. MinBoundCircle, CenterOfMass, and MBCCminusCoM cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. MinBoundCircle, CenterOfMass, and MBCCminusCoM are computed in the XY-plane.',nDim,ii)
        end

        for jj = 1:obj(ii).nTracks
            % --- minimum bounding circle 2D
            [C,R] = obj(ii).minboundcircle(obj(ii).coords{jj}(:,1),obj(ii).coords{jj}(:,2),1);
            MinBoundCircleRadius(jj) = R;
            MinBoundCircleCenter(jj,:) = C;
            % --- center of mass (CoM)
            CenterOfMass(jj,:) = mean(obj(ii).coords{jj}(:,1:2),1); % take mean over all rows
            % get difference between CoM and center of minimum bounding circle
            % in percent; idea: if this difference is small the track is a
            % 'ball of wool', if the difference is large it is an 'uncoiled ball of wool': 
            % in percent of MinBoundCircleRadius (i.e. between 0 to 1)
            MBCCminusCoM(jj) = pdist2(MinBoundCircleCenter(jj,1:2),CenterOfMass(jj,1:2))/R; % default value for pdist2 is Eucledian
        end

        % write to table as table
        if ismember('MinBoundCircleRadius',do); T.MinBoundCircleRadius = MinBoundCircleRadius;...
                T.Properties.VariableUnits{end} = 'pixelsize'; T.Properties.VariableDescriptions{end} = 'MBC radius'; end
        if ismember('MinBoundCircleCenter',do); T.MinBoundCircleCenter = MinBoundCircleCenter;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'MCB center'; end
        if ismember('CenterOfMass',do); T.CenterOfMass = CenterOfMass;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Center of mass'; end
        if ismember('MBCCminusCoM',do); T.MBCCminusCoM = MBCCminusCoM;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'MBCC minus CoM'; end

    end
    
    waitbar(0.6/nObj+(ii-1)/nObj,wb,waitstr)
    %% Entropy
        
    thisFeatures = {'entropy'};

    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Entropy cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Entropy is computed in the XY-plane.',nDim,ii)
        end
        
        if ~ismember('MinBoundCircleRadius',do) || ~ismember('MinBoundCircleCenter',do)
            objc = obj(ii).computeFeatures(false,{'MinBoundCircleRadius','MinBoundCircleCenter'});
            R = objc.features.MinBoundCircleRadius;
            CC = objc.features.MinBoundCircleCenter;
            clear objc
        else
            % use embounding circle radius as width and height
            R = T.MinBoundCircleRadius;
            CC = T.MinBoundCircleCenter;            
        end

        % preallocation
        entropy_ = zeros(obj(ii).nTracks,1);

        for jj = 1:obj(ii).nTracks
            % idea: make it a binary image and get the entropy of the track
            imgSize = (round(R(jj)*obj(ii).resolution)*2)+2;
            I = zeros(imgSize,imgSize);
            % now place the points of the track in this image:

            relCoordX = uint8(round((obj(ii).coords{jj}(:,1) - (CC(jj,1)-R(jj)))*obj(ii).resolution)+1);
            relCoordY = uint8(round((obj(ii).coords{jj}(:,2) - (CC(jj,2)-R(jj)))*obj(ii).resolution)+1);
            I(sub2ind(size(I),relCoordY,relCoordX)) = 1; % avoid for loop
            entropy_(jj) = entropy(uint8(I));
        end

        % write to table as table
        if ismember('entropy',do); T.entropy = entropy_;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Entropy'; end

    end
    
    waitbar(0.8/nObj+(ii-1)/nObj,wb,waitstr)
    %% displacement
    
    thisFeatures = {'trackLength','tortuosity'};
        
    if any(ismember(thisFeatures,do))
    
        % preallocation
    
        trackLength = NaN(obj(ii).nTracks,1);
        tortuosity = NaN(obj(ii).nTracks,1);
        
        % compute
        for jj = 1:obj(ii).nTracks
            trackLength(jj) = sum(sqrt(sum((obj(ii).coords{jj}(1:end-1,:)-obj(ii).coords{jj}(2:end,:)).^2,2)));
            tortuosity(jj) = trackLength(jj) / sqrt(sum((obj(ii).coords{jj}(1,:)-obj(ii).coords{jj}(end,:)).^2));
        end
        
        % write to table as table
        if ismember('trackLength',do); T.trackLength = trackLength;...
                T.Properties.VariableUnits{end} = 'pixelsize'; T.Properties.VariableDescriptions{end} = 'Length'; end
        if ismember('tortuosity',do); T.tortuosity = tortuosity;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Tortuosity'; end
    
    end
    
    waitbar(0.9/nObj+(ii-1)/nObj,wb,waitstr)
    %% PCA
    
    thisFeatures = {'EVec1','EVec2','CVE1','CVE2','EV1angle'};
        
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Elongation (angle) cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Elongation (angle) is computed in the XY-plane.',nDim,ii)
        end
    
        % preallocation
    
        EVec1 = NaN(obj(ii).nTracks,2);
        EVec2 = NaN(obj(ii).nTracks,2);
        CVE1 = NaN(obj(ii).nTracks,1);
        CVE2 = NaN(obj(ii).nTracks,1);
        EV1angle = NaN(obj(ii).nTracks,1);

        % get the principal component of each track (direction of largest variance)
        for jj = 1:obj(ii).nTracks
            covM = cov(obj(ii).coords{jj}(:,1:2));
            [V,D]=eig(covM);   

            % correct all EVecs to be in positive half space: 
            % => if the y-part of EVec1 is negative, flip it by 180deg
            if V(2,1)<0
                V(2,1)=-V(2,1);
                V(1,1)=-V(1,1);
            end
            if V(2,2)<0
                V(2,2)=-V(2,2);
                V(1,2)=-V(1,2);
            end
            CVE(1,1) = D(1,1)/(D(1,1)+D(2,2));
            CVE(2,1) = D(2,2)/(D(1,1)+D(2,2));
            if D(1,1)>D(2,2)
                EVec1(jj,:) = V(:,1);
                EVec2(jj,:) = V(:,2);
                CVE1(jj) = CVE(1,1);
                CVE2(jj) = CVE(2,1);
            else
                EVec1(jj,:) = V(:,2);
                EVec2(jj,:) = V(:,1);
                CVE1(jj) = CVE(2,1);
                CVE2(jj) = CVE(1,1);
            end
            % get angle in degrees for each EVec1 (per definition:
            % angle between basis set vector [1 0] and EVec1)
            EV1angle(jj) = acosd(dot(EVec1(jj,:),[1 0])); % both vectors are of length 1

        end

        % compute normalized CVE1, which ranges from 0 to 1. In the 2D case,
        % CVE1 can be 1 (all points are on a line) or 0.5 (spherical cloud of
        % points). This scales 0.5-1 to 0-1.
        wEV1 = CVE1.*2 - 1;
    
        % write to table as table
        if ismember('EVec1',do); T.EVec1 = EVec1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'EVec1'; end
        if ismember('EVec2',do); T.EVec2 = EVec2;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'EVec2'; end
        if ismember('CVE1',do); T.CVE1 = CVE1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'CVE1'; end
        if ismember('CVE2',do); T.CVE2 = CVE2;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'CVE2'; end
        if ismember('wEV1',do); T.wEV1 = wEV1;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Elongation'; end
        if ismember('EV1angle',do); T.EV1angle = EV1angle;...
                T.Properties.VariableUnits{end} = ''; T.Properties.VariableDescriptions{end} = 'Elongation angle'; end
    
    end
    
    waitbar(1/nObj+(ii-1)/nObj,wb,waitstr)
    %% Voronoi
    
    thisFeatures = {'voronoiSA'};
        
    if any(ismember(thisFeatures,do))
        
        if nDim < 2
            error('Number of dimensions is %i for population %i. Voronoi area cannot be computed.',nDim,ii)
        elseif nDim > 2
            warning('Number of dimensions is %i for population %i. Voronoi area is computed in the XY-plane.',nDim,ii)
        end
        
        if ~ismember('CenterOfMass',do)
            objc = obj(ii).computeFeatures(false,{'CenterOfMass'});
            COM = objc.features.CenterOfMass;
            clear objc
        else
            COM = T.CenterOfMass;
        end
        
        % preallocation
        voronoiSA = zeros(size(COM,1),1);

        % get the convex hull and extend it by 5 pixels; then add the hull points
        % to the data to cut of the voronoi diagram at the extended hull
        hPs = obj(ii).getHullPoints(COM,10);
        startOfhPs = size(COM,1)+1;
        COM = vertcat(COM,hPs);
        [v,c]=voronoin(COM);
        
        for jj = 1:length(c)
            if all(c{jj}~=1)  % If at least one of the indices is 1,
                              % then it is an open region and we can't
                              % patch that. These are introduced as hull points
                              % to create a convex hull in the voronoi diagram
                 if jj<startOfhPs % only patch real data points
                    voronoiSA(jj) = polyarea(v(c{jj},1),v(c{jj},2));
                 end
            end
        end
    
        % write to table as table
        if ismember('voronoiSA',do); T.voronoiSA = voronoiSA;...
                T.Properties.VariableUnits{end} = 'pixelsize.^2'; T.Properties.VariableDescriptions{end} = 'Voronoi area'; end
    
    end
    
    %% Save
    
    % store features
    obj(ii).features = T;
    % set validity
    obj(ii).features_valid = true;
    
end

close(wb)

end

%Template:
%
% thisFeatures = {};
%     
% if any(ismember(thisFeatures,do))
% 
%     % preallocation
% 
% 
% 
%     % write to table as table
%     if ismember('clear',do); T.nTrackPointsAll = nTrackPointsAll; end
% 
% end

