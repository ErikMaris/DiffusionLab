function [out,varargout] = merge(obj,baseObj,shuffleFlag)

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

if nargin < 2
    baseObj = 1;
end

if nargin < 3
    shuffleFlag = false;
end    

nObj = numel(obj);
if nObj < 2
%     warning('Only one object parsed.')
    return
end

objkeys = cell(nObj,1);
objvalues = cell(nObj,1);
for ii = 1:nObj
    objkeys{ii} = keys(obj(ii).dataUnits);
    objvalues{ii} = values(obj(ii).dataUnits);
end

if ~isequal(objkeys{:})
    error('tracks:merge:InconsistentArray', ...
        'Data unit keys are inconsitent.')
end

if ~isequal(objvalues{:})
    error('tracks:merge:InconsistentArray', ...
        'Data unit values are inconsitent.')
end

if unique([obj.dt]) > 1
    error('tracks:merge:InconsistentArray', ...
        'dt values are inconsitent.')
end

if unique([obj.dte]) > 1
    error('tracks:merge:InconsistentArray', ...
        'dte values are inconsitent.')
end

out = tracks;

out.coords = cat(1,obj.coords);
out.time = cat(1,obj.time);
out.features = cat(1,obj.features);
out.doFeatures = obj(baseObj).doFeatures;
out.resolution = obj(baseObj).resolution;
out.segTree = obj(baseObj).segTree;
out.name = obj(baseObj).name;
out.MSD = cat(1,obj.MSD);
out.drift = cat(1,obj.drift);
out.dataUnits = obj(baseObj).dataUnits;
out.unitSystem = obj(baseObj).unitSystem;
out.dt = obj(baseObj).dt;
out.dte = obj(baseObj).dte;

try
    out.fitProps = cat(1,obj.fitProps);
catch
    warning('Cannot merge ''fitProps''.') % if fitprops have different contents
end

% coords zero


y = [];
if shuffleFlag
    nTracksTot = sum(horzcat(obj.nTracks));
    y = randsample(nTracksTot,nTracksTot);
    out.coords = out.coords(y);
    out.time = out.time(y);
    if ~isempty(out.features); out.features = out.features(y,:); end
    if ~isempty(out.drift); out.drift = out.drift(y); end
    if ~isempty(out.MSD); out.MSD = out.MSD(y); end
end

varargout{1} = y;


end


% objProps = {...
%     'drift','MinBoundCircleRadius','MinBoundCircleCenter','CenterOfMass',...
%     'MBCCminusCoM','entropy','Evec1','Evec2','CVE1','CVE2','EV1angle','voronoiSA',...
%     'MSD','dtSD','trackLength','tortuosity','SNR','DTree','fitProps',...
%     'coords','time','allTimes','locVar',...
%     'nTrackPointsAll','firstFrame','onFraction',...
%     'blinkingRate','intTrace','onTrace','onTraceBinEdges','ONOFF','OFFON'};
