function obj = importTrackpy(csvpath,pixel_m,time_s)
% importTrackpy Imports the TrackPy export table for tracks
% 
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj = importTrackpy(csvpath,pixel_m,time_s)
%   csvpath: path to the input file (.csv), empty gives dialog
%   pixel_m: the length unit of the coordinates as string (e.g. 'nm') or as
%   value in meters (e.g. 64e-9 for 64 nm pixels)
%   time_s: the time unit of the time tags as string (e.g. 'ms') or as 
%   value in seconds (e.g. 0.03 for 30 ms)

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
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

% --- input dialog
if isempty(csvpath)
    [file,path] = uigetfile('*.csv', 'Multiselect', 'on');
    if isequal(file,0)
        % do nothing
    else
       csvpath = fullfile(path,file);
    end
end

if iscell(csvpath) % if is cell, the csv are batch processed within 'batchProcessTracks'
    obj = batchProcessTracks(csvpath,pixel_m,time_s);
else
    csvout = csvread(csvpath,1,0);
    obj = importAsTracks(csvout,pixel_m,time_s,csvpath);
end

end

function obj = batchProcessTracks(csvpath,pixel_m,time_s)

numcsv = length(csvpath);

% preallocation
obj = tracks.empty(numcsv,0);

for ii = 1:numcsv
    try
        obj(ii,1) = tracks.importTrackpy(csvpath{ii},pixel_m,time_s);
    catch
        error(['Parsed .csv (' csvpath{ii} ') are not all of the type ''tracks'' in ''' mfilename '''.'])
    end
end

end

function obj = importAsTracks(csvout,pixel_m,time_s,csvpath)

nTracks = max(csvout(:,10)); % maximum written trackID

coords = cell(nTracks,1);
time = cell(nTracks,1);
mass = cell(nTracks,1);
size = cell(nTracks,1);
ecc = cell(nTracks,1);
signal = cell(nTracks,1);
raw_mass = cell(nTracks,1);
ep = cell(nTracks,1);

for ii = 1:nTracks
    take = csvout(:,10) == ii;
    coords{ii} = [csvout(take,2) csvout(take,1)];
    time{ii} = csvout(take,9);
    mass{ii} = csvout(take,3);
    size{ii} = csvout(take,4);
    ecc{ii} = csvout(take,5);
    signal{ii} = csvout(take,6);
    raw_mass{ii} = csvout(take,7);
    ep{ii} = csvout(take,8);
end

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
T = table(mass,size,ecc,signal,raw_mass,ep);
obj.fitProps = [obj.fitProps T];
datasetID = strings(nTracks,1);
[~,filename,~] = fileparts(csvpath);
datasetID(:) = filename;
obj.fitProps.datasetID = datasetID; % overwrite datasetID with filepath
obj.filepath = csvpath;

% obj.fitProps.mass = mass;
% obj.fitProps.size = size;
% obj.fitProps.ecc = ecc;
% obj.fitProps.signal = signal;
% obj.fitProps.raw_mass = raw_mass;
% obj.fitProps.ep = ep;

% obj.coordsZero = 1;
obj.dt = 1; % frame
% obj.dte = 1;

end