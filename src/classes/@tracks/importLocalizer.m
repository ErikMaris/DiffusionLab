function obj = importLocalizer(txtpath,pixel_m,time_s)
% importLocalizer Imports the Localizer output for tracks
% 
% Localizer > Manipulate Particle Tracks > Save tracks to text file...
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj = importLocalizer(txtpath,pixel_m,time_s)
%   txtpath: path to the input file (.txt), empty gives dialog
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
if isempty(txtpath)
    [file,path] = uigetfile('*.txt', 'Multiselect', 'on');
    if isequal(file,0)
        % do nothing
    else
       txtpath = fullfile(path,file);
    end
end

if iscell(txtpath) % if is cell, the csv are batch processed within 'batchProcessTracks'
    obj = batchProcessTracks(txtpath,pixel_m,time_s);
else
    a = importdata(txtpath); % inefficient way to obtain the number of columns
    ncol = size(a.data,2);
    switch ncol
        case 8
            % MLEwG
            obj = importAsTracks8(txtpath,pixel_m,time_s);
        case 10
            % Gaussian fixed width
            obj = importAsTracks10(txtpath,pixel_m,time_s);
        case 12
            % Gaussian
            obj = importAsTracks12(txtpath,pixel_m,time_s);
        case 16
            % Eliptical gaussian
            obj = importAsTracks16(txtpath,pixel_m,time_s);
        otherwise
            error(['Localization method not supported by ''' mfilename '''.'])
    end
end

end

function obj = batchProcessTracks(txtpath,pixel_m,time_s)

numcsv = length(txtpath);

% preallocation
obj = tracks.empty(numcsv,0);

for ii = 1:numcsv
    try
        obj(ii,1) = tracks.importLocalizer(txtpath{ii},pixel_m,time_s);
    catch
        error(['Parsed .csv (' txtpath{ii} ') are not all of the type ''tracks'' in ''' mfilename '''.'])
    end
end

end

function obj = importAsTracks8(txtpath,pixel_m,time_s)

% MATLAB does not seem to like the file format; let do a line-by-line approach
f = fopen(txtpath);
C = textscan(f,'%s %d %s','HeaderLines',1); % skip first header

nTracks = C{2}(1);

coords = cell(nTracks,1);
time = cell(nTracks,1);
coords_err = cell(nTracks,1);
sigma = cell(nTracks,1);
z0 = cell(nTracks,1);
intensity = cell(nTracks,1);

ii = 1;
txtline = fgetl(f);
while ischar(txtline) % end of file
    if strcmp(txtline(1:11),'First frame') % begin track
        trackdata = textscan(f,'%f %f %f %f %f %f %f %f'); 
        coords{ii} = [trackdata{4} trackdata{5}];
        time{ii} = trackdata{1};
        coords_err{ii} = trackdata{7};
        z0{ii} = trackdata{6};
        intensity{ii} = trackdata{2};
        sigma{ii} = sigma{3};
        ii = ii + 1;
    end
    txtline = fgetl(f);
end
fclose(f);

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
obj.filepath = txtpath;
obj.fitProps.coords_err = coords_err;
obj.fitProps.z0 = z0;
obj.fitProps.intensity = intensity;

obj.coordsZero = 1;
obj.dt = 1; % frame
% obj.dte = 1;

end

function obj = importAsTracks10(txtpath,pixel_m,time_s)

% MATLAB does not seem to like the file format; let do a line-by-line approach
f = fopen(txtpath);
C = textscan(f,'%s %d %s','HeaderLines',1); % skip first header

nTracks = C{2}(1);

coords = cell(nTracks,1);
time = cell(nTracks,1);
coords_err = cell(nTracks,1);
z0 = cell(nTracks,1);
z0_err = cell(nTracks,1);
intensity = cell(nTracks,1);
intensity_err = cell(nTracks,1);

ii = 1;
txtline = fgetl(f);
while ischar(txtline) % end of file
    if strcmp(txtline(1:11),'First frame') % begin track
        trackdata = textscan(f,'%f %f %f %f %f %f %f %f %f %f'); 
        coords{ii} = [trackdata{4} trackdata{5}];
        time{ii} = trackdata{1};
        coords_err{ii} = [trackdata{9} trackdata{10}];
        z0{ii} = trackdata{6};
        z0_err{ii} = trackdata{11};
        intensity{ii} = trackdata{2};
        intensity_err{ii} = trackdata{7};
        ii = ii + 1;
    end
    txtline = fgetl(f);
end
fclose(f);

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
obj.filepath = txtpath;
obj.fitProps.coords_err = coords_err;
obj.fitProps.z0 = z0;
obj.fitProps.z0_err = z0_err;
obj.fitProps.intensity = intensity;
obj.fitProps.intensity_err = intensity_err;

obj.coordsZero = 1;

end

function obj = importAsTracks12(txtpath,pixel_m,time_s)

% MATLAB does not seem to like the file format; let do a line-by-line approach
f = fopen(txtpath);
C = textscan(f,'%s %d %s','HeaderLines',1); % skip first header

nTracks = C{2}(1);

coords = cell(nTracks,1);
time = cell(nTracks,1);
coords_err = cell(nTracks,1);
z0 = cell(nTracks,1);
z0_err = cell(nTracks,1);
sigma = cell(nTracks,1);
sigma_err = cell(nTracks,1);
intensity = cell(nTracks,1);
intensity_err = cell(nTracks,1);

ii = 1;
txtline = fgetl(f);
while ischar(txtline) % end of file
    if strcmp(txtline(1:11),'First frame') % begin track
        trackdata = textscan(f,'%f %f %f %f %f %f %f %f %f %f %f %f'); 
        coords{ii} = [trackdata{4} trackdata{5}];
        time{ii} = trackdata{1};
        coords_err{ii} = [trackdata{9} trackdata{10}];
        z0{ii} = trackdata{6};
        z0_err{ii} = trackdata{11};
        sigma{ii} = trackdata{3};
        sigma_err{ii} = trackdata{8};
        intensity{ii} = trackdata{2};
        intensity_err{ii} = trackdata{7};
        ii = ii + 1;
    end
    txtline = fgetl(f);
end
fclose(f);

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
trackID = (1:nTracks)';
T = table(trackID,coords_err,z0,z0_err,sigma,sigma_err,...
    intensity,intensity_err);
obj.fitProps = T;
obj.filepath = txtpath;
% obj.fitProps.coords_err = coords_err;
% obj.fitProps.z0 = z0;
% obj.fitProps.z0_err = z0_err;
% obj.fitProps.sigma = sigma;
% obj.fitProps.sigma_err = sigma_err;
% obj.fitProps.intensity = intensity;
% obj.fitProps.intensity_err = intensity_err;

obj.coordsZero = 1;
obj.dt = 1; % frame
% obj.dte = 1;

end

function obj = importAsTracks16(txtpath,pixel_m,time_s)

% MATLAB does not seem to like the file format; let do a line-by-line approach
f = fopen(txtpath);
C = textscan(f,'%s %d %s','HeaderLines',1); % skip first header

nTracks = C{2}(1);

coords = cell(nTracks,1);
time = cell(nTracks,1);
coords_err = cell(nTracks,1);
z0 = cell(nTracks,1);
z0_err = cell(nTracks,1);
sigma = cell(nTracks,1);
sigma_err = cell(nTracks,1);
intensity = cell(nTracks,1);
intensity_err = cell(nTracks,1);
corrXY = cell(nTracks,1);
corrXY_err = cell(nTracks,1);

ii = 1;
txtline = fgetl(f);
while ischar(txtline) % end of file
    if strcmp(txtline(1:11),'First frame') % begin track
        trackdata = textscan(f,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'); 
        coords{ii} = [trackdata{5} trackdata{6}];
        time{ii} = trackdata{1};
        coords_err{ii} = [trackdata{12} trackdata{13}];
        z0{ii} = trackdata{8};
        z0_err{ii} = trackdata{11};
        sigma{ii} = [trackdata{3} trackdata{4}];
        sigma_err{ii} = [trackdata{10} trackdata{11}];
        intensity{ii} = trackdata{2};
        intensity_err{ii} = trackdata{9};
        corrXY{ii} = trackdata{7};
        corrXY_err{ii} = trackdata{14};
        ii = ii + 1;
    end
    txtline = fgetl(f);
end
fclose(f);

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
T = table(trackID,coords_err,z0,z0_err,sigma,sigma_err,...
    sigma_err,intensity,intensity_err,corrXY,corrXY_err);
obj.fitProps = [obj.fitProps T];
datasetID = strings(nTracks,1);
[~,filename,~] = fileparts(txtpath);
datasetID(:) = filename;
obj.fitProps.datasetID = datasetID; % overwrite datasetID with filepath
obj.filepath = txtpath;

obj.coordsZero = 1;
obj.dt = 1; % frame
% obj.dte = 1;

end
