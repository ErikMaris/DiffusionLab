function obj = importDoM(csvpath,falsePositivesFlag,pixel_m,time_s)
% importDoM Imports the DoM results table (.csv) table for tracks
% 
% Wrapper to read the .csv results table output in DoM. Compatible with DoM
% v.1.1.6. If the results table contains column 23-26, the results are of
% type 'tracks', otherwise of type 'fitresults'. Note that the false
% positives are saved in the results table and should be filtered if not
% taken into account for the track linking. 'csvpath' may be a cell array,
% however the .csv must all be of type 'tracks'.
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% obj = importDoMResultsTable(csvpath,falsePositivesFlag,pixel_m,time_s)
%   csvpath: path to the input file (.csv), empty gives dialog
%   falsePositivesFlag: include false positives?
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
    obj = batchProcessTracks(csvpath,falsePositivesFlag,pixel_m,time_s);
else
    % find the number of columns
    fid=fopen(csvpath);
    tline = fgetl(fid);
    fclose(fid);
    nC = length(find(tline==','))+1;
    if nC == 26
        filespec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
    elseif nC == 27 % in batch processing in imagej, an extra column is added on the first position containing 1:n, which has to be removed  
        filespec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
    else
        error(['Invalid number of columns in csv-file in function ''' mfilename '''.'])
    end
    fileID = fopen(csvpath);
    f = fread(fileID,'*char')';
    fclose(fileID);
    f = strrep(f,'Infinity','Inf'); % make compatible with matlab
    csvout = textscan(f, filespec,'delimiter',',','headerLines',1,'CollectOutput',true); % nan is "Inifinity"
    csvout = [csvout{:}];
    if nC == 27
        csvout = csvout(:,2:end);
    end
    obj = importAsTracks(csvout,falsePositivesFlag,pixel_m,time_s,csvpath);
end

end

function obj = batchProcessTracks(csvpath,falsePositivesFlag,pixel_m,time_s)

numcsv = length(csvpath);

% preallocation
obj = tracks.empty(numcsv,0);

for ii = 1:numcsv
    try
        obj(ii,1) = tracks.importDoM(csvpath{ii},falsePositivesFlag,pixel_m,time_s);
    catch
        error(['Parsed .csv (' csvpath{ii} ') are not all of the type ''tracks'' in ''' mfilename '''.'])
    end
end

end

function obj = importAsTracks(csvout,falsePositivesFlag,pixel_m,time_s,csvpath)

nTracks = max(csvout(:,23)); % last written trackID

% % Preallocate results, which will be an array of structures.
% fitProps = struct(...
%     'intensity', cell(nTracks,1),...
%     'state', cell(nTracks,1),...
%     'Rsq', cell(nTracks,1),...
%     'amp', cell(nTracks,1)...
%     );

coords = cell(nTracks,1);
time = cell(nTracks,1);
coords_err = cell(nTracks,1);
coords_err_nm = cell(nTracks,1);
amp = cell(nTracks,1);
amp_err = cell(nTracks,1);
z0 = cell(nTracks,1);
sigma = cell(nTracks,1);
sigma_nm = cell(nTracks,1);
sigma_err = cell(nTracks,1);
sigma_err_nm = cell(nTracks,1);
intensity = cell(nTracks,1);
state = cell(nTracks,1);
r2fit = cell(nTracks,1);

wb = waitbar(0,'Loading tracks...','Name','Please wait');

for ii = 1:nTracks
    waitbar(ii/nTracks,wb)
    indframe = csvout(:,23) == ii; % trackID is 3rd column
    if falsePositivesFlag % take all
        take = (indframe == 1);
    else % skip false positives
        take = (csvout(:,18) == 0) & (indframe == 1);
    end
    coords{ii} = [csvout(take,1) csvout(take,2)];
    time{ii} = csvout(take,3);
    coords_err{ii} = [csvout(take,5) csvout(take,7)]./(pixel_m.*1e9);
    coords_err_nm{ii} = [csvout(take,5) csvout(take,7)];
    amp{ii} = csvout(take,10);
    amp_err{ii} = csvout(take,11);
    z0{ii} = csvout(take,12);
    sigma{ii} = [csvout(take,14) csvout(take,16)]./(pixel_m.*1e9);
    sigma_nm{ii} = [csvout(take,14) csvout(take,16)];
    sigma_err{ii} = [csvout(take,15) csvout(take,17)]./(pixel_m.*1e9);
    sigma_err_nm{ii} = [csvout(take,15) csvout(take,17)];
    state{ii} = csvout(take,18);
    intensity{ii} = csvout(take,19);
    r2fit{ii} = csvout(take,21);
end

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
trackID = (1:nTracks)';
T = table(trackID,coords_err,coords_err_nm,amp,amp_err,z0,sigma,sigma_nm,...
    sigma_err,sigma_err_nm,intensity,state,r2fit);
obj.fitProps = T;
obj.filepath = csvpath;

% obj.fitProps.coords_err = coords_err;
% obj.fitProps.coords_err_nm = coords_err_nm;
% obj.fitProps.amp = amp;
% obj.fitProps.amp_err = amp_err;
% obj.fitProps.z0 = z0;
% obj.fitProps.sigma = sigma;
% obj.fitProps.sigma_nm = sigma_nm;
% obj.fitProps.sigma_err = sigma_err;
% obj.fitProps.sigma_err_nm = sigma_err_nm;
% obj.fitProps.intensity = intensity;
% obj.fitProps.state = state;
% obj.fitProps.r2fit = r2fit;

obj.coordsZero = 1;
obj.dt = 1; % frame
% obj.dte = 1;

close(wb)

end

%   Expected input csv:
% 	public static final int Col_X=0;
% 	public static final int Col_Y=1;
% 	public static final int Col_FrameN=2;
% 	public static final int Col_Xnm=3;
% 	public static final int Col_loc_errX=4;
% 	public static final int Col_Ynm=5;
% 	public static final int Col_loc_errY=6;
% 	public static final int Col_Znm=7;	
% 	public static final int Col_loc_errZ=8;
% 	public static final int Col_AmplFit=9;
% 	public static final int Col_Amp_error=10;
% 	public static final int Col_BGfit=11;
% 	public static final int Col_BGfit_error=12;
% 	public static final int Col_SD_X=13;
% 	public static final int Col_SD_X_err=14;
% 	public static final int Col_SD_Y=15;
% 	public static final int Col_SD_Y_err=16;
% 	public static final int Col_Fp=17;	
% 	public static final int Col_IntegrInt=18;
% 	public static final int Col_SNR=19;
% 	public static final int Col_chi=20;
% 	public static final int Col_IterN=21;
% 	
% 	public static final int Col_TrackID=22;
% 	public static final int Col_ParticleID=23;
% 	public static final int Col_TrackLength=24;
% 	public static final int Col_TrackReverseN=25;
