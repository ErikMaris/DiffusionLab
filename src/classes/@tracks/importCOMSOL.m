function obj = importCOMSOL(txtpath,nDim,pixel_m,time_s)

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

% --- input dialog
if isempty(txtpath)
    [file,path] = uigetfile('*.txt', 'Multiselect', 'on');
    if isequal(file,0)
        % do nothing
    else
       txtpath = fullfile(path,file);
    end
end

if iscell(txtpath) % if is cell, the txt are batch processed within 'batchProcessTracks'
    obj = batchProcessTracks(txtpath,nDim,pixel_m,time_s);
else
    obj = importAsTracks(txtpath,nDim,pixel_m,time_s,txtpath);
end

end

function obj = batchProcessTracks(txtpath,nDim,pixel_m,time_s)

numtxt = length(txtpath);

% preallocation
obj = tracks.empty(numtxt,0);

for ii = 1:numtxt
    try
        obj(ii,1) = tracks.importCOMSOL(txtpath{ii},nDim,pixel_m,time_s);
    catch
        error(['Parsed .txt (' txtpath{ii} ') could not be imported by the function ''' mfilename '''.'])
    end
end

end

function obj = importAsTracks(filepath,nDim,pixel_m,time_s,txtpath)


%section wise export
out = importdata(filepath);

Ntracks = str2double(regexp(horzcat(out.textdata{:}), '(?<=Elements:[^0-9]*)[0-9]*\.?[0-9]+', 'match')); % find number of elements https://nl.mathworks.com/matlabcentral/answers/166837-pull-out-number-after-specific-string-in-txt-file
Ntimesteps = numel(strfind(fileread(filepath),'x (')); % get number of time steps
data = reshape(out.data(2*Ntracks+1:end),[Ntracks, numel(out.data(2*Ntracks+1:end))/Ntracks/Ntimesteps Ntimesteps]); % reshape
data = permute(data,[3 2 1]); % [timesteps column tracks] column  = [x y z R rho]

coords = squeeze(num2cell(data(:,1:nDim,:),[1 2])); % put coords in cell

time = cell(Ntracks,1); % create time cell array
time(:) = {(1:Ntimesteps)'};

for ii = 1:Ntracks % clear nan
    clear = any(isnan(coords{ii}),2);
    coords{ii}(clear,:) = [];
    time{ii}(clear) = [];
end

% --- create object

obj = tracks(time,coords,pixel_m,time_s);
obj.filepath = txtpath;

obj.coordsZero = 0;
obj.dt = 1; % frame


end
