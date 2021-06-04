function varargout = segment(obj,bool)
% SEGMENT Segments the object based on user-input
% 
% Assumes that the to-be-segmented data is stored in columns
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% [obj=bool] = segment(obj,bool)
%   segments the recursively object based on a bool with true =
%   keep and false = throw
%
% [obj=bool,obj~=bool] = segment(obj,bool)
%   segments the recursively object based on a bool with true =
%   first object and false = second object

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

nObj = size(obj,1);
bool = logical(bool);

if numel(unique(bool)) ~=2
    error('Segmented population should contain at least one track')
end

for ii = 1:nObj
    if nargout > 1
        varargout{2} = obj.segment(~bool); % compute the ~bool
    end
    % --- decide which ones to keep
    % --- loop over all fields and take fields specified in 'keepBoolean'
    objProps = properties(obj(ii));
    objPropsSkip = [obj(ii).findAttrValue('Dependent') obj(ii).findAttrValue('Constant')]; % find dependent and constant object properties
    objProps = objProps(~ismember(objProps,objPropsSkip)); % remove dependent and constant properties
    nBool = numel(bool);
    for jj = 1:numel(objProps)
        % --- field is thresholded if it has a number of rows equal to nTracks
        if strcmp('features',objProps{jj}) && ~isempty(obj(ii).features)
            T = sortrows(obj(ii).features,'trackID'); % bool assumes sorted columns
            obj(ii).(objProps{jj}) = T(bool,:); % take desired rows
        elseif isstruct(obj(ii).(objProps{jj})) 
            if size(obj(ii).(objProps{jj}),1) == nBool % first do deal with structure arrays
                obj(ii).(objProps{jj}) = obj(ii).(objProps{jj})(bool,:);
            else % deal with inhomogeneous structures
                obj(ii).(objProps{jj}) = segmentStruct(obj(ii).(objProps{jj}),bool); % recursive thresholding // needs a special function dealing with structures, since these dont have a segment method
            end
        elseif isnumeric(obj(ii).(objProps{jj})) && size(obj(ii).(objProps{jj}),1) == nBool
            obj(ii).(objProps{jj}) = obj(ii).(objProps{jj})(bool, :);
        elseif iscell(obj(ii).(objProps{jj})) && size(obj(ii).(objProps{jj}),1) == nBool
            obj(ii).(objProps{jj}) = obj(ii).(objProps{jj})(bool, :);
        elseif isobject(obj(ii).(objProps{jj}))
            if ismethod(obj(ii).(objProps{jj}),'segment') % if class supports segmentation, perform thresholding
                obj(ii).(objProps{jj}) = obj(ii).(objProps{jj}).segment(bool); % make sure subclasses have a threshold function
            else
                if size(obj(ii).(objProps{jj}),1) == nBool % check size and segment accordingly
                    obj(ii).(objProps{jj}) = obj(ii).(objProps{jj})(bool, :);
                end % else do nothing
            end
        else
            % skip
        end
    end
    varargout{1} = obj;
end

end

function S = segmentStruct(S,bool)
% SEGMENTSTRUCT Segments a structure based on user-input
% 
% Helper function for segment. Does not allow recursive 
% segmentation of objects
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% [S=bool] = segmentStruct(S,bool)
%   segments the recursively structure based on a bool with 
%   true = keep and false = throw
%

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% 
% 
% -------------------------------------------------------------
% Copyright (C) 2019 J.J. Erik Maris
% Inorganic Chemistry & Catalysis, Utrecht University
fn = fieldnames(S);
nBool = numel(bool);
for jj = 1:numel(fn)
    % --- field is thresholded if it has a number of rows equal to nTracks
    if isstruct(S.(fn{jj})) 
        if size(S.(fn{jj})) == nBool % first do deal with structure arrays
            S.(fn{jj}) = S.(fn{jj})(bool,:);
        else % deal with inhomogeneous structures
            S.(fn{jj}) = segmentStruct(S.(fn{jj}),bool); % recursive thresholding
        end
    elseif isnumeric(S.(fn{jj})) && size(S.(fn{jj}),1) == nBool
        S.(fn{jj}) = S.(fn{jj})(bool, :);
    elseif iscell(S.(fn{jj})) && size(S.(fn{jj}),1) == nBool
        S.(fn{jj}) = S.(fn{jj})(bool, :);
    elseif isobject(S.(fn{jj}))
        if ismethod(S.(fn{jj}),'segment') % if class supports segmentation, perform thresholding
            S.(fn{jj}) = S.(fn{jj}).segment(bool); % make sure subclasses have a threshold function
        else
            if size(S.(fn{jj}),1) == nBool % check size and segment accordingly
                S.(fn{jj}) = S.(fn{jj})(bool, :);
            end
        end
    else
        % skip
    end
end

end