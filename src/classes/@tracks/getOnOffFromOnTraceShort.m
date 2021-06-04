function [on,off] = getOnOffFromOnTraceShort(obj, discardEndsFlag, indices)
% [on,off] = getOnOffFromOnTraceShort(tracks,discardEndsFlag)
%--------------------------------------------------------------------------
% Description
%
% Function extracts the on and off times from the onTraceShort.
%
% Literature:
% 
%--------------------------------------------------------------------------
% Necessary Inputs (name/data type)
% 
% tracks/struct: tracks structure containing the following fields
%   onTraceShort
% discardEndsFlag/boolean: exclude on/off values located at the first or
%   last position of the onTraceShort? The values at the begin and end of 
%   the COMPLETE onTrace are discarded to obtain values only enclosed by 
%   measured events
%
%--------------------------------------------------------------------------
% Outputs (name/data type)
%
% on/cell: on{numMovies}{numTracks} values of on time lengths in
%   chronological order. 
% off/cell: off{numMovies}{numTracks} values of on time lengths in
%   chronological order.
%
%--------------------------------------------------------------------------
% Variable Inputs (flag/ data type /(default)):
% 
% 
%--------------------------------------------------------------------------
% Dependencies
%
%

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

if numel(obj) > 1
    error('SC_blinking:multipleObjParsed', ...
        'Function does not accept object arrays.')
end

if nargin < 2
    discardEndsFlag = true;
end

if nargin < 3
    indices = 1:obj.nTracks;
end

n = numel(indices);

on = cell(n,1);
off = cell(n,1);

% --- get on and off traces and remove begin and end off states
for ii = 1:n
    jj = indices(ii);
    onTraceShort = (obj.onTrace{jj})'; % ' is trick to allow concatenation with empty arrays - which are row vectors
    if discardEndsFlag
        off{jj,1} = onTraceShort(3:2:end-1); % odd vector
        if onTraceShort(1) == 0 % real first value is index 2, which has to be skipped
            on{jj,1} = onTraceShort(4:2:end-1); % even vector
        else
            on{jj,1} = onTraceShort(2:2:end-1); % even vector
        end
    else
        on{jj,1} = onTraceShort(2:2:end);  % even vector
        off{jj,1} = onTraceShort(1:2:end); % odd vector
        if off{jj,1}(1) == 0 % remove odd zero
            off{jj,1} = off{jj,1}(2:end);
        end
    end
end

on = horzcat(on{:})';
off = horzcat(off{:})';

