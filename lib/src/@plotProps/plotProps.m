classdef plotProps
    % PLOTPROPS This class stores and handles properties for plotting 

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
    
    % plotting
    properties % make transient?
        colormap = 'viridis'; % string // stores the name of the colormap https://nl.mathworks.com/help/matlab/ref/colormap.html
        lineWidth = 1; % line thickness
    end
    
    methods
        function obj = plotProps(varargin)
            %plotting Construct an instance of this class
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
        
        function obj = set.lineWidth(obj,value)
            validateattributes(value,{'numeric'},{'scalar','nonnegative','nonzero'},mfilename)
            obj.lineWidth = value;
        end
        
        function obj = set.colormap(obj,value)
            validateattributes(value,{'char','string'},{'nonempty'},mfilename)
            obj.colormap = value;
        end
        
        function cmap = getColormap(obj,n)
            % function gets colormap RGB values
            if nargin < 2
                cmap = feval(obj.colormap);
            else
                cmap = feval(obj.colormap,n);
            end
        end
        
    end
end