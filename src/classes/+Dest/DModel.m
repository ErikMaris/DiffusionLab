classdef (Abstract) DModel
    % DestModel Abstract subclass for defining diffusion model objects
    %
    % Class is supposed to be stored as a variable in another "parent"
    % class which provides the necessary input for the analysis.
    
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
    
    properties (Abstract, Constant)
        name % name of the diffusion estimator
    end
    
    methods (Abstract, Static)
%         Dest = computeResults(Dest,parent)
%         Dest = computeMeanResults(Dest,parent)
        report(Dest,parent)
    end
end

