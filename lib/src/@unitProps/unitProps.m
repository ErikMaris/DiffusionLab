classdef unitProps < classHandling
    %UNITPROPS Superclass stores and converts units
    %
    % Supported units are given in: 
    % https://nl.mathworks.com/help/symbolic/units-list.html
    % always give units in BASE UNITS, otherwise code fails
    %
    % Assumes that all data is provided in the same units stored in this
    % object.
    
    %
    % KNOWN ISSUES:
    % - saving an variable which contains a symunit is orders of mangnitude
    % slower than saving a number. Try to avoid saving symunits in the
    % dictionary of the unitProps.
    %
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

    
    properties (SetAccess = protected)
        unitSystem = baseUnits('SI'); % stores the base unit system for display; modify this to change the display. Hardcoded in function resetUnitSystem
        dataUnits; % containers.Map , dictionary
    end

    methods
        function obj = unitProps(varargin)
            %UNITPROPS constructs an instance of this class
            
            obj.dataUnits = containers.Map('UniformValues',false);
            
            if nargin > 1
                obj = obj.setUnit(varargin);
            end
                
        end
        
        function obj = set.unitSystem(obj,value)
            validateattributes(value,{'sym'},{'vector'},mfilename,'unitSystem')
            obj.unitSystem = value;
        end
        
        function obj = set.dataUnits(obj,value)
            validateattributes(value,{'containers.Map'},{},mfilename,'dataUnits')
            obj.dataUnits = value;
        end
        
        function obj = setUnit(obj,varargin)
            
            nIn = numel(varargin);
            if mod(nIn,2) ~= 0
                error('Variables can only be parsed as Name-Value pairs')
            end
            keySet = {varargin{1:2:end}}; % input variable name
            cellfun(@(x) validateattributes(x,{'char','string'},{'nonempty'}),keySet)
            valueSet = {varargin{2:2:end}}; % it's value as symunit or string
%             cellfun(@(x) obj.validateattributesGUI(x,{'sym'},{'nonempty'}),valueSet)
            for ii = 1:numel(obj)
                for jj = 1:numel(keySet)
                    if isa(valueSet{jj},'sym')
                        obj(ii).dataUnits(keySet{jj}) = valueSet{jj};
                    elseif ischar(valueSet{jj}) || isstring(valueSet{jj})
                        obj(ii).dataUnits(keySet{jj}) = str2symunit(valueSet{jj});
                    else
                        error('Data type not recognized.');
                    end
                end
            end
            
%             if numel(obj) > 1
%                 obj = obj.setArrayProperty('dataUnits',obj(1).dataUnits);
%                 obj = obj.setArrayProperty('unitSystem',obj(1).unitSystem);
%             end
            
        end
        
        function varargout = getUnit(obj,varargin)
            
            varargout = cell(nargin-1,1);
            for ii = 1:nargin-1
                varargout{ii} = obj(1).dataUnits(varargin{ii});
            end
            
        end
        
        function varargout = getUnitFactor(obj,dataUnit)
            
            if numel(obj) > 1
                error('Only one object can be parsed to ''getUnitFactor''.')
            end
            
            if isempty(dataUnit)
                val = 1;
                name = [];
            elseif ~isa(dataUnit,'sym') % input is string or char
                % https://nl.mathworks.com/help/matlab/matlab_prog/matlab-operators-and-special-characters.html
                units = unique(split(dataUnit,{'+','-','.*','*','./','/','.\','\','.^','^','.''',''''}));
                take = cellfun(@(x) str2double(x),units); % find nan to filter out multiplication factors
                units = units(isnan(take)); % and take
                junits = join(units);
                eval(['[' junits{:} ']=obj.getUnit(units{:});']);
                u = rewrite(eval(dataUnit),obj.unitSystem);
                [val, name] = separateUnits(u);
            else
                u = rewrite(dataUnit,obj.unitSystem);
                [val, name] = separateUnits(u);
            end
            
            varargout{1} = double(val);
            if nargout > 1
                l = symunit2str(name);
                l = strrep(l,'mcm',[char(181) 'm']);
                if isempty(l); l = []; end % symunit2str gives empty cell array on empty name. make emtpy double
                varargout{2} = l;
            end
            
        end
        
        function obj = subsUnitSystem(obj,old,new)
            % change unit in unit system from old to new. Old and new can
            % be either a symunit or character
            % or give subsUnitSystem(obj,new)
            
            u = symunit;
            
            if nargin == 3
                if isa(old,'sym') % symunit
                elseif ischar(old) || isstring(old)
                    old = u.(old);
                else
                    error('Old unit must either be a symunit, character, or string')
                end

                if isa(new,'sym') % symunit
                elseif ischar(new) || isstring(new)
                    new = u.(new);
                else
                    error('New unit must either be a symunit, character, or string')
                end
            else
                if isa(old,'sym') % symunit
                elseif ischar(old) || isstring(old)
                    new = u.(old); % this is expected as new input
                else
                    error('New unit must either be a symunit, character, or string')
                end
                
                uSys = obj(1).unitSystem;
                
                % check which of the entries in the unit system is
                % compatible
                bool = checkUnits(uSys==new,'Compatible');
                if sum(bool) ~= 1
                    error('Multiple compatible units available. Cannot change the length unit.')
                end
                
                % get old unit
                old = uSys(bool);
            end
            
            % do substitution
            for ii = 1:numel(obj)
                obj(ii).unitSystem = subs(obj(ii).unitSystem,old,new);
            end
            
        end
        
        function obj = resetUnitSystem(obj)
            
            for ii = 1:numel(obj)
                obj(ii).unitSystem = baseUnits('SI');
            end
            
        end
    end
    
end
        