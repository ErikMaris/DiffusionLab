classdef classHandling
    % CLASSHANDLING This class deals with general operations on objects

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
    %classHandling Superclass deals with general properties and
    %functionality
    % It manages the way objects are displayed in a GUI and the output they
    % give.
    
    properties % (SetAccess = private)
        verbose = false; % logical // sets whether a function disps additional info
%         GUImode = false; % logical // sets whether the object displays info as being in an GUI
    end
    
    % Constructor
    methods
        function obj = classHandling(verbose)
            %generalClassHandling Construct an instance of this class
            % classHandling constructs an instance of this class.
            %   classHandling() gives the default values 'false' for 
            %   both 'verbose' and 'GUImode'.
            %   classHandling(verbose) allows the user to set 
            %   the boolean values for the variables.
            
            if nargin > 0
                obj.verbose = verbose;
%                 obj.GUImode = varargin{2};
            end
            
        end 
    end
    
    % set functions
    methods
        function obj = set.verbose(obj,value)
            validateattributes(value,{'logical'},{'scalar'},mfilename,'verbose')
            obj.verbose = value;
        end
    end
    
    % Public functions
    
    methods (Static)
        function objOut = open(openpath,classname)
            %converts structure s to an object of class classname.
            %assumes classname has a constructor which takes no arguments
            
            
            % --- input dialog
            if isempty(openpath)
                [file,path] = uigetfile('*.mat');
                if isequal(file,0)
                    % do nothing
                else
                   openpath = fullfile(path,file);
                end
            end
            
            fprintf('Loading %s\n',openpath)
            
            o = open(openpath);
            
            loadObj = o.saveObj;
            
            objOut = eval([classname '.empty([numel(loadObj) 0])']);
            objOut(:) = eval(classname);  %create object
            
            for ii = 1:numel(loadObj)
                for fn = fieldnames(loadObj)'    %enumerat fields
                  try
                      objOut(ii).(fn{1}) = loadObj(ii).(fn{1});   %and copy
                  catch
%                       warning('Could not copy field %s', fn{1});
                  end
                end
            end
            
            objOut = reshape(objOut,size(loadObj));
            
        end
        
        function objOut = openStruct(S,classname)
           objOut = eval([classname '.empty([numel(S) 0])']);
            objOut(:) = eval(classname);  %create object
            
            for ii = 1:numel(S)
                for fn = fieldnames(S)'    %enumerat fields
                  try
                      objOut(ii).(fn{1}) = S(ii).(fn{1});   %and copy
                  catch
%                       warning('Could not copy field %s', fn{1});
                  end
                end
            end
            
            objOut = reshape(objOut,size(S));
        end
        
    end
    
    methods
        
        function saveAsStruct(obj,savepath)
            if nargin < 2
                savepath = [];
            end
            
            % --- input dialog
            if isempty(savepath)
                [file,path] = uiputfile('*.mat');
                if isequal(file,0)
                    % do nothing
                else
                   savepath = fullfile(path,file);
                end
            end
            
            fprintf('Saving %s\n',savepath)
            
            warning('off','MATLAB:structOnObject')
            
            saveObj(1) = struct(obj(1));

            if numel(obj) > 1
                saveObj = repmat(saveObj,(size(obj)));
                for ii = 2:numel(obj)
                    saveObj(ii) = struct(obj(ii));
                end
            end
            warning('on','MATLAB:structOnObject')
            
            save(savepath,'saveObj')%,'-v7.3');
            
            fprintf('Done!\n')
            
        end
        
        function save(obj,savename,savepath)
            
            if nargin < 2
                savename = 'saveObj';
            end
            
            if nargin < 3
                savepath = [];
            end
            
            % --- input dialog
            if isempty(savepath)
                [file,path] = uiputfile('*.mat');
                if isequal(file,0)
                    % do nothing
                else
                   savepath = fullfile(path,file);
                end
            end
            
            fprintf('Saving %s\n',savepath)
            
            eval([savename ' = obj;'])
            
            save(savepath,savename,'-v7.3');
            
            fprintf('Done!\n')
            
        end
        
        function props = findNonEmptyProperties(obj)
            
            props = {};
            for ii = 1:numel(obj)
                for fn = fieldnames(obj)'    %enumerat fields
                    if ~isempty(obj(ii).(fn{1}))
                        props = [props fn{1}];
                    end
                end
            end
            
            props = props';
            
        end
        
        function props = findEmptyProperties(obj)
            
            props = {};
            for ii = 1:numel(obj)
                for fn = fieldnames(obj)'    %enumerat fields
                    if isempty(obj(ii).(fn{1}))
                        props = [props fn{1}];
                    end
                end
            end
            
            props = props';
            
        end
        
        function cl_out = findAttrValue(obj,attrName,varargin)
            % example usage:
            % - obj.findAttrValue('Constant')
            % - obj.findAttrValue('Dependent')
            % Options https://nl.mathworks.com/help/matlab/matlab_oop/property-attributes.html
           if ischar(obj)
              mc = meta.class.fromName(obj);
           elseif isobject(obj)
              mc = metaclass(obj);
           end
           ii = 0; numb_props = length(mc.PropertyList);
           cl_array = cell(1,numb_props);
           for  c = 1:numb_props
              mp = mc.PropertyList(c);
              if isempty (findprop(mp,attrName))
                 obj.error('Not a valid attribute name')
              end
              attrValue = mp.(attrName);
              if attrValue
                 if islogical(attrValue) || strcmp(varargin{1},attrValue)
                    ii = ii + 1;
                    cl_array(ii) = {mp.Name};
                 end
              end
           end
           cl_out = cl_array(1:ii);
        end
        
        function fieldVals = getProp(obj,fieldName,noCellFlag)
            % Deals with referencing fields within a structure within an
            % object. Output is cell(nObj,1) if noCellFlag = false,
            % otherwise if noCellFlag = true only the first object is taken
            % and the output is given as is, i.e. without cell.
            
            if nargin < 3
                noCellFlag = false;
            end
            
            fieldName = split(fieldName,'.');
            if numel(fieldName) > 1
                fieldVals = cell(numel(obj),1);
                for ii = 1:numel(obj) % I don't see a way to avoid a loop here :-(
                    try
                        fieldVals{ii} = obj(ii).(fieldName{1});
                        for jj = 2:numel(fieldName)
                            fieldVals{ii} = fieldVals{ii}.(fieldName{jj});
                            if isempty(fieldVals{ii}) % if it doesn't exist return empty
                                fieldVals = [];
                                return
                            end
                        end
                    catch
                        try
                            fieldVals{ii} = eval(['obj' num2str(ii) '.' fieldName{1} ]); 
                        catch
%                             error('SC_classHandling:getProp:badInput', ...
%                                 ['''' [fieldName{:}] ''' is not valid as field name or function.'])
                            fieldVals{ii} = [];
                        end
                    end
                end
            else
                fieldVals = {obj.(fieldName{1})}; % make cell
                fieldVals = fieldVals(:); % make column vector
            end
            
            if noCellFlag % could be faster MEM-wise if cat of all fieldvals is prevented
                fieldVals = fieldVals{1};
            end
            
        end
        
        % functions dealing with display
        
%         function errorGUI(obj,varargin)
%             % ERROR allows dynamic error display dependent on 'GUImode'.
%             % displays the error in the log window (GUImode = off) or in an
%             % error dialog box (GUImode = on)
% 
%             if obj(1).GUImode % take first object for leading properties if multiple parsed
%                 try
%                     error(varargin{:})
%                 catch e %e is a MException struct
%                     errordlg(e.message)
%                     rethrow(e) % required to stop code
%                 end
%             else
%                 error(varargin{:})
%             end
%         end
%         
%         function validateattributesGUI(obj,varargin)
%             % VALIDATEATTRIBUTES allows dynamic error display dependent on 
%             % 'GUImode'. Displays the error in the log window (GUImode = 
%             % off) or in an error dialog box (GUImode = on)
%             
%             if obj(1).GUImode % take first object for leading properties if multiple parsed
%                 try
%                     validateattributes(varargin{:})
%                 catch e %e is a MException struct
%                     errordlg(e.message)
%                     rethrow(e) % required to stop code
%                 end
%             else
%                 validateattributes(varargin{:})
%             end
%         end
%         
%         function varargout = validatestringGUI(obj,varargin)
%             % VALIDATESTRING allows dynamic error display dependent on 
%             % 'GUImode'. Displays the error in the log window (GUImode = 
%             % off) or in an error dialog box (GUImode = on)
% 
%             if obj(1).GUImode % take first object for leading properties if multiple parsed
%                 try
%                     validatestring(varargin{:})
%                 catch e %e is a MException struct
%                     errordlg(e.message)
%                     rethrow(e) % required to stop code
%                 end
%             else
%                 varargout{1} = validatestring(varargin{:});
%             end
%         end
%         
%         function obj = setVerbose(obj,value)
%             % SETVERBOSE allows the user to set the value of 'verbose'.
%             %   SETVERBOSE(obj) switches the value of 'verbose'.
%             %   SETVERBOSE(obj,value) sets the value of 'verbose' to either
%             %   'true' or 'false'.
%             
%             
%             if nargin < 2
%                 if obj.verbose == true
%                     obj.verbose = false;
%                 else
%                     obj.verbose = true;
%                 end
%             else
%                 switch value
%                     case true
%                         obj.verbose = true;
%                     case false
%                         obj.verbose = false;
%                     otherwise
%                         obj.error('classHandling:setVerbose:incorrectType',...
%                             ['Input ''' value ''' is not ' ...
%                             'recognised by function ''setVerbose'' in '''...
%                             mfilename '''.'])
%                 end
%             end
%         end
%                 
%         function obj = setGUImode(obj,value)
%             % SETGUIMODE allows the user to set the value of 'GUImode'.
%             %   SETGUIMODE(obj) switches the value of 'GUImode'.
%             %   SETGUIMODE(obj,value) sets the value of 'GUImode' to either
%             %   'true' or 'false'.
%             
%             if nargin < 2
%                 if obj.GUImode == true
%                     obj.GUImode = false;
%                 else
%                     obj.GUImode = true;
%                 end
%             else
%                 switch value
%                     case true
%                         obj.GUImode = true;
%                     case false
%                         obj.GUImode = false;
%                     otherwise
%                         obj.error('classHandling:setGUImode:incorrectType',...
%                             ['Input ''' value ''' is not ' ...
%                             'recognised by function ''setGUImode'' in '''...
%                             mfilename '''.'])
%                 end
%             end
%         end
        
        function obj = setArrayProperty(obj,property,value)
            % set the property to value for all objects in obj array
            [obj.(property)] = deal(value);
        end
        
%         obj = threshold(obj,thresholdField,minThreshold,maxThreshold);
%         keepBoolean = thresholdBoolean(obj,thresholdField,minThreshold,maxThreshold);
        
    end
    
end

