classdef (Abstract) DEstimator
    % DestModel Abstract subclass for defining diffusion estimation objects
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
        estName % name of the diffusion estimator
    end
    
    properties (Abstract, Constant, Hidden)
        % parameters for GUI
        inputdlgPrompt % inputdlg prompt
        inputdlgVariables % names of variables corresponding to prompt
        plotdlgNames; % plot function names
        plotdlgDescription; % description for gui
        plotdlgType; % track, population, bootstrap
    end
    
    methods(Abstract)
        obj = compute(obj);
    end
    
    methods
        
        function varargout = segment(obj,bool)
            % SEGMENT Segments the obeject based on user-input
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
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University


            nObj = size(obj,1);
            bool = logical(bool);

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
                    if isstruct(obj(ii).(objProps{jj})) 
                        if size(obj(ii).(objProps{jj}),1) == nBool % first do deal with structure arrays
                            obj(ii).(objProps{jj}) = obj(ii).(objProps{jj})(bool,:);
                        else % deal with inhomogeneous structures
                            obj(ii).(objProps{jj}) = obj.segmentStruct(obj(ii).(objProps{jj}),bool); % recursive thresholding // needs a special function dealing with structures, since these dont have a segment method
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
        
        function obj = inputdlgVars(obj,dims)
            % INPUTDLGVARS Allows the user to set variables via a GUI
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % obj = inputdlgVars(obj)
            %   Asks for user-input on variables specified in diffusion
            %   estimator
            % obj = inputdlgVars(obj,dims)
            %   allows the user to set the size of the GUI dialog in dims
            %   [height width]
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % 
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            if nargin < 2
                dims = [1 70];
            end
            
            % set variables with GUI
            prompt =  obj.inputdlgPrompt; % inputdlg prompt
            vars = obj.inputdlgVariables;
            if isempty(prompt) || isempty(vars)
                msg = sprintf('No options available for %s',obj.estName);
                msgbox(msg,'Notification')
                return
            end
            nVars = numel(vars);
            definput = cell(nVars,1);
            type = NaN(nVars,1);
            for ii = 1:nVars
                if islogical(obj.(vars{ii}))
                    definput{ii} = obj.(vars{ii});
                    type(ii) = 1;
                elseif isnumeric(obj.(vars{ii})) && ismatrix(obj.(vars{ii}))
                    definput{ii} = mat2str(obj.(vars{ii}));
                    type(ii) = 2;
                elseif isnumeric(obj.(vars{ii}))
                    definput{ii} = num2str(obj.(vars{ii}));
                    type(ii) = 3;
                elseif isobject(obj.(vars{ii}))
                    definput{ii} = class(obj.(vars{ii}));
                    type(ii) = 4;
                else
                    definput{ii} = obj.(vars{ii});
                    type(ii) = 5;
                end
            end
            definput = string(definput);
            dlgtitle = sprintf('Input variables for %s',obj.estName);
            uservalues = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(uservalues)
                return % user pressed cancel
            end
            for ii = 1:nVars
                switch type(ii)
                    case 1
                        if ischar(uservalues{ii})
                            obj.(vars{ii}) = strcmpi(uservalues{ii}, 'true'); % case insensitive
                        else % isnumeric
                            obj.(vars{ii}) = logical(str2double(uservalues{ii}));
                        end
                    case 2
                        obj.(vars{ii}) = str2num(uservalues{ii});
                    case 3
                        obj.(vars{ii}) = str2num(uservalues{ii});
                    case 4
                        try
                            obj.(vars{ii}) = eval(uservalues{ii});
                        catch
                            error('Input is not recognised.')
                        end
                    case 5
                        obj.(vars{ii}) = uservalues{ii};
                end
            end
        end
        
        function plotdlg(obj,parent,trackIdx)
            % PLOTDLG Wrapper to for plotting via the GUI
            %
            % Assumes standard settings with type=tracks 
            % obj.plot(parent,[],trackID).
            % 
            % -------------------------------------------------------------
            % -                         USAGE                             -
            % -------------------------------------------------------------
            % 
            % plotdlg(obj,parent,trackIdx)
            %   trackIdx sets the track to be plotted if required by the
            %   plot function.
            
            % -------------------------------------------------------------
            % -                         HISTORY                           -
            % -------------------------------------------------------------
            % 
            % 
            % 
            % -------------------------------------------------------------
            % Copyright (C) 2019 J.J. Erik Maris
            % Inorganic Chemistry & Catalysis, Utrecht University
            
            if isempty(obj.plotdlgNames)
                msgbox('No plots available for this diffusion estimator','Notification')
                return
            end
            
            if ~parent.Dest_valid
                msgbox('Please (re)compute diffusion','Notification')
                return
            end
            
            [idx,tf] = listdlg('PromptString','Select a plot:',...
                           'SelectionMode','single',...
                           'ListString',obj.plotdlgDescription,...
                           'ListSize',[250 300]);
            
            if tf == 0 || isempty(idx)
                return % user pressed cancel
            end
                
            switch obj.plotdlgType{idx}
                case 'population'
                    figure
                    obj.(obj.plotdlgNames{idx})(parent)
                case 'bootstrap'
                    figure
                    obj.(obj.plotdlgNames{idx})(parent)
                case 'track'
                    figure
                    obj.(obj.plotdlgNames{idx})(parent,[],trackIdx)
                otherwise
                    error('Plot type %s is not recognised in %s.',obj.plot4GUItype{plotIdx},mfilename)
            end
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
                 error('Not a valid attribute name')
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
            
    end

    methods (Static)
        
        [x,fval,exitflag,output]=fminsearchbnd(fun,x0,LB,UB,options,varargin)
        
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
                        S.(fn{jj}) = Dest.DestModel.segmentStruct(S.(fn{jj}),bool); % recursive thresholding
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
        
    end
    
end

