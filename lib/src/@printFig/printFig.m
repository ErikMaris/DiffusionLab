classdef printFig
    % PRINTFIG This class deals with figure printing
    %   GUI available 'ExportFigures_App'

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
    
    properties % make transient?
        savepath = ''; % char, string // savepath; if name then the current MATLAB path will be taken
        projectname = ''; % char, string // name appended before file name
%         overwriteFlag = false; % logical // overwriting permitted? if not, a number is appended
        fh = []; % figure handle object
        paperPosition = [0 0 10 8]; % vector // last two values specify the size of the printed figure [0 0 W H]
        printFlag = true; % logical // save when printFig is called?
        plotExt = {'-dsvg','-dmeta','fig','-dpng'}; % cell // save figure with extensions (.fig & as defined in https://nl.mathworks.com/help/matlab/ref/print.html)
        paperUnits = 'centimeters'; % specify units for print size
        printDPI = 300; % scalar // resolution in DPI for printing
        fontName = 'Arial';
    end
    
    methods
        function obj = printFig(varargin)
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
        
        function obj = set.savepath(obj,value)
            if isa(value,'char') || isa(value,'string')
            else
                error('Expected type ''char'' or ''string'', but found a ''%s'' instead.',...
                    class(value))
            end
            obj.savepath = value;
        end
        
        function obj = set.fh(obj,value)
            validateattributes(value,{'matlab.ui.Figure','numeric'},{},mfilename)
            obj.fh = value;
        end
        
        function obj = set.paperPosition(obj,value)
            validateattributes(value,{'numeric'},{'vector','numel',4},mfilename)
            obj.paperPosition = value;
        end
        
        function obj = set.printFlag(obj,value)
            validateattributes(value,{'logical'},{'nonempty'},mfilename)
            obj.printFlag = value;
        end
        
        function obj = set.plotExt(obj,value)
            if iscell(value)
                for ii = 1:numel(value)
                    validatestring(value{ii},{'-dpdf','-deps','-depsc','-deps2','-depsc2',...
                        '-dmeta','-dsvg','-dps','-dpsc','-dps2','-dpsc2',...
                        '-djpeg','-dpng','-dtiff','-dtiffn','-dmeta',...
                        '-dbmpmono','-dbmp','-dbmp16m','-dbmp256','-dhdf','-dpbm',...
                        '-dpbmraw','-dpcxmono','-dpcx24b','-dpcx256','-dpcx16','-dpgm',...
                        '-dpgmraw','-dppm','-dppmraw','fig','svg'},mfilename,'plotExt',ii);
                end
            else
                validatestring(value,{'-dpdf','-deps','-depsc','-deps2','-depsc2',...
                            '-dmeta','-dsvg','-dps','-dpsc','-dps2','-dpsc2',...
                            '-djpeg','-dpng','-dtiff','-dtiffn','-dmeta',...
                            '-dbmpmono','-dbmp','-dbmp16m','-dbmp256','-dhdf','-dpbm',...
                            '-dpbmraw','-dpcxmono','-dpcx24b','-dpcx256','-dpcx16','-dpgm',...
                            '-dpgmraw','-dppm','-dppmraw','fig','svg'},mfilename);
            end
            obj.plotExt = value;
        end
        
        function obj = set.paperUnits(obj,value)
            validatestring(value,{'pixels','normalized','inches','centimeters','points','characters'});
            obj.paperUnits = value;
        end
        
        function obj = set.printDPI(obj,value)
            validateattributes(value,{'numeric'},{'nonzero','nonempty'});
            obj.printDPI = value;
        end
        
        function cmap = getColormap(obj,n)
            % function gets colormap RGB values
            if nargin < 2
                cmap = feval(obj.colormap);
            else
                cmap = feval(obj.colormap,n);
            end
        end
        
        function print(obj,savename,printTitleFlag)
            % function saves figure; is the OOP variant of the function
            % 'printFig'.

            if ~obj.printFlag
                return
            end
            
            if nargin < 3
                printTitleFlag = true;
            end

            if isempty(obj.fh)
                Fh = gcf; % get current figure information
            else
                Fh = obj.fh;
            end

            if ischar(obj.plotExt) || isstring(obj.plotExt)
                Extensions = {obj.plotExt};
            else
                Extensions = obj.plotExt;
            end
            
            Fh.PaperUnits = obj.paperUnits; % specify units for print size
            Fh.PaperPosition = obj.paperPosition; % size of graph for printing
            Fh.PaperPositionMode = 'manual'; % from the website don't know
            

            thisSavepath = fullfile(obj.savepath,[obj.projectname savename]);
            
            
%             if ~obj.overwriteFlag
%                 % check if the file already exists
%                 if exist(thisSavepath, 'file')
%                     ii = 1;
%                     newSavepath = thisSavepath;
%                     % keep trying new appendices
%                     while exist(newSavepath, 'file')
%                         % safety to escape while loop
%                         if ii > 1e6
%                             error('Cannot find an unique filename.\nSaving of %s aborted.',thisSavepath)
%                         end
%                         newSavepath = [thisSavepath '_' num2str(ii)];
%                         ii = ii + 1;
%                     end
%                     thisSavepath = newSavepath;
%                 end
%             end
            
            set(gca, 'FontName', obj.fontName) % numst be before x,y,zlabels
            
            % after title has been taken, do option to turn it off
            if ~printTitleFlag
                ax = gca;
                titlename = ax.Title.String;
                ax.Title.String = '';
            end

            printFnc = 'print';
            for ii = 1:numel(Extensions)
                thisExt = Extensions{ii};
                % if vector, set Painters renderer
                if ismember(thisExt,{'-dpdf','-deps','-depsc','-deps2','-depsc2',...
                        '-dmeta','-dsvg','-dps','-dpsc','-dps2','-dpsc2'})
                    Fh.Renderer = 'Painters';
                    printFnc = 'print';
                elseif ismember(thisExt,{'-djpeg','-dpng','-dtiff','-dtiffn','-dmeta',...
                        '-dbmpmono','-dbmp','-dbmp16m','-dbmp256','-dhdf','-dpbm',...
                        '-dpbmraw','-dpcxmono','-dpcx24b','-dpcx256','-dpcx16','-dpgm',...
                        '-dpgmraw','-dppm','-dppmraw'})
                    Fh.Renderer = 'OpenGL';
                    printFnc = 'print';
                elseif ismember(thisExt,{'fig'})
                    printFnc = 'savefig';
                elseif ismember(thisExt,{'svg'})
                    printFnc = 'fig2svg';
                end

                switch printFnc
                    case 'print'
                        print(['-r' num2str(obj.printDPI)],Fh,thisSavepath,thisExt);
                    case 'savefig'
                        savefig(Fh,thisSavepath,'compact')
                    case 'fig2svg'
                        % we need to convert the size to pixels to feed the
                        % size to fig2svg
                        fig2svg([thisSavepath '-i.svg'],Fh,[],[],0,...
                            [convertunit(obj.paperPosition(3), obj.paperUnits, 'pixels') ...
                            convertunit(obj.paperPosition(4), obj.paperUnits, 'pixels')],[],0)
                end
            end
            
            % in the end, we restore the title
            if ~printTitleFlag
                ax.Title.String = titlename;
            end

        end
        
    end
end

function rvalue = convertunit(value, from, to, parentheight)
% helper function from fig2svg
  % From SVG 1.1. Specification:
  % "1pt" equals "1.25px" (and therefore 1.25 user units)
  % "1pc" equals "15px" (and therefore 15 user units)
  % "1mm" would be "3.543307px" (3.543307 user units)
  % "1cm" equals "35.43307px" (and therefore 35.43307 user units)
  % "1in" equals "90px" (and therefore 90 user units)
  % Modification by Jonathon Harding:
  % MATLAB however, assumes a variable number of pixels per inch, and
  % assuming that the pixels match is dangerous.
  if nargin < 4
    parentheight = 1.25; % Default
  end
  
try
    ScreenPixelsPerInch = round(get(0, 'ScreenPixelsPerInch'));
catch
    ScreenPixelsPerInch = 96;
end
resolutionScaling = ScreenPixelsPerInch/96;

switch lower(from) % convert from input unit to points
    case 'pixels', rvalue = value*72/ScreenPixelsPerInch/resolutionScaling;
    case 'points', rvalue = value;
    case 'centimeters', rvalue = value/2.54*72;
    case 'inches', rvalue = value*72; % 72 points = 1 inch
    case 'normalized', rvalue = value*(parentheight*0.8);
    otherwise, error(['Unknown unit ', from, '.']);
end
switch lower(to) % convert from points to specified unit
    case 'pixels', rvalue = rvalue*ScreenPixelsPerInch/72;
    case 'points' % do nothing
    case 'centimeters', rvalue = rvalue*2.54/72;
    case 'inches', rvalue = rvalue/72; % 72 points = 1 inch
    case 'normalized', rvalue = value/(parentheight*0.8);
    otherwise, error(['Unknown unit ', to, '.']);
end

end

%             -dpdf
%             -deps
%             -depsc
%             -deps2
%             -depsc2
%             -dmeta
%             -dsvg
%             -dps
%             -dpsc
%             -dps2
%             -dpsc2
%             -djpeg
%             -dpng
%             -dtiff
%             -dtiffn
%             -dmeta
%             -dbmpmono
%             -dbmp
%             -dbmp16m
%             -dbmp256
%             -dhdf
%             -dpbm
%             -dpbmraw
%             -dpcxmono
%             -dpcx24b
%             -dpcx256
%             -dpcx16
%             -dpgm
%             -dpgmraw
%             -dppm
%             -dppmraw