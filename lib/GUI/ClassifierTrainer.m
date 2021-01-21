function varargout = ClassifierTrainer(varargin)
% CLASSIFIERTRAINER MATLAB code for ClassifierTrainer.fig
%      CLASSIFIERTRAINER, by itself, creates a new CLASSIFIERTRAINER or raises the existing
%      singleton*.
%
%      H = CLASSIFIERTRAINER returns the handle to a new CLASSIFIERTRAINER or the handle to
%      the existing singleton*.
%
%      CLASSIFIERTRAINER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFIERTRAINER.M with the given input arguments.
%
%      CLASSIFIERTRAINER('Property','Value',...) creates a new CLASSIFIERTRAINER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClassifierTrainer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClassifierTrainer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClassifierTrainer

% Last Modified by GUIDE v2.5 17-Feb-2020 12:05:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClassifierTrainer_OpeningFcn, ...
                   'gui_OutputFcn',  @ClassifierTrainer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ClassifierTrainer is made visible.
function ClassifierTrainer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClassifierTrainer (see VARARGIN)

% Choose default command line output for ClassifierTrainer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

Track = varargin{2}; % training data
popNr = varargin{3};
parent = varargin{4}; % handles parent

if numel(Track) ~= 1
    error('Number of parsed tracks should be one.')
end

handles.pushbutton_loadclassifier.UserData = Track;
% handles.handles.pushbutton_nosegmentation.UserData = Track; % output
handles.uipanel_general.UserData = popNr;
handles.pushbutton_saveclassifier.UserData = parent;

% UIWAIT makes ClassifierTrainer wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ClassifierTrainer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.pushbutton_nosegmentation.UserData;
% The figure can be deleted now
delete(handles.figure1);

% --- Helper functions

function errorBox(handles,err)
% if there is an error invoked, we output it in a dialog box
if handles.pushbutton_saveclassifier.UserData.checkbox_debug.Value % get from parent
    errordlg(getReport(err,'extended','hyperlinks','off'),'Error');
else
    errordlg(getReport(err,'basic','hyperlinks','off'),'Error');
end

function updateCategories(handles)

TrackC = handles.pushbutton_newclassifier.UserData;

numCat = numel(TrackC);
catNames = {TrackC.name};

handles.listbox_categories.String = catNames; % set names of categories

for ii = 1:4
    thisName = ['pushbutton_category' num2str(ii)];
    thisName2 = ['text_category' num2str(ii)];
    if ii > numCat
        handles.(thisName).Visible = 'off';
        handles.(thisName2).Visible = 'off';
    else
        handles.(thisName).Visible = 'on';
        handles.(thisName2).Visible = 'on';
        handles.(thisName).String = catNames{ii};
    end

end

function plotClassification(handles)

idx = handles.pushbutton_startclassification.UserData;
nManualClassifications = handles.edit_mintrackpoints.UserData;

cla(handles.axes1);

% if we're done with all classifications, stop
if idx > nManualClassifications
    finishedClassification(handles)
    return
end

% else

% set text
txt = [num2str(idx) '/' num2str(nManualClassifications)];
handles.text_classification.String = {'Classification in progress:',txt};

% plot

Track = handles.pushbutton_loadclassifier.UserData;
sampleIDX = handles.text_mintrackpoints.UserData;
index = sampleIDX(idx);
[unitFactor,unitLabel] = Track.getUnitFactor('pixelsize');
track = Track.coords{index};

x = track(:,1).*unitFactor;
y = track(:,2).*unitFactor;

plot(handles.axes1, x, y, '-x', ...
    'LineWidth',1.5);



MBCC = Track.MinBoundCircleCenter(sampleIDX(idx),:).*unitFactor;
MBCR = Track.MinBoundCircleRadius(sampleIDX(idx)).*unitFactor;
CoM = Track.CenterOfMass(sampleIDX(idx),:).*unitFactor;

plotrange_m = str2double(handles.edit_plotrange.String); % use userspecified range
if ~handles.checkbox_scaleauto.Value && ~isnan(plotrange_m)
    addLim = plotrange_m/2;
else % do auto based on MBCR
    addLim = MBCR*1.1;
end
xln = [MBCC(1)-addLim MBCC(1)+addLim];
yln = [MBCC(2)-addLim MBCC(2)+addLim];
xlim(xln);
ylim(yln);

hold on;
plot(handles.axes1,MBCC(1),MBCC(2),'ro');
plot(handles.axes1,CoM(1),CoM(2),'go');
hold off;

xlabel(['x (' unitLabel ')'])
ylabel(['y (' unitLabel ')'])
axis equal

function setClassification(handles,nCat)

idx = handles.pushbutton_startclassification.UserData;
handles.text_classification.UserData(idx) = nCat;

idx = idx + 1;
handles.pushbutton_startclassification.UserData = idx;



plotClassification(handles)

function finishedClassification(handles)

handles.pushbutton_startclassification.UserData = [];

cla(handles.axes1)

TrackC = handles.pushbutton_newclassifier.UserData;
Track = handles.pushbutton_loadclassifier.UserData;
sampleIDX = handles.text_mintrackpoints.UserData;
nManualClassifications = handles.edit_mintrackpoints.UserData;
sampleIDX = sampleIDX(1:nManualClassifications);
categorization = handles.text_classification.UserData;

Track = Track.resetUnitSystem; % set back to SI

% txt = cell(numel(TrackC),1);
for ii = 1:numel(TrackC)
    take = categorization == ii;
    coords = Track.coords(sampleIDX(take));
    time = Track.time(sampleIDX(take));
    % --- convert to SI
    uFc = Track.getUnitFactor('pixelsize');
    uFt = Track.getUnitFactor('dt');
    coords = cellfun(@(x) x.*uFc,coords,'UniformOutput',false);
    time = cellfun(@(x) x.*uFt,time,'UniformOutput',false);
    % --- add
    TrackC(ii) = TrackC(ii).addTrack(time,coords);
end

% have to recompute properties, due to addition new tracks.
wb = waitbar(0,'Computing center','Name','Computing all track properties');
TrackC = TrackC.computeCenter;
waitbar(0.2,wb,'Computing entropy')
TrackC = TrackC.computeEntropy;
waitbar(0.4,wb,'Computing PCA')
TrackC = TrackC.computePCA;
waitbar(0.6,wb,'Computing MSD')
TrackC = TrackC.computeMSD;
% waitbar(0.8,wb,'Computing D')
% TrackC = TrackC.setArrayProperty('Dest',Dest.MSD);
% TrackC = TrackC.computeDiffusion;
close(wb)

handles.pushbutton_newclassifier.UserData = TrackC;

updateText(handles,'Finished classification')

pushbutton_plotTrainingset_Callback([], [], handles) % plot track


function updateText(handles,msg)

TrackC = handles.pushbutton_newclassifier.UserData;
    
nTracks = [TrackC.nTracks];

txt = cell(numel(TrackC),1);
for ii = 1:numel(TrackC)
    txt{ii} = [TrackC(ii).name ': ' num2str(nTracks(ii))];
end

txt2 = sprintf('Trainings set contains %i tracks.',sum([TrackC.nTracks]));

handles.text_classification.String = [{msg},...
    {txt2},txt{:}];

% --- Executes on selection change in listbox_categories.
function listbox_categories_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_categories (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_categories contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_categories


% --- Executes during object creation, after setting all properties.
function listbox_categories_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_categories (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_addcategory.
function pushbutton_addcategory_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addcategory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    TrackC = handles.pushbutton_newclassifier.UserData;
    
    if ~isempty(TrackC)
        numCat = numel(TrackC);

        if numCat >= 4
            error('Maximum 4 categories are allowed.')
        end

        TrackC = [TrackC tracks()]; % units pixelsize and dt are initialised to SI by default
        TrackC(end).name = ['Category ' num2str(numCat+1)];

        handles.pushbutton_newclassifier.UserData = TrackC;
        
        updateCategories(handles)
        updateText(handles,'Added category')
    end

catch err
    errorBox(handles,err)
end




% --- Executes on button press in pushbutton_deletecategory.
function pushbutton_deletecategory_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletecategory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    TrackC = handles.pushbutton_newclassifier.UserData;
    
    if ~isempty(TrackC)

        numCat = numel(TrackC);

        if numCat == 1
            error('Minimum of 1 categories is required.')
        end

        thisCat = handles.listbox_categories.Value;
        TrackC(thisCat) = [];
        handles.listbox_categories.Value = 1;

        handles.pushbutton_newclassifier.UserData = TrackC;
        
        updateCategories(handles)
        updateText(handles,'Deleted category')
        
    end

catch err
    errorBox(handles,err)
end


% --- Executes on button press in pushbutton_newclassifier.
function pushbutton_newclassifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_newclassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% The classifier is a tracks object array, in which each population is a
% classified category. In this way, the classifier can be easily plotted
% and trainings data can be added. However, to allow for extension, the
% units should be the same, hence, the classifier is always defined in SI
% units.

try
    TrackC = tracks(); % units pixelsize and dt are initialised to SI by default. Thus, first convert added tracks to SI units.
    TrackC.name = 'Category 1';
    TrackC.dt = 1;
    TrackC.dte = 0;

    handles.pushbutton_newclassifier.UserData = TrackC;

    updateCategories(handles)
    updateText(handles,'Started new classifier')
    
catch err
    errorBox(handles,err)
end


% --- Executes on button press in pushbutton_loadclassifier.
function pushbutton_loadclassifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadclassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % --- load Track
    [fname,pname,idx] =uigetfile( ...
        {'*.mat','MATLAB Files (*.mat)'},...
       'Please select the file containing the track classifier');
   
%     u = symunit;
    switch idx
        case 0 % cancel
            return
        case 1
            fprintf('Loading %s\n', fullfile(pname,fname));  
            in = load(fullfile(pname,fname)); % Loads 'Track'
            fn = fieldnames(in);
            nTracksObj = 0;
            for ii = 1:numel(fn)
                if isa(in.(fn{ii}),'tracks')
                    TrackC = in.(fn{ii});
                    nTracksObj = nTracksObj + 1;
                end
            end
            if nTracksObj ~= 1
                error('Expected one ''tracks'' object but found %d instead.',...
                    nTracksObj);
            end
            handles.pushbutton_newclassifier.UserData = TrackC;
    end
    
    updateCategories(handles)
    updateText(handles,'Loaded classifier')
        
catch err
    errorBox(handles,err)
end
        

% --- Executes on button press in pushbutton_startclassification.
function pushbutton_startclassification_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_startclassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    TrackC = handles.pushbutton_newclassifier.UserData;
    if ~isempty(TrackC)
        Track = handles.pushbutton_loadclassifier.UserData;
        nManualClassifications = str2double(get(handles.edit_noclassify,'String'));
        if ~isnan(nManualClassifications) 
            doneIDX = handles.text_mintrackpoints.UserData;
            
            tooShortIDX = Track.nTrackPoints < str2double(get(handles.edit_mintrackpoints,'String'));
            tooShortIDX(doneIDX) = true;
            goodIDX = find(~tooShortIDX);
            
            if numel(goodIDX) >= nManualClassifications
                sampleIDX = randsample(goodIDX,nManualClassifications);
                handles.text_mintrackpoints.UserData = vertcat(sampleIDX,doneIDX); % store the done idx without classifying them, because numel(sampleIDX) == nClassifications
                handles.pushbutton_startclassification.UserData = 1; % start at one
                handles.edit_mintrackpoints.UserData = nManualClassifications;
                handles.text_classification.UserData = nan(nManualClassifications,1);
                plotClassification(handles);
            else
                error('The number of remaining tracks (%i) is smaller than the training set (%i).',numel(goodIDX),nManualClassifications);
            end
        end
    end
catch err
    errorBox(handles,err)
end



function edit_noclassify_Callback(hObject, eventdata, handles)
% hObject    handle to edit_noclassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_noclassify as text
%        str2double(get(hObject,'String')) returns contents of edit_noclassify as a double


% --- Executes during object creation, after setting all properties.
function edit_noclassify_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_noclassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mintrackpoints_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mintrackpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mintrackpoints as text
%        str2double(get(hObject,'String')) returns contents of edit_mintrackpoints as a double


% --- Executes during object creation, after setting all properties.
function edit_mintrackpoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mintrackpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_category1.
function pushbutton_category1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_category1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    idx = handles.pushbutton_startclassification.UserData;
    if isempty(idx)
        return
    end

    setClassification(handles,1)
catch err
    errorBox(handles,err)
end



% --- Executes on button press in pushbutton_category2.
function pushbutton_category2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_category2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    idx = handles.pushbutton_startclassification.UserData;
    numCat = numel(handles.pushbutton_newclassifier.UserData);
    if isempty(idx) || numCat < 2
        return
    end

    setClassification(handles,2)
catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_category3.
function pushbutton_category3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_category3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    idx = handles.pushbutton_startclassification.UserData;
    numCat = numel(handles.pushbutton_newclassifier.UserData);
    if isempty(idx) || numCat < 3
        return
    end

    setClassification(handles,3)
catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_category4.
function pushbutton_category4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_category4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    idx = handles.pushbutton_startclassification.UserData;
    numCat = numel(handles.pushbutton_newclassifier.UserData);
    if isempty(idx) || numCat < 4
        return
    end

    setClassification(handles,4)
catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_saveclassifier.
function pushbutton_saveclassifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveclassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    TrackC = handles.pushbutton_newclassifier.UserData;
    if ~isempty(TrackC)
        [fname,pname] = uiputfile('.mat','Save tracks classifier');
        wb = waitbar(0.5,'Saving classifier...','Name','Please wait');
        save(fullfile(pname,fname),'TrackC','-v7.3');
        close(wb)
    end
catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_renamecategory.
function pushbutton_renamecategory_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_renamecategory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    
    TrackC = handles.pushbutton_newclassifier.UserData;
    
    if ~isempty(TrackC)

        thisCat = handles.listbox_categories.Value;

        prompt = {'Enter category name:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {TrackC(thisCat).name};
        newName = inputdlg(prompt,dlgtitle,dims,definput);
        TrackC(thisCat).name = newName{1};

        handles.pushbutton_newclassifier.UserData = TrackC;

        updateCategories(handles)
        updateText(handles,'Renamed category')
    
    end

catch err
    errorBox(handles,err)
end


% --- Executes on button press in pushbutton_plotTrainingset.
function pushbutton_plotTrainingset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotTrainingset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    
    TrackC = handles.pushbutton_newclassifier.UserData;
    
    if ~isempty(TrackC)
        
        if any([TrackC.nTracks] == 0)
            
            nTracks = [TrackC.nTracks];
            idx = nTracks == 0; % cannot plot empty tracks
            
            TrackC(idx) = [];
            
        end
            
        figure
        hps = TrackC.plotTracks;

        % manually set legend
        nCat = numel(TrackC);
        legendTxt = cell(nCat,1);
        legendN = cumsum([TrackC.nTracks]);
        legendHps = zeros(nCat,1);
        for ii = 1:nCat
            legendTxt{ii} = TrackC(ii).name;
            legendHps(ii) = hps(legendN(ii));
        end
        legend(legendHps,legendTxt{:})
    
    end

catch err
    errorBox(handles,err)
end


% --- Executes on button press in pushbutton_cancelclassification.
function pushbutton_cancelclassification_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancelclassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
        
    classifactionIDX = handles.pushbutton_startclassification.UserData;
    
    if ~isempty(classifactionIDX)
        
        handles.pushbutton_startclassification.UserData = []; % only variable checked to mark progress classification
    
        handles.edit_mintrackpoints.UserData = []; % is not checked, doesnt have to be reset
        handles.text_classification.UserData = []; % is not checked, doesnt have to be reset
        
        updateText(handles,'Cancelled classification')
        
    end

catch err
    errorBox(handles,err)
end


% --- Executes on button press in pushbutton_comptree.
function pushbutton_comptree_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_comptree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    
    TrackC = handles.pushbutton_newclassifier.UserData;
    
    if ~isempty(TrackC)
        maxnumsplits = str2double(handles.edit_maxnumsplits.String);
        if ~isnan(maxnumsplits)
            
            wb = waitbar(0,'Preparing data','Name','Computing all decision tree');

            TrackProps_wanted = handles.pushbutton_saveclassifier.UserData.popupmenu_TrackProp1.UserData; % from diffusionlab
            % name, location in obj.(location), unit for unitFactor, column of data

            Track = handles.pushbutton_loadclassifier.UserData;

            take = zeros(size(TrackProps_wanted,1),1);
            for ii = 1:size(TrackProps_wanted,1)
                try
                    prop = Track.getProp(TrackProps_wanted{ii,2},true);
                    if ~isempty(prop) && size(prop,2) >= TrackProps_wanted{ii,4}
                        take(ii) = 1;
                    end
                catch
                    % do nothing
                end
            end
            
            TrackProps = TrackProps_wanted(find(take),:);
            TrackProps = TrackProps(~strncmpi(TrackProps(:,2),'Dest',4),:); % throw out diffusion constant
%             TrackProps = TrackProps(~strcmp(TrackProps(:,1),'Localization variance'),:); % throw out diffusion constant
%             TrackProps = TrackProps(~strcmp(TrackProps(:,1),'Diffusion SNR'),:); % throw out diffusion constant
            
%             if isempty(TrackProps)
%                 close(wb)
%                 error('No track properties computed for tracks.')
%             end
                

            nTracks = [TrackC.nTracks];
            nProps = size(TrackProps,1);


            X = nan(sum(nTracks),nProps);
            
            for ii = 1:size(TrackProps,1) % loop over categories
                data = vertcat(TrackC.(TrackProps{ii,2})); % check extraction from cell array
                data = data(:,TrackProps{ii,4});
                X(:,ii) = data;
            end
            
            Y = cell(sum(nTracks),1);
            c = 1;
            for ii = 1:numel(TrackC)

                Y(c:c+nTracks(ii)-1) = {TrackC(ii).name};
                c = c + nTracks(ii);

            end
            
            waitbar(0.5,wb,'Computing decision tree')

            ctree = fitctree(X,Y,'MaxNumSplits',maxnumsplits,'PredictorNames',TrackProps(:,1)); % create classification tree
            
            view(ctree,'mode','graph');

            handles.pushbutton_comptree.UserData = ctree;
            handles.pushbutton_segment.UserData = TrackProps;
            
            close(wb)
            
        end
        
    end

catch err
    errorBox(handles,err)
end



function edit_maxnumsplits_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxnumsplits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxnumsplits as text
%        str2double(get(hObject,'String')) returns contents of edit_maxnumsplits as a double


% --- Executes during object creation, after setting all properties.
function edit_maxnumsplits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxnumsplits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plottree.
function pushbutton_plottree_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plottree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    ctree = handles.pushbutton_comptree.UserData;
    
    if ~isempty(ctree)
        view(ctree,'mode','graph');
    end

catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_segment.
function pushbutton_segment_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    ctree = handles.pushbutton_comptree.UserData;
    if isempty(ctree)
        error('Please compute decision tree first.')
    end
    
    TrackProps = handles.pushbutton_segment.UserData;
    TrackC = handles.pushbutton_newclassifier.UserData;
    Track = handles.pushbutton_loadclassifier.UserData;

    X = nan(Track.nTracks,size(TrackProps,1));

    for ii = 1:size(TrackProps,1) % loop over categories
        data = vertcat(Track.(TrackProps{ii,2})); % check extraction from cell array
        data = data(:,TrackProps{ii,4});
        data = data.*Track.getUnitFactor(TrackProps{ii,3}); % data is Track is not in SI units! Only our classifier is!
        X(:,ii) = data;
    end
    
    TrackClass = predict(ctree,X);
    
    uniqueClass = unique(TrackClass);
    nClass = numel(uniqueClass);
    
    Track_out = tracks.empty(nClass,0);
    
    for ii = 1:nClass
%         take = strcmp(TrackClass,uniqueClass(ii));
        Track_out(ii,1) = Track.segment(strcmp(TrackClass,uniqueClass{ii}));
    end
    
%     classNames = {TrackC.name};
    popNr = handles.uipanel_general.UserData;
    
    Track_out(1).segTree = Track_out(1).segTree.add(popNr,uniqueClass{:});
    Track_out = Track_out.setArrayProperty('segTree',Track_out(1).segTree);
    
    handles.pushbutton_nosegmentation.UserData = Track_out;
    
    figure1_CloseRequestFcn(handles.figure1, [], handles)
    
    
catch err
    errorBox(handles,err)
end

% --- Executes on button press in pushbutton_nosegmentation.
function pushbutton_nosegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nosegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


try
    Track = handles.pushbutton_loadclassifier.UserData; % give original data back
    handles.pushbutton_nosegmentation.UserData = Track;
    handles.pushbutton_nosegmentation.UserData = [];
    figure1_CloseRequestFcn(handles.figure1, [], handles)
catch err
    errorBox(handles,err)
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end

catch err
    errorBox(handles,err)
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Character
    case 'z'
        pushbutton_category1_Callback([], [], handles)
    case 'x'
        pushbutton_category2_Callback([], [], handles)
    case 'c'
        pushbutton_category3_Callback([], [], handles)
    case 'v'
        pushbutton_category4_Callback([], [], handles)
    otherwise
        
end


% --- Executes on button press in pushbutton_back.
function pushbutton_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    idx = handles.pushbutton_startclassification.UserData;
    if isempty(idx)
        return
    end
    
    if handles.pushbutton_startclassification.UserData > 1
        handles.pushbutton_startclassification.UserData = handles.pushbutton_startclassification.UserData - 1; % go back
    end

    plotClassification(handles)
catch err
    errorBox(handles,err)
end


% --- Executes on button press in checkbox_scaleauto.
function checkbox_scaleauto_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_scaleauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_scaleauto

plotClassification(handles)

function edit_plotrange_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plotrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plotrange as text
%        str2double(get(hObject,'String')) returns contents of edit_plotrange as a double

plotClassification(handles)

% --- Executes during object creation, after setting all properties.
function edit_plotrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plotrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_cancelclassification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancelclassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_addall.
function pushbutton_addall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.text_classification.UserData(:) = handles.listbox_categories.Value; % set all to selected category
finishedClassification(handles)
