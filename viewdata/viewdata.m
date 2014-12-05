function varargout = viewdata(varargin)
% VIEWDATA MATLAB code for viewdata.fig
%      VIEWDATA, by itself, creates a new VIEWDATA or raises the existing
%      singleton*.
%
%      H = VIEWDATA returns the handle to a new VIEWDATA or the handle to
%      the existing singleton*.
%
%      VIEWDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWDATA.M with the given input arguments.
%
%      VIEWDATA('Property','Value',...) creates a new VIEWDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewdata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewdata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewdata

% Last Modified by GUIDE v2.5 10-Jul-2014 12:51:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewdata_OpeningFcn, ...
                   'gui_OutputFcn',  @viewdata_OutputFcn, ...
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

function makeplot(handles, hObject)

% 0 = Don't display all catch contours
% 1 = Display all catch contours
displayCatch = 0;

% 0 = Don't display sonogram
% 1 = Display sonogram
displaySonogram = 0;

% 0 = Don't display contourmap
% 1 = Display contourmap
displayContourmap = 0;

% 0 = Don't display syllable bars
% 1 = Display syllable bars
displaySyllableBars = 0;

handles.data = handles.data(~cellfun('isempty',handles.data));
handles.latencies = handles.latencies(handles.latencies~=0);
handles.onsets = handles.onsets(handles.onsets~=0);

data = handles.data{handles.current};

set(handles.fname, 'String', data.filename);

% Huge block for actually plotting the data
if displayCatch % Display the catch contours
    for i=1:size(handles.catchdata,2)
        mydata = handles.catchdata(:,i)';
        mydata = mydata - mydata(round(data.b * length(mydata) / length(data.emg)));
        plot(linspace(data.stimx(1), data.stimx(end),numel(mydata)), mydata, 'b');
        hold on;
    end
end

if displaySonogram
    % Old: data.accessory.basefile.timerange
    imagesc([min(data.stimx), max(data.stimx)], data.accessory.basefile.freqrange - data.diff, abs(data.sonogram));
    set(gca, 'YDir', 'normal');
    hold on;
end

if displayContourmap
    imagesc([min(data.stimx), max(data.stimx)], data.accessory.basefile.freqrange - data.diff, abs(data.contourmap));
    set(gca, 'YDir', 'normal');
    hold on;
end

plot(data.stimx, data.stim, 'r');
hold on;
plot(data.stimx, data.confint, 'g');
hold on;
plot([data.intersections data.intersections], [min(data.stim) max(data.stim)], 'b');
hold on;
plot([data.stimonset data.stimonset], [min(data.stim) max(data.stim)], 'k');
hold off;

% disp(['Stim onset ', 'Intersection'])
% disp([data.stimonset, data.intersections]);

xlabel('time [ms]');
ylabel('frequency [kHz]');

% Choose default command line output for viewdata
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function [handles] = cca(directory,handles)

% 0 = Don't use spine interpolation
% 1 = Use spline interpolation
makeSplines = 1;

filename = 'ZZZCOMPILED.mat';

f.all_dirs{1}.name = 'contours';
if ~exist('directory','var') || isempty(directory)
    mydir = dir;
else
    mydir = dir(directory);
end
counter = 0;
for i=1:numel(mydir)
    if ~mydir(i).isdir && ~strcmp(mydir(i).name, filename)
        counter = counter + 1;
        files{counter} = mydir(i).name;
    end
end
myfile = load(files{1});
f.all_dirs{1}.files = files;
if ~exist(filename,'file')
    [contours sonos accessory ids] = compilecontours(f);
    save(filename,'contours','sonos','accessory','ids');
else
    load(filename);
end

% Time and freqency ranges for actual data
timescale = accessory.basefile.timerange(2) - accessory.basefile.timerange(1);
freqscale = accessory.basefile.freqrange;

contours = cellfun(@(x) x*(freqscale(2)-freqscale(1))/size(myfile.contourmap,1)+freqscale(1),contours,'UniformOutput',false);
ids = logical(ids);

catchdata = cell2mat(contours(~ids)');
stimdata = cell2mat(contours(ids)');
timeranges = accessory.timerange(ids);
emgdata = cell2mat(accessory.emg(ids)');
sonograms = accessory.sonogram(ids);
contourmaps = accessory.contourmap(ids);

stimfiles = accessory.filenames(ids);

conf_int = zeros(size(catchdata,1),1); % Main confidence interval
conf_int2 = zeros(size(catchdata,1),1);

% Calculate confidence interval using 95th quantile
for i=1:size(catchdata,1)
    % Remove zero-valued elements
    mydata = catchdata(i,:)';
    mydata = mydata(mydata~=freqscale(1));
    if numel(mydata) < 5
        mydata = catchdata(i,:)';
    end
    pd = fitdist(mydata,'Normal');
    % Which one of these is better?
    conf_int2(i) = pd.mu + 2*pd.std; % Mean plus 2 standard deviations
    conf_int(i) = quantile(mydata,0.975); % 95th quantile
end

splineres = 0:0.01:numel(conf_int); % Resolution of splines

if makeSplines
    conf_int = spline(0:numel(conf_int)-1,conf_int,splineres);
end

% Figure out where stim data becomes greater than conf_int
for i=1:size(stimdata,2)
    invalid = false; % Marker for invalid contours
    a = emgdata(:,i);
    a = find(a>mean(a));
    b = a(1); % Onset of stimulation
    a = b; % First stimulation point in emg
    mystim = stimdata(:,i);
    if makeSplines
        mystim = spline(0:numel(mystim)-1,mystim,splineres);
    end
    
    scalar = (length(mystim)/length(emgdata(:,i)));

    my_conf_int = zeros(size(catchdata,1),1); % Main confidence interval
    pt = round(a * size(catchdata,1) / length(emgdata(:,i)));
    if pt > size(catchdata, 1) || pt <= 0
        continue;
    end
    
    % Calculate confidence interval using 95th quantile
    for k=1:size(catchdata,2)
        catchdata(:,k) = catchdata(:,k) - catchdata(pt,k);
    end
    diff = mystim(round(a * length(mystim) / length(emgdata(:,i))));
    mystim = mystim - diff;
    for k=1:size(catchdata,1)
        % Remove zero-valued elements
        mydata = catchdata(k,:)';
        if numel(mydata) < 5
            mydata = catchdata(k,:)';
        end
        my_conf_int(k) = quantile(mydata,0.975); % 95th quantile
    end
    
    my_conf_int = spline(0:numel(my_conf_int)-1,my_conf_int,splineres);
    
%     if mystim(round(a*scalar)) >= my_conf_int(round(a*scalar))
%         invalid = true;
%     end
    
    while invalid == false && mystim(round(a*scalar)) <= my_conf_int(round(a*scalar))
        % Work forwards to find where the line crosses
        a = a + 1/scalar;
        if round(a*scalar) > numel(my_conf_int)
            invalid = true;
            break;
        end
    end
    
    slope = zeros(1,numel(mystim)-2);
    for j=6:numel(mystim)-5
        diff_a = my_conf_int(j-5) - mystim(j-5);
        diff_b = my_conf_int(j+5) - mystim(j+5);
        slope(j) = diff_a - diff_b;
    end
    
    if invalid
        continue;
    end
    
    timerange = timeranges{i};
    timescale = timerange(2) - timerange(1);
    
    data{i}.intersections = a * timescale / length(emgdata(:,i)) + timerange(1);
    data{i}.latencies = (a - b) * timescale / length(emgdata(:,i)) + timerange(1);
    
    % Store variables in struct
    data{i}.diff = diff;
    data{i}.stimx = linspace(timerange(1),timerange(2),length(mystim));
    data{i}.emgx = linspace(timerange(1),timerange(2),length(emgdata(:,i)));
    data{i}.stim = mystim;
    data{i}.confint = my_conf_int;
    data{i}.emg = emgdata(:,i);
    data{i}.stimonset = b * timescale / length(emgdata(:,i)) + timerange(1);
    data{i}.slope = slope;
    data{i}.sonogram = sonograms{i};
    data{i}.contourmap = contourmaps{i};
    data{i}.accessory = accessory;
    data{i}.filename = stimfiles{i};
    data{i}.b = b;
end

handles.catchdata = catchdata;

data = data(~cellfun('isempty',data)); % Remove empty cells

latencies = zeros(numel(data),1);
latencies2 = zeros(numel(data),1);
onsets = zeros(numel(data),1);
for i = 1:numel(data)
    latencies(i) = data{i}.latencies;
    onsets(i) = data{i}.intersections - data{i}.latencies;
end

handles.data = data;
handles.latencies = latencies;
handles.latencies2 = latencies2;
handles.onsets = onsets;

% --- Executes just before viewdata is made visible.
function viewdata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewdata (see VARARGIN)

[handles] = cca(pwd,handles);

handles.current = 1; % current point in data

makeplot(handles,hObject);

% Choose default command line output for viewdata
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes viewdata wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = viewdata_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
foldername = uigetdir;
if foldername
    handles = cca(foldername,handles);
    handles.current = 1; % current point in data

    makeplot(handles,hObject);

    % Choose default command line output for viewdata
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, filterindex] = uiputfile({'*.mat','MAT files (*.mat)'});
if filename
    data = handles.data;
    latencies = handles.latencies;
    onsets = handles.onsets;
    save(fullfile(pathname,filename),'data', 'latencies', 'onsets');
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes during object creation, after setting all properties.
function allplots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on button press in deletebutton.
function deletebutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data{handles.current} = [];
handles.latencies(handles.current) = 0;
handles.onsets(handles.current) = 0;
if handles.current == length(handles.data)
    handles.current = handles.current - 1;
end
makeplot(handles, hObject);

% --- Executes on button press in prevbutton.
function prevbutton_Callback(hObject, eventdata, handles)
% hObject    handle to prevbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current > 1
    handles.current = handles.current - 1;
    makeplot(handles, hObject);
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current < length(handles.data)
    handles.current = handles.current + 1;
    makeplot(handles, hObject);
end


% --------------------------------------------------------------------
function AboutButton_Callback(hObject, eventdata, handles)
% hObject    handle to AboutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = msgbox('This program allows users to view and delete contours that may have some problem with them. It is intended as an accessory for analyzecontours.','About')
