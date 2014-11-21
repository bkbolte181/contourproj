function varargout = analyzecontours(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzecontours_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzecontours_OutputFcn, ...
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

function finalHandles = loadFiles(handles, hObject)
handles.filename = handles.all_files(handles.current_file).name;
f = fullfile(handles.pathname, strcat(handles.filename, '.not.mat')); % Path to song metadata

if exist(f, 'file')
    handles.f = load(f);
end

for i=1:length(handles.all_files)
    files_list{i} = handles.all_files(i).name;
end
set(handles.files_menu,'String',files_list);
set(handles.files_menu,'Value',handles.current_file);

guidata(hObject, handles);

finalHandles = handles; % Return handles

function generate_spectrogram(handles, hObject, contourPoint)
channel = get(handles.channel_val,'String');
if strcmp(handles.filename(end-3:end),'cbin') % Handle cbin files
    [dat fs] = evsoundin(handles.pathname, handles.filename, channel); % look into modifying this
elseif strcmp(handles.filename(end-2:end),'mat') % Handle mat files
    file = load(fullfile(handles.pathname, handles.filename));
    dat = file.data;
    fs = file.fs;
    if isfield(handles,'f')
        handles = rmfield(handles,'f');
    end
end

if ~exist('dat','var') || isempty(dat) || numel(dat) == 0
    error('generate_spectrogram: Failed to load file');
end

buffer = 30; % Buffer before and after syllable - change this if desired, must be integer
handles.buffer = buffer;
handles.fs = fs;
handles.polynomialmode = false;

if length(get(handles.edit_syllable,'String')) > 0
    if handles.current_syllable == 0
        handles.current_syllable = 1;
    end
else
    handles.current_syllable = 0;
    handles.in = 0;
    handles.out = 0;
end

if ~isfield(handles,'f');
    handles.current_syllable = 0;
elseif handles.current_syllable > 0
    syllable = get(handles.edit_syllable, 'String');
    counter = 0;
    for i=(1:numel(handles.f.onsets))
        if sum(strfind(syllable,handles.f.labels(i)))
            counter = counter + 1;
            if counter == handles.current_syllable
                handles.in = handles.f.onsets(i);
                handles.out = handles.f.offsets(i);
                in = round((handles.in-buffer)*fs/1000); % start point
                out = round((handles.out+buffer)*fs/1000); % end point
                dat = dat(max(in-buffer, 1):min(out+buffer, numel(dat))); % one syllable
            end
            % emgdata = emg(max(in, 1):min(out, numel(emg))); % associated emg data
        end
    end
    handles.max_counter = counter;
end

% Get variables from edit boxes
N = str2double(get(handles.edit_nfft, 'String'));
OVERLAP = str2double(get(handles.edit_overlap, 'String'));
Threshold = str2array(get(handles.edit_arthreshold, 'String'));
if length(Threshold) ~= 3
    Threshold = Threshold(1);
end
sigmas = str2array(get(handles.edit_timescales, 'String'));
min_freq = str2double(get(handles.edit_minfreq, 'String'));
max_freq = str2double(get(handles.edit_maxfreq, 'String'));
angles = str2array(get(handles.edit_angles, 'String'))*pi;

new = {N, OVERLAP, Threshold, sigmas, angles};

if exist('contourPoint','var') && ~any(any(contourPoint))
    reset = 'yes';
end

if ~isfield(handles,'N') || ~isequal(dat,handles.dat) || ~any(any(handles.contourmap)) || exist('reset','var') || N ~= handles.N || OVERLAP ~= handles.OVERLAP || ~isequal(Threshold,handles.Threshold) || ~isequal(sigmas,handles.sigmas) || min_freq ~= handles.min_freq || max_freq ~= handles.max_freq || ~isequal(angles,handles.angles)
    [handles.contourmap, handles.sonogram, f, t] = zf_contour_func(dat, fs, N, OVERLAP, Threshold, sigmas, [min_freq max_freq], angles);
    [handles.N, handles.OVERLAP, handles.Threshold, handles.sigmas, handles.angles] = deal(new{:});
    handles.times = t;
    handles.freqs = f;
end

[a,b] = size(handles.contourmap);
timescale = 1000*b*(N-OVERLAP)/fs;
if isfield(handles,'in') && isfield(handles,'out') && handles.out ~= 0
    min_time = handles.in-buffer;
    max_time = handles.out+buffer;
    conv = b/(max_time - min_time); % Ratio from ms to pix
    inpoint = round(buffer * conv);
    outpoint = b - inpoint;
    handles.mapbounds = [inpoint outpoint];
    % Draw vertical lines around contour
    if get(handles.syllable_bars, 'Value');
        handles.contourmap(:,inpoint) = max(max(handles.contourmap));
        handles.sonogram(:,inpoint) = max(max(handles.sonogram));
        handles.contourmap(:,outpoint) = max(max(handles.contourmap));
        handles.sonogram(:,outpoint) = max(max(handles.sonogram));
    end
else
    min_time = handles.times(1);
    max_time = handles.times(end);
    conv = b/(max_time - min_time); % Ratio from ms to pix
    inpoint = round(min_time * conv);
    outpoint = b - inpoint;
    handles.mapbounds = [inpoint outpoint];
end
handles.timescale = [min_time max_time];
vibrancy = get(handles.vibrancy_slider, 'Value');
FONTSIZE = 12;

vibrancy = 20^(vibrancy-1);

if isfield(handles,'dat') && isequal(dat,handles.dat) && isequal(handles.min_freq, min_freq) && isequal(handles.max_freq, max_freq)
    v = axis;
end

% timescale=1000*(b*Nshift)/SAMPLING;

% Code allowing you to click on a contour
if exist('contourPoint','var') && ~isempty(contourPoint) && any(any(contourPoint))
    [a b] = size(handles.contourmap);
    contourPoint(:,1) = round(contourPoint(:,1) * a / (max_freq - min_freq));
    contourPoint(:,2) = round((contourPoint(:,2) - min_time) * b / (max_time - min_time));
    try
        xy = sub2ind(size(handles.contourmap),contourPoint(:,1),contourPoint(:,2));
    catch
        myerror = strcat({'Invalid subscripts for sub2ind => Contourmap size: '}, mat2str(size(handles.contourmap)), {' Coordinates given: '}, mat2str([contourPoint(:,1) contourPoint(:,2)]));
        myerror = myerror{1};
        error(myerror);
    end
    CC = bwconncomp(handles.contourmap);
    for i = 1:numel(CC.PixelIdxList)
        if ~sum(ismember(xy,CC.PixelIdxList{i}))
            handles.contourmap(CC.PixelIdxList{i}) = 0;
        end
    end
end

% Draws the contour image
if get(handles.toggle_contour, 'Value'); % Generate contours
    imagesc([min_time max_time],[min_freq max_freq],log(handles.contourmap+vibrancy));
else % Generate regular spectrogram
    imagesc([min_time max_time],[min_freq max_freq], log(abs(handles.sonogram)+vibrancy));
end
set(gca,'YDir','normal');
xlabel('time [ms]', 'FontSize', FONTSIZE);
ylabel('frequency [kHz]', 'FontSize', FONTSIZE);

colormap('hot');

if isfield(handles,'dat') && isequal(dat,handles.dat) && isequal(handles.min_freq, min_freq) && isequal(handles.max_freq, max_freq)
    axis(v);
else
    handles.min_freq = min_freq;
    handles.max_freq = max_freq;
    handles.dat = dat;
end

% Choose default command line output for analyzecontours
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function f = getfields(handles)
f.times = handles.times;
f.freqs = handles.freqs;
f.songdata = handles.dat;
f.data = getimage(gcf);
f.contourmap = handles.contourmap;
f.sonogram = handles.sonogram;
f.mapbounds = handles.mapbounds;
f.N = str2double(get(handles.edit_nfft, 'String'));
f.OVERLAP = str2double(get(handles.edit_overlap, 'String'));
f.Threshold = str2array(get(handles.edit_arthreshold, 'String'));
if length(f.Threshold) ~= 3
    f.Threshold = f.Threshold(1);
end
f.sigmas = str2array(get(handles.edit_timescales, 'String'));
min_freq = str2double(get(handles.edit_minfreq, 'String'));
max_freq = str2double(get(handles.edit_maxfreq, 'String'));
f.freqrange = [min_freq max_freq];
f.timerange = handles.timescale;
if isfield(handles,'in') && isfield(handles,'out')
    f.in = handles.in;
    f.out = handles.out;
else
    f.in = handles.timescale(1);
    f.out = handles.timescale(2);
end
f.angles = str2array(get(handles.edit_angles, 'String'))*pi;
f.fs = handles.fs;
f.degree = handles.degree;
if isfield(handles,'f')
    f.notmatfile = handles.f;
end
% THIS IS SPECIFICALLY FOR KYLE'S EMG DATA
recfname = strcat(handles.filename(1:end-4),'rec');
if exist(recfname,'file')
    f.recfile = recfname;
    f.rec = readrecf(recfname);
    THRESHOLD = 1000; % Stim trial if peak emg is above this value
    in = round((handles.in-handles.buffer)*handles.fs/1000);
    out = round((handles.out+handles.buffer)*handles.fs/1000);
    outsound = round(handles.out*handles.fs/1000);
    [emg emgds] = evsoundin('',handles.filename,'obs1'); % EMG DATA ON CH2
    f.emg = emg(in:out);
    if max(emg(in:outsound)) > THRESHOLD
        f.type = 'stim';
    else
        f.type = 'catch';
    end
end


% --- Executes just before analyzecontours is made visible.
function analyzecontours_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyzecontours (see VARARGIN)

% Prompt user for input file

handles.all_files = [];
handles.polynomialmode = false;

while length(handles.all_files) < 1
    pathname = uigetdir;
    if pathname == 0
        error('Did not load a directory');
    end
    handles.all_files = [dir(fullfile(pathname,'*.cbin')) dir(fullfile(pathname,'*.compiled.mat'))];
    if length(handles.all_files) < 1
        disp('No valid files found');
    end
end
handles.current_file = 1;
handles.pathname = pathname;

if ~isfield(handles,'degree')
    handles.degree = 3; % Default value for polynomial degree
end

handles = loadFiles(handles,hObject);

set(handles.status_bar, 'String', '');

handles.current_syllable = 0;

generate_spectrogram(handles, hObject);

% UIWAIT makes analyzecontours wait for user response (see UIRESUME)
% uiwait(handles.maincanvas);


% --- Outputs from this function are returned to the command line.
function varargout = analyzecontours_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in plot_popup.
function plot_popup_Callback(hObject, eventdata, handles)
% hObject    handle to plot_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_popup

val = get(hObject, 'Value');
str = get(hObject, 'String');
switch str{val}
    case 'peaks'
        handles.current_data = handles.peaks;
    case 'membrane'
        handles.current_data = handles.membrane;
    case 'sinc'
        handles.current_data = handles.sinc;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in graph_pushbutton.
function graph_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to graph_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

generate_spectrogram(handles, hObject);


% --------------------------------------------------------------------
function open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt user for input file

handles.all_files = [];

pathname = '/';

while pathname ~= 0 % Ends if you click "cancel"
    pathname = uigetdir;
    handles.all_files = [dir(fullfile(pathname,'*.cbin')) dir(fullfile(pathname,'*.compiled.mat'))];
    if length(handles.all_files) > 0
        break;
    else
        disp('No valid files found');
    end
end
handles.current_file = 1;
handles.pathname = pathname;

handles = loadFiles(handles,hObject);

generate_spectrogram(handles, hObject);


function edit_nfft_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nfft as text
%        str2double(get(hObject,'String')) returns contents of edit_nfft as a double


% --- Executes during object creation, after setting all properties.
function edit_nfft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_overlap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_overlap as text
%        str2double(get(hObject,'String')) returns contents of edit_overlap as a double


% --- Executes during object creation, after setting all properties.
function edit_overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_arthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_arthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_arthreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_arthreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_arthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_timescales_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timescales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timescales as text
%        str2double(get(hObject,'String')) returns contents of edit_timescales as a double


% --- Executes during object creation, after setting all properties.
function edit_timescales_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timescales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_minfreq as a double


% --- Executes during object creation, after setting all properties.
function edit_minfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_maxfreq as a double


% --- Executes during object creation, after setting all properties.
function edit_maxfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_angles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_angles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_angles as text
%        str2double(get(hObject,'String')) returns contents of edit_angles as a double


% --- Executes during object creation, after setting all properties.
function edit_angles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_angles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggle_contour.
function toggle_contour_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_contour

generate_spectrogram(handles, hObject);


function edit_vibrancy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vibrancy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vibrancy as text
%        str2double(get(hObject,'String')) returns contents of edit_vibrancy as a double


% --- Executes during object creation, after setting all properties.
function edit_vibrancy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vibrancy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_syllable_Callback(hObject, eventdata, handles)
% hObject    handle to edit_syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_syllable as text
%        str2double(get(hObject,'String')) returns contents of edit_syllable as a double


% --- Executes during object creation, after setting all properties.
function edit_syllable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prev_syllable.
function prev_syllable_Callback(hObject, eventdata, handles)
% hObject    handle to prev_syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.polynomialmode % Need to save this variable locally because it is set to false automatically later on
    drawpolynomial = true;
else
    drawpolynomial = false;
end
if handles.current_syllable > 1
    handles.current_syllable = handles.current_syllable - 1;
else
    if handles.current_file > 1
        handles.current_file = handles.current_file - 1;
        handles = loadFiles(handles,hObject);
        if isfield(handles,'f')
            if length(handles.current_syllable) > 1
                handles.current_syllable = 1;
            end
            syllable = get(handles.edit_syllable, 'String');
            counter = 0;
            for i=(1:numel(handles.f.onsets))
                if sum(strfind(syllable,handles.f.labels(i)))
                    counter = counter + 1;
                end
            end
            handles.current_syllable = counter;
        end
    end
end
generate_spectrogram(handles, hObject);
if drawpolynomial
    fit_poly_ClickedCallback(hObject, eventdata, handles, handles.degree);
end

% --- Executes on button press in next_syllable.
function next_syllable_Callback(hObject, eventdata, handles)
% hObject    handle to next_syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.polynomialmode % Need to save this variable locally because it is set to false automatically later on
    drawpolynomial = true;
else
    drawpolynomial = false;
end
if ~isfield(handles,'f')
    nextfile_Callback(hObject, eventdata, handles);
elseif isfield(handles,'max_counter') && handles.current_syllable < handles.max_counter
    handles.current_syllable = handles.current_syllable + 1;
    generate_spectrogram(handles, hObject);
else
    nextfile_Callback(hObject, eventdata, handles);
end
if drawpolynomial
    fit_poly_ClickedCallback(hObject, eventdata, handles, handles.degree);
end

% --------------------------------------------------------------------
function point_selector_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to point_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function mainaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function mainaxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainaxes


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function maincanvas_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to maincanvas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function vibrancy_slider_Callback(hObject, eventdata, handles)
% hObject    handle to vibrancy_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
generate_spectrogram(handles, hObject);


% --- Executes during object creation, after setting all properties.
function vibrancy_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vibrancy_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --------------------------------------------------------------------
function savedata_ClickedCallback(hObject, eventdata, handles)
% SAVE DATA METHOD
% hObject    handle to savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
def_name = handles.filename;
if strcmp(def_name(end-3:end),'cbin')
    def_name = def_name(1:end-5);
end
def_name = strcat(def_name,'.mat');
[filename, pathname] = uiputfile('*.mat','Save data', def_name);
handles.save_path = pathname;
if ~isequal(filename,0) && ~isequal(pathname,0)
    f = getfields(handles);
    [p S] = img_to_cart(f, handles.degree);
    f.polynomial = p;
    f.polynomialinfo = S;
    save(fullfile(pathname,filename),'-struct','f');
    set(handles.status_bar, 'String', strcat({'Saved data to '}, filename));
else
    set(handles.status_bar, 'String', 'Did not save data');
end


% --------------------------------------------------------------------
function select_points_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to select_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x, y] = ginput;

points = [y,x];

generate_spectrogram(handles, hObject, points);


% --------------------------------------------------------------------
function reload_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

generate_spectrogram(handles, hObject, 0);


% --- Executes on button press in prevfile.
function prevfile_Callback(hObject, eventdata, handles)
% hObject    handle to prevfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_file > 1
    handles.current_file = handles.current_file - 1;
    handles = loadFiles(handles,hObject);
    if length(handles.current_syllable) > 0
        handles.current_syllable = 1;
    end
    generate_spectrogram(handles, hObject);
end


% --- Executes on button press in nextfile.
function nextfile_Callback(hObject, eventdata, handles)
% hObject    handle to nextfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_file < length(handles.all_files)
    handles.current_file = handles.current_file + 1;
    handles = loadFiles(handles,hObject);
    if length(handles.current_syllable) > 0
        handles.current_syllable = 1;
    end
    generate_spectrogram(handles, hObject);
end


% --- Executes on selection change in files_menu.
function files_menu_Callback(hObject, eventdata, handles)
% hObject    handle to files_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files_menu
val = get(hObject,'Value');
str = get(hObject,'String');
handles.current_file = val;
handles = loadFiles(handles,hObject);
if length(handles.current_syllable) > 0
    handles.current_syllable = 1;
end
generate_spectrogram(handles,hObject);


% --- Executes during object creation, after setting all properties.
function files_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over files_menu.
function files_menu_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to files_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in syllable_bars.
function syllable_bars_Callback(hObject, eventdata, handles)
% hObject    handle to syllable_bars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of syllable_bars
generate_spectrogram(handles, hObject, 0);


% --------------------------------------------------------------------
function fit_poly_ClickedCallback(hObject, eventdata, handles, degree)
% hObject    handle to fit_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update = false; % Default don't update handles object
if ~exist('degree','var') || isempty(degree) % Allows degree to be passed as arg
    degree = inputdlg('Polynomial Degree: ', 'Enter Degree', 1, {num2str(handles.degree)});
    if length(degree) > 0
        degree = degree{1}; % inputdlg returns cell type
    end
end
if length(degree) > 0
    if isnumeric(degree) || str2num(degree)
        handles.polynomialmode = true;
        guidata(hObject, handles);
        FONTSIZE = 12;
        if ~isnumeric(degree)
            degree = round(str2num(degree));
        else
            degree = round(degree);
        end
        if handles.degree ~= degree
            handles.degree = degree;
            guidata(hObject, handles); % Update handles
        end
        f = getfields(handles);
        [p S] = img_to_cart(f, handles.degree); % Spline implementation: img_to_spline(f, handles.degree);
        set(handles.status_bar, 'String', strcat({'R-Squared: '},num2str(S.rsquared)));
        polydata = S.polydata;
        polydata = (polydata - f.freqrange(1))*S.ylen/(f.freqrange(2)-f.freqrange(1));
        polydata = round(polydata);
        S.x = S.x * S.xlen / (f.out - f.in) + f.mapbounds(1);
        S.x = round(S.x);
        mydata = f.data;
        level = -1 * max(max(mydata)); % Value of points along polynomial
        for i=1:numel(polydata)
            if polydata(i) > 0 && polydata(i) < length(mydata(:,i))
                mydata(polydata(i),S.x(i)) = level;
            end
        end
        imagesc(f.timerange,f.freqrange,mydata);
        set(gca,'YDir','normal');
        xlabel('time [ms]', 'FontSize', FONTSIZE);
        ylabel('frequency [kHz]', 'FontSize', FONTSIZE);
    else
        set(handles.status_bar, 'String', 'Invalid degree for polynomial.');
    end
else
    set(handles.status_bar, 'String', 'Did not make polynomial.');
end

% --- Executes on button press in save_and_next.
function save_and_next_Callback(hObject, eventdata, handles)
% hObject    handle to save_and_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save data
if ~isfield(handles,'save_path')
    handles.save_path = uigetdir;
end
if handles.polynomialmode % Need to save this variable locally because it is set to false automatically later on
    drawpolynomial = true;
else
    drawpolynomial = false;
end
filename = handles.filename;
if strcmp(filename(end-3:end),'cbin')
    filename = filename(1:end-5);
end
f = getfields(handles);
[p S] = img_to_cart(f, handles.degree);
f.polynomial = p;
f.polynomialinfo = S;
filename = strcat(filename,'.',num2str(int32(f.in)),'.DATA','.mat');
save(fullfile(handles.save_path,filename),'-struct','f');
[path name ext] = fileparts(handles.save_path);
set(handles.status_bar, 'String', strcat({'Saved data to '}, name));
% Move to next image
next_syllable_Callback(hObject, eventdata, handles);
if drawpolynomial
    fit_poly_ClickedCallback(hObject, eventdata, handles, handles.degree);
end

function channel_val_Callback(hObject, eventdata, handles)
% hObject    handle to channel_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_val as text
%        str2double(get(hObject,'String')) returns contents of channel_val as a double


% --- Executes during object creation, after setting all properties.
function channel_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
