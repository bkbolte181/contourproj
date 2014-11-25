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

function handles = load_files(handles, hObject)
% Load a particular file, updating the handles structure

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

function handles = generate_spectrogram(handles, hObject, contourPoint)
% Function to actually generate the spectrogram

channel = get(handles.channel_val,'String');

% Load the sound using an external method
[dat, fs] = load_sound(handles.pathname, handles.filename,channel);
handles.dat = dat;

buffer = 30; % Buffer before and after syllable - change this if desired, must be integer
handles.buffer = buffer;
handles.fs = fs;
handles.polynomialmode = false;

% Get the current syllable
if length(get(handles.edit_syllable,'String')) > 0
    if handles.current_syllable == 0
        handles.current_syllable = 1;
    end
else
    handles.current_syllable = 0;
    handles.in = 0;
    handles.out = 0;
end

% Load data from the specified file
if ~isfield(handles,'f');
    % If the file field is empty, set the current syllable to zero
    handles.current_syllable = 0;
elseif handles.current_syllable > 0
    % This is the bit of code that handles the edit_syllable input
    syllable = get(handles.edit_syllable, 'String');
    
    % Count the number of syllables
    counter = 0;
    
    % Loop through each syllable in the .cbin file
    for i=(1:numel(handles.f.onsets))
        % If the syllable's label is one of the labels specified by the
        % user, increment counter
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

handles.freqrange = [min_freq, max_freq];

new = {N, OVERLAP, Threshold, sigmas, angles};

if exist('contourPoint','var') && ~any(any(contourPoint))
    reset = 'yes';
end

% Huge list of conditions. If one of them is satisfied
[handles.contourmap, handles.sonogram, f, t] = zf_contour_func(dat, fs, N, OVERLAP, Threshold, sigmas, [min_freq max_freq], angles);
[handles.N, handles.OVERLAP, handles.Threshold, handles.sigmas, handles.angles] = deal(new{:});
handles.times = t;
handles.freqs = f;

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
handles.timescale = [min_time, max_time];
vibrancy = get(handles.vibrancy_slider, 'Value');
handles.FONTSIZE = 12;

vibrancy = 20^(vibrancy-1);
handles.vibrancy = vibrancy;

% Choose default command line output for analyzecontours
handles.output = hObject;
handles = draw_spectrogram(handles);

% Update handles structure
guidata(hObject, handles);

function handles = draw_spectrogram(handles)
%DRAW_SPECTROGRAM Method for actually drawing the spectrogram

% Draws the contour image
if get(handles.toggle_contour, 'Value'); % Generate contours
    imagesc(handles.timescale,handles.freqrange,log(handles.contourmap+handles.vibrancy));
else % Generate regular spectrogram
    imagesc(handles.timescale,handles.freqrange, log(abs(handles.sonogram)+handles.vibrancy));
end
set(gca,'YDir','normal');
xlabel('time [ms]', 'FontSize', handles.FONTSIZE);
ylabel('frequency [kHz]', 'FontSize', handles.FONTSIZE);

colormap('hot');

function f = getfields(handles)
% Get all available fields and store in a structured array 'f'
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

handles = load_files(handles,hObject);

set(handles.status_bar, 'String', '');

handles.current_syllable = 0;

handles = generate_spectrogram(handles, hObject);
handles = draw_spectrogram(handles);

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

handles = generate_spectrogram(handles, hObject);
handles = draw_spectrogram(handles);

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

handles = load_files(handles,hObject);

handles = generate_spectrogram(handles, hObject);
handles = draw_spectrogram(handles);

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

handles = draw_spectrogram(handles);

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
        handles = load_files(handles,hObject);
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
handles = generate_spectrogram(handles, hObject);
handles = draw_spectrogram(handles);
if drawpolynomial
    fit_poly_ClickedCallback(hObject, eventdata, handles, handles.degree);
end

% --- Executes on button press in next_syllable.
function next_syllable_Callback(hObject, eventdata, handles)
% hObject    handle to next_syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'f')
    nextfile_Callback(hObject, eventdata, handles);
elseif isfield(handles,'max_counter') && handles.current_syllable < handles.max_counter
    handles.current_syllable = handles.current_syllable + 1;
    handles = generate_spectrogram(handles, hObject);
    handles = draw_spectrogram(handles);
else
    % ALTERATION
    % This is currently set to scroll backwards - change to
    % nextfile_Callback to scroll forwards
    prevfile_Callback(hObject, eventdata, handles);
end


% --- Executes on slider movement.
function vibrancy_slider_Callback(hObject, eventdata, handles)
% hObject    handle to vibrancy_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = generate_spectrogram(handles, hObject);
handles = draw_spectrogram(handles);


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
function reload_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = generate_spectrogram(handles, hObject, 0);

handles = draw_spectrogram(handles);


% --- Executes on button press in prevfile.
function prevfile_Callback(hObject, eventdata, handles)
% hObject    handle to prevfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_file > 1
    handles.current_file = handles.current_file - 1;
    handles = load_files(handles,hObject);
    if length(handles.current_syllable) > 0
        handles.current_syllable = 1;
    end
    handles = generate_spectrogram(handles, hObject);
    handles = draw_spectrogram(handles);
end


% --- Executes on button press in nextfile.
function nextfile_Callback(hObject, eventdata, handles)
% hObject    handle to nextfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_file < length(handles.all_files)
    handles.current_file = handles.current_file + 1;
    handles = load_files(handles,hObject);
    if length(handles.current_syllable) > 0
        handles.current_syllable = 1;
    end
    handles = generate_spectrogram(handles, hObject);
end
handles = draw_spectrogram(handles);

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
handles = load_files(handles,hObject);
if length(handles.current_syllable) > 0
    handles.current_syllable = 1;
end
handles = generate_spectrogram(handles,hObject);
handles = draw_spectrogram(handles);

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

% --- Executes on button press in syllable_bars.
function syllable_bars_Callback(hObject, eventdata, handles)
% hObject    handle to syllable_bars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of syllable_bars
handles = generate_spectrogram(handles, hObject, 0);
handles = draw_spectrogram(handles);

% --- Executes on button press in save_and_next.
function save_and_next_Callback(hObject, eventdata, handles)
% hObject    handle to save_and_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save and continue button

% Save data
if ~isfield(handles,'save_path')
    % If no current save destination, prompt user for new save destination
    handles.save_path = uigetdir;
end

% Get the file's name from the handles structure
filename = handles.filename;

% Remove the '.cbin' part from the end of the file's name
if strcmp(filename(end-3:end),'cbin')
    filename = filename(1:end-5);
end

% Basically get all needed fields from the handles structure
f = getfields(handles);

% Draw polynomial with the given file
[p S] = img_to_cart(f, handles.degree);
f.polynomial = p;
f.polynomialinfo = S;

% Make a new file name with the old one and save it
filename = strcat(filename,'.',num2str(int32(f.in)),'.DATA','.mat');
save(fullfile(handles.save_path,filename),'-struct','f');

% Update the status bar with the relevant information
[path name ext] = fileparts(handles.save_path);
set(handles.status_bar, 'String', strcat({'Saved data to '}, name));

% Move to next image
next_syllable_Callback(hObject, eventdata, handles);

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

% Some functions required so that the GUI doesn't throw errors
function mainaxes_CreateFcn(hObject, eventdata, handles)
% Nothing here

function maincanvas_WindowButtonDownFcn(hObject, eventdata, handles)
% Nothing here

function edit_syllable_Callback(hObject, eventdata, handles)
% Nothing here

function edit_overlap_Callback(hObject, eventdata, handles)
% Nothing here

function edit_nfft_Callback(hObject, eventdata, handles)
% Nothing here

function edit_arthreshold_Callback(hObject, eventdata, handles)
% Nothing here

function edit_minfreq_Callback(hObject, eventdata, handles)
% Nothing here
