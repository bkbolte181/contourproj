function varargout = ChangeSettings(varargin)
% CHANGESETTINGS M-file for ChangeSettings.fig
%      CHANGESETTINGS, by itself, creates a new CHANGESETTINGS or raises the existing
%      singleton*.
%
%      H = CHANGESETTINGS returns the handle to a new CHANGESETTINGS or the handle to
%      the existing singleton*.
%
%      CHANGESETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHANGESETTINGS.M with the given input arguments.
%
%      CHANGESETTINGS('Property','Value',...) creates a new CHANGESETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChangeSettings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChangeSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ChangeSettings

% Last Modified by GUIDE v2.5 12-Nov-2004 14:03:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChangeSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @ChangeSettings_OutputFcn, ...
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


% --- Executes just before ChangeSettings is made visible.
function ChangeSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChangeSettings (see VARARGIN)

% Choose default command line output for ChangeSettings
handles.output = hObject;

tmpstruct = varargin{1};
handles.MINDUR = tmpstruct.mindur;
handles.MININT = tmpstruct.minint;
handles.SEGTH  = tmpstruct.segth;
handles.SM_WIN = tmpstruct.sm_win;
guidata(hObject,handles);
handles=guidata(hObject);

set(handles.MinDurBox,'String',num2str(handles.MINDUR));
set(handles.MinIntBox,'String',num2str(handles.MININT));
set(handles.SegThBox,'String',num2str(handles.SEGTH));
set(handles.SmWinBox,'String',num2str(handles.SM_WIN));

tmpstruct.mindur = handles.MINDUR;
tmpstruct.minint = handles.MININT;
tmpstruct.segth  = handles.SEGTH;
tmpstruct.sm_win = handles.SM_WIN;
handles.SAVESETTINGS = tmpstruct;
guidata(hObject, handles);

% UIWAIT makes ChangeSettings wait for user response (see UIRESUME)
uiwait(handles.ChangeSettings);

% --- Outputs from this function are returned to the command line.
function varargout = ChangeSettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if (isfield(handles,'SAVESETTINGS'))
    tmpstruct = handles.SAVESETTINGS;
    tmpstruct.DOIT = 1;
    varargout{1} = tmpstruct;
else
    tmpstruct = [];
    tmpstruct.DOIT = 0;
    varargout{1} = tmpstruct;
end
delete(handles.ChangeSettings);
return;

function MinDurBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinDurBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinDurBox as text
%        str2double(get(hObject,'String')) returns contents of MinDurBox as a double

tmpstruct = handles.SAVESETTINGS;
tmpstruct.mindur = str2num(get(handles.MinDurBox,'String'));
if (isnan(tmpstruct.mindur))
    set(handles.MinDurBox,'String',num2str(handles.MINDUR));
    tmpstruct.mindur = handles.MINDUR;
end
handles.SAVESETTINGS = tmpstruct;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function MinDurBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinDurBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

function SmWinBox_Callback(hObject, eventdata, handles)
% hObject    handle to SmWinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmWinBox as text
%        str2double(get(hObject,'String')) returns contents of SmWinBox as a double

tmpstruct = handles.SAVESETTINGS;
tmpstruct.sm_win = str2num(get(handles.SmWinBox,'String'));
if (isnan(tmpstruct.sm_win))
    set(handles.SmWinBox,'String',num2str(handles.SM_WIN));
    tmpstruct.sm_win = handles.SM_WIN;
end
handles.SAVESETTINGS = tmpstruct;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function SmWinBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmWinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

function SegThBox_Callback(hObject, eventdata, handles)
% hObject    handle to SegThBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SegThBox as text
%        str2double(get(hObject,'String')) returns contents of SegThBox as a double

tmpstruct = handles.SAVESETTINGS;
tmpstruct.segth = str2num(get(handles.SegThBox,'String'));
if (isnan(tmpstruct.segth))
    set(handles.SegThBox,'String',num2str(handles.SEGTH));
    tmpstruct.segth = handles.SEGTH;
end
handles.SAVESETTINGS = tmpstruct;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function SegThBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SegThBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

function MinIntBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinIntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinIntBox as text
%        str2double(get(hObject,'String')) returns contents of MinIntBox as a double

tmpstruct = handles.SAVESETTINGS;
tmpstruct.minint = str2num(get(handles.MinIntBox,'String'));
if (isnan(tmpstruct.minint))
    set(handles.MinIntBox,'String',num2str(handles.MININT));
    tmpstruct.minint = handles.MININT;
end
handles.SAVESETTINGS = tmpstruct;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function MinIntBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinIntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes on button press in OKBtn.
function OKBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OKBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;
return;

% --- Executes on button press in CancelBtn.
function CancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=rmfield(handles,'SAVESETTINGS');
guidata(hObject,handles);
uiresume;
return;
