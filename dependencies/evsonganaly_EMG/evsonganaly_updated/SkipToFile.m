function varargout = SkipToFile(varargin)
% SKIPTOFILE M-file for SkipToFile.fig
%      SKIPTOFILE, by itself, creates a new SKIPTOFILE or raises the existing
%      singleton*.
%
%      H = SKIPTOFILE returns the handle to a new SKIPTOFILE or the handle to
%      the existing singleton*.
%
%      SKIPTOFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SKIPTOFILE.M with the given input arguments.
%
%      SKIPTOFILE('Property','Value',...) creates a new SKIPTOFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SkipToFile_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SkipToFile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help SkipToFile

% Last Modified by GUIDE v2.5 09-Nov-2004 14:37:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SkipToFile_OpeningFcn, ...
                   'gui_OutputFcn',  @SkipToFile_OutputFcn, ...
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

% --- Executes just before SkipToFile is made visible.
function SkipToFile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SkipToFile (see VARARGIN)

% Choose default command line output for SkipToFile
handles.output = hObject;

%load up the filenames into the window
%need the handles.FILENAMES structure array from evsonganaly
tempvar = varargin{1};

if isstruct(tempvar) & isfield(tempvar,'fname')
    % i.e. if SkipToFile was called by evsonganaly, which sends a struct
    handles.FILENAMES={tempvar(:).fname}';
elseif iscell(tempvar)
    % i.e. if SkipToFile was called by uievtafsim, which sends a cell array
    handles.FILENAMES=tempvar;
end

handles.NFILE = varargin{2};
handles.ORIG_NFILE = handles.NFILE;
handles.IsCancel = 0;
guidata(hObject,handles);

set(handles.FilesListBox,'Value',1);
set(handles.FilesListBox,'String',handles.FILENAMES);
set(handles.FilesListBox,'Value',handles.NFILE);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SkipToFile wait for user response (see UIRESUME)
uiwait(handles.SkipToFile);
return;

% --- Outputs from this function are returned to the command line.
function varargout = SkipToFile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.NFILE;
varargout{2} = handles.IsCancel;
delete(handles.SkipToFile);
return;

% --- Executes on selection change in FilesListBox.
function FilesListBox_Callback(hObject, eventdata, handles)
% hObject    handle to FilesListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FilesListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FilesListBox

stype=get(handles.SkipToFile,'SelectionType');
if (strcmp(stype,'open'))
	handles.IsCancel=0;
	handles.NFILE = get(handles.FilesListBox,'Value');
	guidata(hObject,handles);
	uiresume(handles.SkipToFile);
elseif (strcmp(stype,'normal'))
	handles.NFILE = get(handles.FilesListBox,'Value');
	guidata(hObject,handles);
end
return;

% --- Executes during object creation, after setting all properties.
function FilesListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilesListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
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

handles.IsCancel=0;
handles.NFILE=get(handles.FilesListBox,'Value');
guidata(hObject,handles);
uiresume(handles.SkipToFile);
return;

% --- Executes on button press in CancelBtn.
function CancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.IsCancel=1;
handles.NFILE=handles.ORIG_NFILE;
guidata(hObject,handles);
uiresume(handles.SkipToFile);
return;

% --- Executes on button press in FindLastLabeledFile.
function FindLastLabeledFile_Callback(hObject, eventdata, handles)
% hObject    handle to FindLastLabeledFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%go through the files and try to find the last file with a .not.mat
%file associated with it
% makes it easier to continue your labeling session

newval = 1;
startval = get(handles.FilesListBox,'Value');
if (startval<1)
    startval=1;
end
fnames = handles.FILENAMES;
for ii=startval:length(fnames)
    fn=fnames{ii};
    if (~exist([fn,'.not.mat'],'file'))
        newval = ii;
        break;
    end
end
set(handles.FilesListBox,'Value',newval);
guidata(hObject,handles);
return;
