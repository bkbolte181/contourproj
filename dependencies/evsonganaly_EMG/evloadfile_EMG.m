function varargout = evloadfile_EMG(varargin)
% EVLOADFILE M-file for evloadfile.fig
%      EVLOADFILE, by itself, creates a new EVLOADFILE or raises the existing
%      singleton*.
%
%      H = EVLOADFILE returns the handle to a new EVLOADFILE or the handle to
%      the existing singleton*.
%
%      EVLOADFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVLOADFILE.M with the given input arguments.
%
%      EVLOADFILE('Property','Value',...) creates a new EVLOADFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evloadfile_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evloadfile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help evloadfile

% Last Modified by GUIDE v2.5 24-May-2012 11:09:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evloadfile_OpeningFcn, ...
                   'gui_OutputFcn',  @evloadfile_OutputFcn, ...
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


% --- Executes just before evloadfile is made visible.
function evloadfile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evloadfile (see VARARGIN)

% Choose default command line output for evloadfile
handles.output = hObject;

%set some defaults here
CURDIR=pwd;
%set(handles.DirTxtBox,'String',CURDIR);
%set(handles.DirTxtBox,'String','/Users/evren/tempdata');
dirtmp = pwd;
set(handles.DirTxtBox,'String',dirtmp);
set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Max'));
set(handles.ExtTxtBox,'String','.cbin');
set(handles.OKBtn,'Value',get(handles.OKBtn,'Min'));
set(handles.CancelBtn,'Value',get(handles.CancelBtn,'Min'));
set(handles.FindDirBtn,'Value',get(handles.FindDirBtn,'Min'));
set(handles.FindFileBtn,'Value',get(handles.FindFileBtn,'Min'));
set(handles.FileTypePulDown,'Value',2);
set(handles.evloadfile,'Selected','on');
set(handles.ChanSpec,'String','0');
set(handles.EMGChanSpec,'String','1');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes evloadfile wait for user response (see UIRESUME)
uiwait(handles.evloadfile);
return;

% --- Outputs from this function are returned to the command line.
function varargout = evloadfile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output=handles.FILENAMELIST;
varargout{1} = handles.output;
varargout{2} = get(handles.ChanSpec,'String');
varargout{3} = get(handles.EMGChanSpec,'String');
delete(handles.evloadfile);
return;

% --- Executes on button press in FindFileBtn.
function FindFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to FindFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% looks for a batch file or a single data file
[fname,pathname]=uigetfile('*','Pick a FILE :');
FILENAME=fullfile(pathname,fname);
set(handles.FileTxtBox,'String',FILENAME);
set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Min'));
[atmp,btmp,ext]=fileparts(fname);
if (strcmp(ext,'.wav')|strcmp(ext,'.ebin')|strcmp(ext,'.cbin')|strcmp(ext,'.bbin'))
    set(handles.FileTypePulDown,'Value',3);
else
    set(handles.FileTypePulDown,'Value',1);
end
set(handles.FindDirBtn,'Value',get(handles.FindDirBtn,'Min'));
guidata(hObject, handles);
return;

function FileTxtBox_Callback(hObject, eventdata, handles)
% hObject    handle to FileTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileTxtBox as text
%        str2double(get(hObject,'String')) returns contents of FileTxtBox as a double
return;

% --- Executes during object creation, after setting all properties.
function FileTxtBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes on selection change in FileTypePulDown.
function FileTypePulDown_Callback(hObject, eventdata, handles)
% hObject    handle to FileTypePulDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FileTypePulDown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileTypePulDown

if (get(hObject,'Value')==2)
    set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Max'));   
else
    set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Min'));
end
return;

% --- Executes during object creation, after setting all properties.
function FileTypePulDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileTypePulDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function DirTxtBox_Callback(hObject, eventdata, handles)
% hObject    handle to DirTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirTxtBox as text
%        str2double(get(hObject,'String')) returns contents of DirTxtBox as
%        a double

set(handles.FileTypePulDown,'Value',2);
set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Max'));
return;

% --- Executes during object creation, after setting all properties.
function DirTxtBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in FindDirBtn.
function FindDirBtn_Callback(hObject, eventdata, handles)
% hObject    handle to FindDirBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FindDirBtn
curdir=pwd;
pathname=uigetdir('','Pick a Directory :');
set(handles.DirTxtBox,'String',pathname);
set(handles.UseDirChkBox,'Value',get(handles.UseDirChkBox,'Max'));
set(handles.FileTypePulDown,'Value',2);
set(handles.FindDirBtn,'Value',get(handles.FindDirBtn,'Min'));
guidata(hObject, handles);
return;

% --- Executes on button press in UseDirChkBox.
function UseDirChkBox_Callback(hObject, eventdata, handles)
% hObject    handle to UseDirChkBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseDirChkBox
if (get(hObject,'Value')==get(hObject,'Max'))
    set(handles.FileTypePulDown,'Value',2);
else
    set(handles.FileTypePulDown,'Value',1);
end
return;

% --- Executes on button press in OKBtn.
function OKBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OKBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of OKBtn

psep = filesep;
if (strcmp(computer,'MAC')&(isunix))
    psep = '/';
end
if (get(handles.UseDirChkBox,'Value')==get(handles.UseDirChkBox,'Max'))
    % Get all the files with the given extension in the given directory
    dirname=get(handles.DirTxtBox,'String');        
    if (~strcmp(dirname(end),psep))
        dirname=[dirname,psep];
    end
    
    extname=get(handles.ExtTxtBox,'String');
    p = findstr(extname,'.');
    if (length(p)==0)
        if (strcmp(extname(1),'*'))
            extname = ['*.',extname(2:end)];
        else
            extname = ['*.',extname];
        end
    else
        extname = ['*',extname(p(end):end)];
    end
    dirfiles = dir([dirname,extname]);
    FILENAMELIST=[];
    for ind=1:length(dirfiles)
            FILENAMELIST(ind).fname = [dirname,dirfiles(ind).name];
    end
else
    % Use the batch file or the single file given in the dir
    filename=get(handles.FileTxtBox,'String');
    if (get(handles.FileTypePulDown,'Value')==1)
        %batch file
        batchfilename=filename;
        [fid,msg] = fopen(batchfilename,'r');
        if (fid==-1)
            errordlg(msg);
            return;
        end
        cnt=0;
        while (1)
            fn=fgetl(fid);
            if (~ischar(fn))
                break;
            end
            [dirn,fntmp,ext]=fileparts(fn);

            if (length(dirn)==0)
                dirn = fileparts(batchfilename);
                dirn=[dirn,psep];
            end
            
            if (~strcmp(dirn(end),psep))
                dirn=[dirn,psep];
            end
            cnt = cnt + 1;
            FILENAMELIST(cnt).fname = [dirn,fntmp,ext];
        end
        fclose(fid);
    else
        %single file
        FILENAMELIST(1).fname=filename;
    end
end

handles.FILENAMELIST=FILENAMELIST;
guidata(hObject, handles);
uiresume(handles.evloadfile);
return;

% --- Executes on button press in CancelBtn.
function CancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CancelBtn
handles.FILENAMELIST=[];
guidata(hObject, handles);
uiresume(handles.evloadfile);
return;

function ExtTxtBox_Callback(hObject, eventdata, handles)
% hObject    handle to ExtTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExtTxtBox as text
%        str2double(get(hObject,'String')) returns contents of ExtTxtBox as a double
return;

% --- Executes during object creation, after setting all properties.
function ExtTxtBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExtTxtBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

function ChanSpec_Callback(hObject, eventdata, handles)
% hObject    handle to ChanSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChanSpec as text
%        str2double(get(hObject,'String')) returns contents of ChanSpec as a double
return;

% --- Executes during object creation, after setting all properties.
function ChanSpec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChanSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;



function EMGChanSpec_Callback(hObject, eventdata, handles)
% hObject    handle to EMGChanSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EMGChanSpec as text
%        str2double(get(hObject,'String')) returns contents of EMGChanSpec as a double


% --- Executes during object creation, after setting all properties.
function EMGChanSpec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMGChanSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
