function [data, fs] = load_sound(directory_path, file_name, channel)
%LOAD_SOUND General handler for loading sounds
% Currently supports:
%   -cbin
%   -mat with structs:
%     -> data: song data
%     -> fs: 

if strcmp(file_name(end-3:end),'cbin')
    % Load the data from the cbin file
    [data fs] = evsoundin(directory_path, file_name, channel);
elseif strcmp(handles.filename(end-2:end),'mat')
    file = load(fullfile(directory_path, file_name));
    if ~isfield(file,'data') || ~isfield(file,'fs')
        error('Error loading sound: .mat file must have the attributes "data" (song data) and "fs" (sampling frequency)');
    end
    data = file.data;
    fs = file.fs;
end

% Check to make sure the file loaded correctly
if ~exist('data','var') || ~exist('fs','var')
    error('load_sound: Error loading data');
elseif isempty(data) || numel(data) == 0
    error('load_sound: Data file empty');
end

end

