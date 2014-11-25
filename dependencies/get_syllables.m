function [songs] = get_syllables(varargin)
%GET_SYLLABLES Get all syllables associated with a file or set of files
% Declare:
%   -> Structure or cell containing files
%
%
%
%
%

files = [];

% Parse input arguments
[files, buffer, directory, channel] = parse_argin(varargin);

% Get the data by using the parse_file method for each file
songs = cell(1, numel(files));
for i=1:length(files)
    file = files{i};
    songs{i} = parse_file(file, buffer, directory, channel);
    disp(sprintf('Processing %s...', file));
end

% Prompt user to save the structured data
[filename, pathname] = uiputfile({'*.mat','MATLAB File'},'Save Syllable Data',fullfile(pwd,'songs.mat'));

if filename
    save(fullfile(pathname, filename), 'songs');
else
    disp('Did not save file.');
end

end

function [song] = parse_file(file, buffer, directory, channel);
% Function for parsing each file
%
%
%
%
%

[data, fs] = load_sound(directory, file, channel);

% Store some attributes in the struct
song.fs = fs;
song.file = file;
song.directory = directory;
song.channel = channel;
song.buffer = buffer;

% Load the associated matfile
matfilename = [file, '.not.mat'];
if ~exist(matfilename, 'file')
    disp(sprintf('load_sound:File %s does not exist. Skipping...', matfilename));
    return
end
matfile = load(fullfile(directory, matfilename));

song.labels = matfile.labels;

syllables = cell(1,length(matfile.labels));
for i=1:length(matfile.labels)
    % Get onset and offset in terms of samples, given milliseconds
    on = round(matfile.onsets(i) * fs / 1000);
    off = round(matfile.offsets(i) * fs / 1000);
    
    % Store timing information as fields in the syllable structure
    syllables{i}.on = on;
    syllables{i}.on_ms = matfile.onsets(i);
    syllables{i}.on_buffer = on - buffer;
    syllables{i}.on_buffer_ms = matfile.onsets(i) - buffer*1000/fs;
    syllables{i}.off = off;
    syllables{i}.off_ms = matfile.offsets(i);
    syllables{i}.off_buffer = off + buffer;
    syllables{i}.off_buffer_ms = matfile.offsets(i) + buffer*1000/fs;
    syllables{i}.label = matfile.labels(i);
    
    % Get the song data associated with a syllable
    syllables{i}.data = data(max(on-buffer,1):min(off+buffer,length(data)));
end

song.syllables = syllables;

end

function [files, buffer, directory, channel] = parse_argin(varargin)
% Helper function for parsing the input arguments
%
%
%
%
%
%

v = deal(varargin{:});

% Set default values
buffer = 0;
directory = pwd;
channel = 'obs0';

% The first variable can be the filename
if numel(v) > 0
    f = v{1};
    if isstruct(f)
        files = load_struct(f);
    elseif iscell(f)
        files = f;
    else
        files = load_struct(dir('*.cbin'));
    end
else
    files = load_struct(dir('*.cbin'));
end

% This is the main method for actually processing the arguments
for i=1:length(v)
    if ischar(v{i})
        % Parse "buffer" argument
        if strcmp(v{i}, 'buffer')
            if i <= numel(v) + 1
                if ischar(v{i+1})
                    buffer = round(str2num(v{i+1}));
                elseif isnumeric(v{i+1})
                    buffer = round(v{i+1});
                else
                    error('parse_argin:"buffer" must be a numeric or string value');
                end
            else
                error('parse_argin:Did not specify "buffer"');
            end
        end
        
        % Parse "directory" argument
        if strcmp(v{i}, 'directory')
            if i <= numel(v) + 1
                if ischar(v{i+1}) && isdir(v{i+1})
                    directory = v{i+1};
                else
                    error('parse_argin:"directory" must be a string');
                end
            else
                error('parse_argin:Did not specify "directory"');
            end
        end
        
        % Parse "channel" argument
        if strcmp(v{i}, 'channel')
            if i <= numel(v) + 1
                if ischar(v{i+1})
                    channel = v{i+1};
                else
                    error('parse_argin:"channel" must be a string');
                end
            else
                error('parse_argin:Did not specify "channel"');
            end
        end
    end
end

% Post-process checking
if buffer < 0
    error('parse_argin:You cannot specify a negative value for "buffer"');
end

end

function [files] = load_struct(d)
% Load data from struct, for example the result of the "dir" function
% It is assumed that the name of each file is the "name" attribute of the
% struct
%
%
%
%
%

files = cell(1,numel(d));
for i=1:numel(d)
    files{i} = d(i).name;
end

end
