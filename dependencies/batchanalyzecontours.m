% -----
% About
% -----

% Run analyzecontours in batch mode

clear;

% ----------
% Parameters
% ----------

% STIM PARAMS
Stim = struct(...
    'threshold',    [100 0 0], ... % Thresholds
    'angles',       [0, 0.1*pi, 0.2*pi, 0.9*pi], ... % Angles for contour
    'type',         'stim', ... % Type
    'end', -1);

% CATCH PARAMS
Catch = struct(...
    'threshold',    [100 0 0], ... % Thresholds
    'angles',       [0, 0.1*pi, 0.2*pi, 0.9*pi], ... % Angles for contour
    'type',         'catch', ... % Type
    'end', -1);

% Other
window_size = 2048;             % Size of STFT window
overlap = 2040;                 % Overlap in STFT
timescales = [1 1];
freq_range = [1.5 4];
syllable = 'bc';                % Target syllable
channel_with_song = 'obs0';     % Channel of cbin with song data
emg_channel = 'obs2';           % Channel with EMG data
buffer = 40;                    % Buffer before and after syllable, in ms
stim_threshold = 1000;          % Threshold to consider a trial a 'stim'
simple = 1;                     % Instead of finding contours, sum the power spectrum

cbin_dir = uigetdir;
save_dir = uigetdir;

% ------------
% Main routine
% ------------

cbin_files = dir(fullfile(cbin_dir,'*.cbin'));     % Load all the cbin files
for i=1:numel(cbin_files)
    % Load in the .cbin song data
    cbin_name = cbin_files(i).name;
    
    [data, fs] = evsoundin(cbin_dir,cbin_name,channel_with_song);
    
    % Load in the .cbin emg data
    [emg_data, emg_fs] = evsoundin(cbin_dir,cbin_name,emg_channel);
    
    % Load in the .not.mat file data
    notmat_name = strcat(cbin_name,'.not.mat');
    loc = fullfile(cbin_dir,notmat_name);
    
    if ~exist(loc,'file')
        disp(sprintf('%s does not exist', loc));
        continue;
    end
    
    notmat_file = load(fullfile(cbin_dir,notmat_name));
    
    % Find the desired syllable label
    syl_locations = regexp(notmat_file.labels,syllable)+length(syllable)-1;
    
    for j=1:numel(syl_locations)
        disp([i j]);
        % Find the onset and offset of the syllable, minus and plus buffer 
        onset = round((notmat_file.onsets(syl_locations(j)) - buffer) * fs / 1000);
        offset = round((notmat_file.offsets(syl_locations(j)) + buffer) * fs / 1000);
        
        % Get the song and emg data of this point
        syl_data = data(onset:offset);
        emg = emg_data(onset:offset);
        
        % What the file will be named
        save_file_as = strcat(cbin_name, '.onset-', num2str(notmat_file.onsets(syl_locations(j))), 'ms.mat');
        
        % Determine if it is a stim or catch trial
        if max(emg) > stim_threshold
            m = Stim;
        else
            m = Catch;
        end
        type = m.type;
        
        % Two different methods for finding the contour
        if simple
            [contourmap, sonogram, f, t] = simple_contour_from_syllable(syl_data, fs, window_size, overlap);
        else
            % Use the contour finding method
            [contourmap, sonogram, f, t] = zf_contour_func(syl_data, fs, window_size, overlap, m.threshold, timescales, freq_range, m.angles);
        end
        
        songdata = syl_data;
        
        timerange = [t(1), t(end)];
        freqrange = freq_range; %[min(f), max(f)];
        notmatfile = notmat_file;
        parameters = m;
        
        % Store timing information in a struct
        timing.onset = onset;
        timing.offset = offset;
        timing.onset_no_buffer = round((notmat_file.onsets(syl_locations(j))) * fs / 1000);
        timing.offset_no_buffer = round((notmat_file.offsets(syl_locations(j))) * fs / 1000);
        timing.ms_onset = notmat_file.onsets(syl_locations(j)) - buffer;
        timing.ms_offset = notmat_file.offsets(syl_locations(j)) + buffer;
        timing.ms_onset_no_buffer = notmat_file.onsets(syl_locations(j));
        timing.ms_offset_no_buffer = notmat_file.offsets(syl_locations(j));
        
%         % REMOVE THE BUFFER FROM THE BEGINNING OF EACH SYLLABLE
%         % Save the dimensions of the sonogram as variables
%         [xlen, ylen] = size(sonogram);
% 
%         % Find the first "voiced" point in the song
%         has_voice = 0;
%         k = 0;
%         while ~has_voice
%             k = k + 1;
%             if isvoiced(sonogram(:,k), contourmap(:,k))
%                 has_voice = 1;
%             end
%         end
% 
%         timerange(1) = timerange(1) + (timerange(2) - timerange(1)) * k / ylen;
% 
%         sonogram = sonogram(:,k:end);
%         contourmap = contourmap(:,k:end);
% 
%         l = length(songdata);
%         diff = timing.ms_offset - timing.ms_onset;
% 
%         songdata = songdata(round(l * k / ylen):end);
%         timing.ms_onset = timing.ms_onset + round(diff * k / ylen);
        
        save(fullfile(save_dir, save_file_as), 'contourmap', 'sonogram', 'emg', 'type', 'timerange', 'freqrange', 'simple', 'notmatfile', 'parameters', 'timing', 'songdata');
    end
end