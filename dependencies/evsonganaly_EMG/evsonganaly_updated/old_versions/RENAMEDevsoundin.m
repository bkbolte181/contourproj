function [dat,fs]=RENAMEDevsoundin(fullfname,filetype)
%[dat,Fs]=ReadDataFile(fullfname,filetype);
%
% Reads data from soundfile into matlab vector
%
% INPUTS:
% 
%   fullname - is the full file name with path and extension
%   filetype - is not neecessary if the extension is .cbin, .wav, or .ebin
%   
% OUTPUTS:
%   dat - is a matrix, with all the colums of data
%         wav files will only have one column
%   fs  - The sampling rate if available from the file specified in SOUNDFILE.
%           If the sampling rate is not available, then Fs = -1 is returned



if (exist('filetype'))
    if (strcmp(filetype(1:3),'obs'))
        ext = '.cbin';
        
        
    else
    
        [pth,fn,ext]=fileparts(fullname);
        ext = lower(ext);
end

    if (strcmp(ext,'.cbin'))
        [dat,fs] = ReadObsData(soundfile_full,chan_spec);
    elseif (strcmp(ext,'.ebin'))
        [dat,fs] = ReadEvTaf(soundfile_full,chan_spec);
    elseif (strcmp(ext,'.wav'))
        [dat,fs] = wavread(fullname);
    end



return;