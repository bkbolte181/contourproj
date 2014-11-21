function [smooth,spec,t,f]=SmoothData(rawsong,Fs,DOFILT,nfft,olap,sm_win,F_low,F_High)
% [smooth,spec,t,f]=evsmooth(rawsong,Fs,DOFILT,nfft,olap,sm_win,F_low,F_High,DOFILT);
% returns the smoothed waveform/envelope + the spectrum
%


if (~exist('F_low'))
	F_low  = 500.0;
end
if (~exist('F_high'))
	F_high = 10000.0;
end
if (~exist('nfft'))
	nfft = 512;
end
if (~exist('olap'))
	olap = 0.8;
end
if (~exist('sm_win'))
	sm_win = 2.0;%ms
end

%if(~exist('SPTH'))
%    SPTH=0.01;
%end

if (~exist('DOFILT'))
    DOFILT=1;
end

filter_type = 'hanningfir';

if (DOFILT==1)
    %filtsong=bandpass(rawsong,Fs,F_low,F_high,filter_type);
%disp('SmoothData.m called - now using filtfilt')
filtsong=bandpass_filtfilt(rawsong,Fs,F_low,F_high,filter_type);
else
    filtsong=rawsong;
end
    
if (nargout > 1)
    spect_win = nfft;
    if 1
disp('m_files\evsonganaly\SmoothData.m - allowing overlap for nice-looking spect - could overrun memory')
    noverlap = floor(olap*nfft);
else
disp('hack in m_files\evsonganaly\SmoothData.m - noverlap set to zero')
    noverlap=0;
end
[spec,f,t] = specgram(filtsong,nfft,Fs,spect_win,noverlap);
    spec = abs(spec);
    %p = find(abs(spec)<=SPTH);
    %spec(p) = SPTH;
end

squared_song = filtsong.^2;

%smooth the rectified song
len = round(Fs*sm_win/1000);
h   = ones(1,len)/len;
smooth = conv(h, squared_song);
offset = round((length(smooth)-length(filtsong))/2);
smooth=smooth(1+offset:length(filtsong)+offset);
return;
