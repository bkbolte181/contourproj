function [onsets, offsets]=evsegment(smooth, Fs, min_int, min_dur, threshold)
% [ons,offs]=evsegment(smooth,Fs,min_int,min_dur,threshold);
% segment takes smoothed filtered song and returns vectors of note
% onsets and offsets values are in ms 

h=[1 -1];

%threshold input

notetimes=double(smooth>threshold);
%a=find(notetimes);
%a(1)

%extract index values for note onsets and offsets
trans=conv(h,notetimes);
t_onsets  = find(trans>0);
t_offsets = find(trans<0);

onsets = t_onsets;offsets=t_offsets;
if ((length(onsets)<1)|(length(offsets)<1))
	onsets=[];offsets=[];
	return;
end

if (length(t_onsets) ~= length(t_offsets))
	disp('number of note onsets and offsets do not match')
else         
	% use first threshold (set low on purpose) to find general areas of 
	% activity in song then use a percentage of the local smooth value to
	% get the actual on and offset times
	%for ii = 1:length(t_onsets)
	%	tmpdat=smooth(t_onsets(ii):t_offsets(ii));
	%	mxv = max(tmpdat);
	%	tmpdat2=tmpdat>0.05*mxv;
	%	ppos = find(tmpdat2>0);
	%	onsets(ii)  = ppos(1)  -1 + t_onsets(ii);
	%	offsets(ii) = ppos(end)-1 + t_onsets(ii);
	%end

	%eliminate short intervals
	temp_int=(onsets(2:length(onsets))-offsets(1:length(offsets)-1))*1000/Fs;
	real_ints=temp_int>min_int;
	onsets=[onsets(1); nonzeros(onsets(2:length(onsets)).*real_ints)];
	offsets=[nonzeros(offsets(1:length(offsets)-1).*real_ints); offsets(length(offsets))];
    
	%eliminate short notes
	temp_dur=(offsets-onsets)*1000/Fs;
	real_durs=temp_dur>min_dur;
	onsets=[nonzeros((onsets).*real_durs)];
	offsets=[nonzeros((offsets).*real_durs)];
    
	%convert to ms: peculiarities here are to prevent rounding problem
	% if t_ons is simply replaced with onsets, everything gets rounded
	onsets = onsets*1000/Fs;
	offsets = offsets*1000/Fs;
end

