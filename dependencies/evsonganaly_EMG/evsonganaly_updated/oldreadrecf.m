function recdata=readrecf(fname)
% recdata=readrecf(fname);
%

pp = findstr(fname,'.rec');
if (length(pp)<1)
	pp2 = findstr(fname,'.');
	if (length(pp2)<1)
		recf = [fname,'.rec'];
	else
		recf = [fname(1:pp2(end)),'rec'];
	end
else
	recf = fname;
end

if (~exist(recf,'file'))
	disp(['Rec file : ',recf,' does not exist']);
	recdata = [];
	return;
end

fid = fopen(recf,'r');
while (1)
	fl = fgetl(fid);
	if (~ischar(fl))
		break;
	end
	if (strcmp(upper(fl(1:6)),'ADFREQ'))
		pp = findstr(fl,'=');
		recdata.adfreq = str2num(fl(pp(end)+1:end));
	elseif (strcmp(upper(fl(1:5)),'CHANS'))
		pp = findstr(fl,'=');
		recdata.nchan = str2num(fl(pp(end)+1:end));
	elseif (strcmp(upper(fl(1:7)),'SAMPLES'))
		pp = findstr(fl,'=');
		recdata.nsamp = str2num(fl(pp(end)+1:end));
	elseif (strcmp(upper(fl(1:5)),'CATCH'))
		pp = findstr(fl,'=');
		recdata.iscatch = str2num(fl(pp(end)+1:end));
	elseif (strcmp(upper(fl(1:7)),'TRIGGER'))
		tmpvec = [];
		while (1)
			fl = fgetl(fid);
			if (~ischar(fl))
				break;
			end
			if (length(str2num(fl(1)))==0)
				break;
			end
			tmpvec = [tmpvec;str2num(fl)];
		end
		recdata.ttimes = tmpvec;
		if (~ischar(fl))
			break;
		end
	end
end
fclose(fid);
return;
