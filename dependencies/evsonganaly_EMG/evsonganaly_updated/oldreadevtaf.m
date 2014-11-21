function [data,fs]=readevtaf(fname,chanspec)
%[data,fs]=readevtaf(fname,chanspec);
%

if (exist(fname,'file'))
	fid=fopen(fname,'r','b');
	data=fread(fid,inf,'float');
	fclose(fid);

	% load evtaf defaults
	fs = 32000.0;
	Nchan = 2;
	Ndata = length(data)/Nchan;

	pos = strfind(fname,'.ebin');
	if (length(pos)==0)
		recfile = [fname,'.rec'];
	else
		recfile = [fname(1:pos(end)),'rec'];
	end
	if (~exist(recfile,'file'))
		warning(['Could not file rec file: ',recfile,...
		         '\n Assuming standard fs and 1 data channel']);
	else
		frec = fopen(recfile,'r');
		recln=fscanf(frec,'%c');
		fclose(frec);
		rm=recln;
		while (length(rm)>0)
			[ln,rm]=strtok(rm,13); % cr = delimiter
			pos = strfind(ln,'=');
			if (findstr(upper(ln),'ADFREQ'))
				if (length(pos)==0)
				   pos=strfind(upper(ln),'ADFREQ')+5;
			        end
				fs=str2num(ln((pos(1)+1):end));
			elseif (findstr(upper(ln),'CHANS'))
				if (length(pos)==0)
				   pos=strfind(upper(ln),'CHANS')+4;
			        end
				Nchan=str2num(ln((pos(1)+1):end));
			elseif (findstr(upper(ln),'SAMPLES'))
				if (length(pos)==0)
				   pos=strfind(upper(ln),'SAMPLES')+6;
			        end
				Ndata=str2num(ln((pos(1)+1):end));
			end
        end
    end
    
    if (strcmp(chanspec(end),'r'))
        whichchan = Nchan - str2num(chanspec(1:end-1));
    else
        whichchan = str2num(chanspec);
    end

	if (Ndata*Nchan~=length(data))
		warning(['Data size does not match REC file\n',...
			         'Returning EMPTY MATRIX!']);
                data = [];
		fs = 0.0;
	end

	if (exist('whichchan'))
		data = data(whichchan:Nchan:end);
	else
		datatmp = zeros([Ndata,Nchan]);
		for ii = 1:Nchan
			datatmp(:,ii) = data(ii:Nchan:end);
		end
		data=datatmp;
	end
else
	warning(['Could not find file: ',fname,...
	         '\n skipping it']);
end
return;
