function [data,fs]=readevtaf(fname,chanspec)
%[data,fs]=readevtaf(fname,chanspec);
% chanel spec starts at 1 or subtracts from Nchan '0r'  = Nchan - 0;

if (~exist('chanspec'))
    chanspec = '0r';
end

if (exist(fname,'file'))
	fid=fopen(fname,'r','b');
	data=fread(fid,inf,'float');
	fclose(fid);

	% load evtaf defaults
	fs = 32000.000;
	Nchan = 2;
	Ndata = length(data)/Nchan;

	pos = strfind(fname,'.ebin');
	if (length(pos)==0)
		recfile = [fname,'.rec'];
	else
		recfile = [fname(1:pos(end)),'rec'];
	end
	if (~exist(recfile,'file'))
		warning(['Could not find rec file: ',recfile,...
		         ' - Assuming standard fs = 32kHz and 2 data channel']);
    else
        rdata=readrecf(fname);
        fs=rdata.adfreq;
        Nchan=rdata.nchan;
        Ndata=rdata.nsamp;
    end
    
    if (strcmp(chanspec(end),'r'))
        whichchan = Nchan - str2num(chanspec(1:end-1));
        if ((whichchan < 1)|(whichchan>Nchan))
            warning(['Chan spec not right using ''0r''']);
            whichchan = Nchan;
        end
    else
        whichchan = str2num(chanspec);
    end

	if (Ndata*Nchan~=length(data))
        warning(['Data size does not match REC file - ',...
                  'Returning vector bases only on Nchan!']);             
        Ndata = length(data)/Nchan;
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
         data = [];
         fs=0.0;
end
return;
