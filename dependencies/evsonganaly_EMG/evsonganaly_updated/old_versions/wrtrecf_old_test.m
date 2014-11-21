function wrtrecf_old_test(fname,recdata)
% wrtrecf(fname,recdata);
% recdata is a structure with all the rec file fields

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
recf
fid = fopen(recf,'w');
if (isfield(recdata,'header'))
    for ii=1:length(recdata.header)
        fprintf(fid,'%s\n',recdata.header{ii});
    end
    fprintf(fid,'\n');
end

if (isfield(recdata,'adfreq'))
    fprintf(fid,'ADFREQ = %12.7e\n',recdata.adfreq);
end

if (isfield(recdata,'outfile'))
    fprintf(fid,'Output Sound File =%s\n',recdata.outfile);
end

if (isfield(recdata,'nchan'))
    fprintf(fid,'Chans = %d\n',recdata.nchan);
end

if (isfield(recdata,'nsamp'))
    fprintf(fid,'Samples = %d\n',recdata.nsamp);
end

if (isfield(recdata,'iscatch'))
    fprintf(fid,'Catch = %d\n',recdata.iscatch);
end

if (isfield(recdata,'tbefore'))
    fprintf(fid,'T BEFORE = %12.7e\n',recdata.tbefore);
end

if (isfield(recdata,'tafter'))
    fprintf(fid,'T AFTER = %12.7e\n',recdata.tafter);
end

if (isfield(recdata,'thresh'))
    fprintf(fid,'THRESHOLDS = ');
    for ii = 1:length(recdata.thresh)
            fprintf(fid,'\n%15.7e',recdata.thresh(ii));
    end
end
fprintf(fid,'\n');

if (isfield(recdata,'ttimes'))
    %fprintf(fid,'Trigger times = ');
    fprintf(fid,'Feedback information:\n\n');

    for ii = 1:length(recdata.ttimes)
        if (recdata.ttimes(ii)>=0)
            disp(['trying to write ' num2str(round(recdata.ttimes(ii)))])
%            fprintf(fid,'%15.7e \n',recdata.ttimes(ii));
            fprintf(fid,'\n%15.7e msec : FB',recdata.ttimes(ii));
        end
    end
end


fclose(fid);
return;
