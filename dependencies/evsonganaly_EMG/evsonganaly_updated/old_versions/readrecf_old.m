function recdata=readrecf_old(fname)
% recdata=readrecf(fname);
%
%     adfreq: 32000
%      nchan: 2
%      nsamp: 199684
%    iscatch: 0
%     ttimes: []
%    tbefore:
%     tafter:
%     thresh:
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
    %disp(['Rec file : ',recf,' does not exist']);
    recdata = [];
    return;
end

fid = fopen(recf,'r');
flsv=[];
while (1)
    if (length(flsv)==0)
        fl = fgetl(fid);
    else
        fl = flsv;
        flsv=[];
    end
    if (~ischar(fl))
        break;
    end
    
    if length(fl)==0
        continue;
    end

    if (strcmp(upper(fl(1:5)),'CATCH'))
        pp = findstr(fl,'=');
        recdata.iscatch = str2num(fl(pp(end)+1:end));

    elseif (strcmp(upper(fl(1:5)),'CHANS'))
        pp = findstr(fl,'=');
        recdata.nchan = str2num(fl(pp(end)+1:end));
    elseif (strcmp(upper(fl(1:6)),'ADFREQ'))
        pp = findstr(fl,'=');
        recdata.adfreq = str2num(fl(pp(end)+1:end));
    elseif (strcmp(upper(fl(1:7)),'SAMPLES'))
        pp = findstr(fl,'=');
        recdata.nsamp = str2num(fl(pp(end)+1:end));
    elseif (strcmp(upper(fl(1:7)),'T AFTER'))
        pp = findstr(fl,'=');
        recdata.tafter = str2num(fl(pp(end)+1:end));
    elseif (strcmp(upper(fl(1:8)),'T BEFORE'))
        pp = findstr(fl,'=');
        recdata.tbefore = str2num(fl(pp(end)+1:end));
    elseif (strcmp(upper(fl(1:10)),'THRESHOLDS'))
        tmpvec = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(str2num(fl))==0)
                flsv = fl;
                break;
            end
            tmpvec = [tmpvec;str2num(fl)];
        end
        recdata.thresh = tmpvec;
        if (~ischar(fl))
            break;
        end
    elseif (strcmp(lower(fl(1:20)),'feedback information'))
        tmpvec = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(fl)<2)
                continue;
            end
            tpos = findstr(fl,'msec');
            tmpvec = [tmpvec;str2num(fl(1:tpos(1)-1))];
        end
        recdata.ttimes = tmpvec;
        if (~ischar(fl))
            break;
        end
    end
end
fclose(fid);
return;
