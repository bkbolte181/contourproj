function recdata=readrecf_old_test(fname)
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
%    outfile:
%     header:

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

    if (length(findstr(upper(fl),'CATCH'))>0)
        pp = findstr(fl,'=');
        recdata.iscatch = str2num(fl(pp(end)+1:end));

    elseif (length(findstr(upper(fl),'CHANS'))>0)
        pp = findstr(fl,'=');
        recdata.nchan = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'ADFREQ'))>0)
        pp = findstr(fl,'=');
        recdata.adfreq = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'SAMPLES'))>0)
        pp = findstr(fl,'=');
        recdata.nsamp = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'T AFTER'))>0)
        pp = findstr(fl,'=');
        recdata.tafter = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'T BEFORE'))>0)
        pp = findstr(fl,'=');
        recdata.tbefore = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'OUTPUT SOUND FILE'))>0)
        pp = findstr(fl,'=');
        recdata.outfile = fl(pp(end)+1:end);
    elseif (length(findstr(upper(fl),'THRESHOLDS'))>0)
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
    elseif (length(findstr(lower(fl),'feedback information'))>0)
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
%            tmpvec = [tmpvec;str2num(fl(1:tpos(1)-1))];
 tmpvec=[];
        end
        recdata.ttimes = tmpvec;
        if (~ischar(fl))
            break;
        end
    elseif (length(findstr(lower(fl),'trigger times'))>0)
        tmpvec = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(fl)<2)
                continue;
            end
            tmpvec = [tmpvec;str2num(fl)];
        end
        recdata.ttimes = tmpvec*1e3;
        if (~ischar(fl))
            break;
        end
    elseif (length(findstr(upper(fl),'FILE CREATED:'))>0)
        % save the next 4 lines and this line as the header
        headr{1} = fl;
        for ijk=1:4
            fl=fgetl(fid);
            headr{ijk+1}=fl;
        end
        recdata.header=headr;
    end
end
fclose(fid);

% disp('HACK - SETTING TBEFORE AND TAFTER TO ZERO')
%     recdata.tbefore=0;
%     recdata.tafter=0;

return;
