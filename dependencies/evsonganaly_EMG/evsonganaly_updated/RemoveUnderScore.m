%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutString=RemoveUnderScore(InString)
% replaces all _ with \_ for proper display

TmpStr=InString;
pos = findstr(InString,'_');
for ind=1:length(pos)
    indu = pos(ind);
    if (indu==1)
        TmpStr=['\',TmpStr];
    else
        TmpStr = [TmpStr(1:(indu-1)),'\',TmpStr(indu:end)];
    end
    pos = pos + 1;
end
OutString=TmpStr;
return;