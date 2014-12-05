function tf = isvoiced(spec_slice, cont_slice)
% Determine if a particular slice of a spectrogram is voiced or unvoiced by
% analyzing if it is concentrated around a particular contour

tf = 0;
cloc = find(cont_slice);
if cloc
    % Find peaks of this slice of sonogram
    [pts, loc] = findpeaks(spec_slice);
    
    % Find which peaks are close to contour
    m = find(abs(loc-cloc(1))<=3);
    
    % Value of points close to contour
    pt = pts(m);
    
    if numel(pt) == 1
        pt = pt(1);
        if numel(pts) == 1
            % Only one peak, at contour
            tf = 1;
        elseif pt > mean(pts)
            tf = 1;
        end
    end
end