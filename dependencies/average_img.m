function avgcont = average_img( imgdata, freqrange, timerange, threshold )
% Returns the average contour of some image data

[y x] = size(imgdata); % y: vert length x: hori length

avgcont = zeros(length(imgdata),1); % Average contour line
zaxis = zeros(length(imgdata),1); % Total power along contour line - not returned

imgdata(imgdata<0) = 0; % Negative values come from resampling and are negligible

for i=1:length(imgdata(1,:))
    if any(imgdata(:,i))
        weights = imgdata(:,i);
        vals = find(imgdata(:,i));
        vals = getconnected(vals);
        weights = weights(vals); % remove zeros so that weights and val are same size
        weights = weights / sum(weights);
        stddev(i) = sqrt(var(vals,weights));
        avgcont(i) = sum(vals.*weights);
        zaxis(i) = sum(weights);
    end
end

if nargin >= 2
    % Adjust according to freqrance
    avgcont = avgcont * (freqrange(2) - freqrange(1)) / y + freqrange(1);
    stddev = stddev * (freqrange(2) - freqrange(1)) / y;
end

if nargin >= 3
    % Resample according to freq and time ranges
    avgcont = resample(avgcont, diff(timerange), length(avgcont));
    stddev = resample(stddev, diff(timerange), length(stddev));
    zaxis = resample(zaxis,diff(timerange), length(zaxis));
end

if nargin >= 4
    avgcont = avgcont(zaxis > threshold);
    stddev = stddev(zaxis > threshold);
end

    function res = diff(tea)
        % Returns absolute difference between two elements in a two-element
        % array
        res = round(abs(tea(2) - tea(1)));
    end

    function res = getconnected(v)
        for j = 1:numel(v)-1
            if v(j+1) ~= v(j)+1
                v = v(1:j);
                break;
            end
        end
        res = v;
    end

end