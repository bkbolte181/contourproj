% general idea to generate info for compilecontours to use
clear
f.all_dirs{1}.name = 'contours';
mydir = dir;
for i=3:numel(mydir)
    files{i-2} = mydir(i).name;
end
counter = 1;
while ~exist('myfile','var') || isempty(myfile)
    if exist(files{counter}, 'file') && ~isdir(files{counter})
        myfile = load(files{counter});
    else
        files{counter} = [];
        counter = counter + 1;
    end
end
files = files(~cellfun('isempty',files));
f.all_dirs{1}.files = files;
[contours sonos accessory ids] = compilecontours(f);

% Time and freqency ranges for actual data
timescale = accessory.basefile.timerange(2) - accessory.basefile.timerange(1);
freqscale = accessory.basefile.freqrange;

% Seems like a hack; probably a better way to do this but I don't know how
ccounter = 0;
scounter = 0;
for i=1:numel(contours)
    contours{i} = contours{i} * (freqscale(2)-freqscale(1)) / size(myfile.contourmap,1) + freqscale(1);
    if ~ids(i)
        ccounter = ccounter + 1;
        catchdata{ccounter} = contours{i};
    else
        scounter = scounter + 1;
        stimdata{scounter} = contours{i};
        emgdata{scounter} = accessory.emg{i};
    end
end

catchdata = cell2mat(catchdata); % Columns are contours, rows are time points
stimdata = cell2mat(stimdata);
emgdata = cell2mat(emgdata);

conf_int = zeros(size(catchdata,1),1); % Main confidence interval
conf_int2 = zeros(size(catchdata,1),1);

% Calculate confidence interval using 95th quantile
for i=1:size(catchdata,1)
    % Remove zero-valued elements
    mydata = catchdata(i,:)';
    mydata = mydata(mydata~=freqscale(1));
    if numel(mydata) < 5
        mydata = catchdata(i,:)';
    end
    pd = fitdist(mydata,'Normal');
    % Which one of these is better?
    conf_int2(i) = pd.mu + 2*pd.std; % Mean plus 2 standard deviations
    conf_int(i) = quantile(mydata,0.975); % 95th quantile
end

% For viewing quantile vs stddev methods
% plot(conf_int, 'r');
% hold on;
% plot(quantiledata,'g');
% for i=1:length(catchdata(1,:))
%     plot(catchdata(:,i));
%     hold on;
% end

% Figure out where stim data becomes greater than conf_int
for i=1:size(stimdata,2)
    invalid = false; % Marker for invalid contours
    a = emgdata(:,i);
    a = find(a>mean(a));
    b = a(1);
    a = a(end); % Last stimulation point in emg
    c = a;
    mystim = stimdata(:,i);
    scalar = (length(mystim)/length(emgdata(:,i)));
    while mystim(round(a*scalar)) > conf_int(round(a*scalar))
        % Work forwards to find where the line crosses
        a = a - 1;
        if mystim(round(a*scalar)) == freqscale(1)
            % Throw this data point out, bad contour
            invalid = true;
        end
    end
    while mystim(round(c*scalar)) > conf_int2(round(c*scalar))
        % Work forwards to find where the line crosses
        c = c - 1;
    end
    if invalid
        intersections(i) = -1;
        latencies(i) = -1;
        latencies2(i) = -1;
        continue;
    end
    a = a + 0.5;
    c = c + 0.5;
    data{i}.intersections = a * timescale / length(emgdata(:,i));
    data{i}.latencies = (a - b) * timescale / length(emgdata(:,i));
    data{i}.latencies2 = (c - b) * timescale / length(emgdata(:,i));
    
    % Store variables in struct
    data{i}.stimx = linspace(0,timescale,length(mystim));
    data{i}.emgx = linspace(0,timescale,length(emgdata(:,i)));
    data{i}.stim = mystim;
    data{i}.confint = conf_int;
    data{i}.emg = emgdata(:,i);
end

data = data(~cellfun('isempty',data)); % Remvoe empty cells
% clearvars -except data;

latencies = zeros(numel(data),1);
latencies2 = zeros(numel(data),1);
for i = 1:numel(data)
    latencies(i) = data{i}.latencies;
    latencies2(i) = data{i}.latencies2;
end