function [ contourmaps, sonograms, accessory, ids ] = compilecontours( f )
% Compiles contours into one image, based on directory information passed
% as a struct

if ischar(f)
    f = load(f);
end

global sum_ps;

for i=1:numel(f.all_dirs)
    mydir = f.all_dirs{i};
    mydirname = f.all_dirs{i}.name;
    counter = 0;
    for j=1:numel(mydir.files)
        disp([j numel(mydir.files)]);
        counter = counter + 1;
        myfile = load(mydir.files{j});
        if isfield(myfile, 'sum_ps')
            sum_ps = myfile.sum_ps;
        else
            sum_ps = 0;
        end
        accessory.sum_ps = sum_ps;
        if j==1 % Do this on the first iteration
            accessory = struct;
            accessory.basefile = myfile;
            contourmaps = cell(numel(mydir.files),1);
            sonograms = cell(numel(mydir.files),1);
            accessory.emg = cell(numel(mydir.files),1);
            ids = zeros(numel(mydir.files),1);
            if sum_ps
                basesize = length(myfile.contourmap);
            else
                basesize = size(myfile.contourmap,2); % Size to which all images will be resampled
            end
            baseemgsize = length(myfile.emg);
        end
%         freqrange = g.freqrange;
%         timerange = [g.in g.out];
        if sum_ps
            contourmaps{counter} = changeToSize(myfile.contourmap', basesize);
        else
            myfile.contourmap = myfile.contourmap ./ max(max(abs(myfile.contourmap)));
            contourmaps{counter} = changeToSize(average_img(myfile.contourmap),basesize);
        end
        myfile.sonogram = myfile.sonogram ./ max(max(abs(myfile.sonogram)));
        sonograms{counter} = changeToSize(average_img(myfile.sonogram),basesize);
        accessory.emg{counter} = changeToSize(myfile.emg,baseemgsize);
        if strcmp(myfile.type,'stim')
            ids(counter) = 1; % mark all stim trials as ones
        end
        accessory.sonogram{counter} = myfile.sonogram;
        accessory.contourmap{counter} = myfile.contourmap;
        accessory.timerange{counter} = myfile.timerange;
        [path, name, ext] = fileparts(myfile.notmatfile.fname);
        accessory.filenames{counter} = name;
    end
end

end

function dataOut = changeToSize(dataIn, targetSize)
dataOut = resample(dataIn, targetSize, length(dataIn));
% if length(dataIn) < targetSize
%     dataOut = padarray(dataIn,[targetSize-length(dataIn) 0],'post');
% else
%     dataOut = dataIn(1:targetSize,:);
% end

end