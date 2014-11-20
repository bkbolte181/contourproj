function [contourmap, varargout] = zf_contour_func(s,SAMPLING,N,OVERLAP,Threshold,sigmas,FREQLIM,angles)

%% Overview
% This program computes a contour map of some data <s>, as decribed in
% "Sparse Contour Representations of Sound", Lim, Cunningham and Gardner
% 2012.

%% Variables
% <s>           Signal data
%                   e.g. s = wavread('finch_sound_1')
% <SAMPLING>    Sampling frequency in Hz
%                   e.g. SAMPLING = 44100 for wav
%                   SAMPLING = 32000 for .cbins I think
% <N>           Size of window when computing FFT for each segment of the STFT
% <OVERLAP>     Overlap between windows when computing STFT
%                   Each "hop" in STFT is N - OVERLAP samples
% <Threshold>   Percentile specificity
%                   With only 1 value: Specificity of area
%                   With 3 values: Specificity of [area, mean intensity,
%                   intensity variation] in that order
% <sigmas>      Range of time scales for computing each FFT.
%                   sigma is the standard deviation on the Gaussian window
%                   function when computing the FFT.
%                   e.g. sigmas = [1.3 1.6 1.9]
% <FREQLIM>     Frequency limit, upper range of plot in kHz
%                   e.g. FREQLIM = 12 upper lim is 12000 Hz
% <angles>      Angles along which data converges to contours

%% Recommendations for Speed / Memory / Quality
% Smaller sounds are obviously faster to process.
% If N - OVERLAP is large (e.g. N is large and OVERLAP is small) then
%   specgram will run much faster, but be less detailed.
% Having fewer values of sigma greatly increases speed. If possible, try to
%   find a value of sigma that works best and use it.
% Increasing ARThreshold reduces the number of contours in the final image,
%   reducing memory usage

%% Example

% s = wavread('your_syllable.wav');
% [contourmap, sonogram] = zf_contour_func(s, 44100, 1024, 1010, 99, [1.3 1.6
%                                     1.9], 12, 'Area');

%% Parameters

% Specify default parameters
if ~exist('N', 'var') || isempty(N)
    N = 256;
end
if ~exist('OVERLAP','var') || isempty(OVERLAP)
    OVERLAP = max(N - 256, 0);
end
if ~exist('Threshold','var') || isempty(Threshold)
    ARThreshold = 99;
    POThreshold = 99;
    VAThreshold = 99;
elseif length(Threshold) == 1
    ARThreshold = Threshold;
    POThreshold = 0;
    VAThreshold = 0;
elseif length(Threshold) == 3
    ARThreshold = Threshold(1);
    POThreshold = Threshold(2);
    VAThreshold = Threshold(3);
else
    ARThreshold = Threshold(1);
    POThreshold = Threshold(1);
    VAThreshold = Threshold(1);
end
if ~exist('sigmas','var') || isempty(sigmas)
    sigmas = [1.5 1.8 2.1];
end
if ~exist('FREQLIM','var') || isempty(FREQLIM)
    FREQLIM = 12;
end
if ~exist('param','var') || isempty(param)
    param = 'Area';
end
if ~exist('angles','var') || isempty(angles)
    angles = pi/8:pi/8:pi; % Compute contours in all these angles. Recommended to leave this alone.
end

Nshift = N - OVERLAP; % This is the "hop" between each FFT in the STFT
sonogramindex=ceil(numel(sigmas)/2); % Time scale to use for the sonogram (e.g. which angle in <angles>)

NtScale = length(sigmas);
Nangle = length(angles);

global sonogramfinal

%% Bulk of algorithm

for sigmacount=1:NtScale
    tScale = sigmas(sigmacount);
    t = -N/2+0.5:N/2-0.5; % Window length for Gaussian function
    tScale=tScale*SAMPLING/1000; % Convert from ms to samples
    w = exp(-(t/tScale).^2); % Gaussian function
    dw = (w).*(t/(-0.5 * tScale^2)); % Derviative of Gaussian function
    
    [q freqs times] = spectrogram(s, w, OVERLAP, N, SAMPLING); % Gaussian windowed spectrogram
    q2 = spectrogram(s, dw, OVERLAP, N, SAMPLING); % Deriv Gaussian windowed spectrogram
    
    px = floor((N/2)*1000*FREQLIM*2/SAMPLING);
    q = q(max(px(1),1):min(px(2),length(q)),:);
    q2 = q2(max(px(1),1):min(px(2),length(q2)),:);
    
    dx = (q2./q)/(2*pi); % Displacement according to remapping algorithm - see paper
    sonogram = abs(q);
    if sigmacount == sonogramindex,
        sonogramfinal = sonogram; % Keep this sonogram for later
    end
    
    for angle_variable=1:Nangle
        theta = angles(angle_variable); % See paper for logic behind this
        
        m = -1*(imag(dx*exp(1j*theta))<0) + (imag(dx*exp(1j*theta))>0);
        [gx gy] = gradient(m); % Computes gradient of m, e.g. partial derivatives wrt freq and time
        BW = ((gx * (-1 * cos(theta+pi/2)) + gy * sin(theta + pi/2))>0);
        
        CC = bwconncomp(BW); % bwconncomp gets each contour connected component of BW
        
        % cc_pix = regionprops(CC,sonogram,'Area','MeanIntensity','PixelValues');
        cc_pix = regionprops(CC,sonogram,'Area','MeanIntensity','PixelValues');
        areav = zeros(size(cc_pix));
        powerv = zeros(size(cc_pix));
        varv = zeros(size(cc_pix));
        
        for i = 1:length(cc_pix)
            areav(i) = cc_pix(i).Area; % Area taken up by contour
            powerv(i) = cc_pix(i).MeanIntensity; % Mean intensity of pixel values
        end
        
        % Return the indices of the longest/strongest contours
        a = find(areav>=prctile(areav,ARThreshold));
        b = find(powerv>=prctile(powerv,POThreshold));
        
        for i = 1:length(cc_pix)
            if any(i==a) && any(i==b)
                varv(i) = var(cc_pix(i).PixelValues);
                % varv(i) = cc_pix(i).MinIntensity - cc_pix(i).MaxIntensity; % Another approach
            end
        end
        
        c = find(varv);
        c = c(varv(c)<=prctile(varv(c),100-VAThreshold));
        
        tempv = zeros(size(dx));
        
        for i = 1:length(c)
            ind = CC.PixelIdxList(c(i));
            tempv(ind{1}) = 1;
        end
        
        BWallk{sigmacount}{angle_variable}=sparse(tempv);
    end
end
    
% Superimpose all contours from all angles and time scales to produce a
% single image. This final image contains no amplitude information.

consensus = zeros(size(BWallk{1}{1}));
consensuspower=consensus;

for sigmacount=1:NtScale-1,
    temp=zeros(size(consensus));
    for angle_variable=1:Nangle,
        if angle_variable==1
            cv = BWallk{sigmacount}{1} + BWallk{sigmacount+1}{angle_variable} + BWallk{sigmacount}{Nangle};
        else
            cv = BWallk{sigmacount}{angle_variable} + BWallk{sigmacount+1}{angle_variable} + BWallk{sigmacount}{angle_variable-1};
        end
        consensus = consensus + (cv>0);
        temp = temp + consensus;
    end
    consensuspower = temp.*sonogram;
end

contourmap = consensuspower;
if nargout == 2
    varargout{1} = sonogramfinal;
elseif nargout == 4
    varargout{1} = sonogramfinal;
    varargout{2} = freqs;
    varargout{3} = times;
end
