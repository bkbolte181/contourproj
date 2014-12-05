function [contour, sonogram, f, t] = simple_contour_from_syllable(data, SAMPLING, N, OVERLAP)
%SIMPLE_CONTOURS_FROM_SYLLABLE Use simple contours, basically

tScale = 1;

t = -N/2+0.5:N/2-0.5;
tScale=tScale*SAMPLING/1000;

w = exp(-(t/tScale).^2); % Gaussian function
[sonogram f t] = spectrogram(data, w, OVERLAP, N, SAMPLING); % Gaussian windowed spectrogram

sonogram = abs(sonogram);

contour = simple_contours(sonogram, [100, 160]);

end

