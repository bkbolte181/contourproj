function [ polynomial, extrainfo ] = img_to_cart( f, degree )
% Convert image data to Cartesian coords, then fits a polynomial of degree
% <degree> (default is 3rd degree) and graphs it. Returns the polynomial as
% output. The polynomial is fit through the weighted average of the contour
% points given as input; however, if there is only one point at a
% particular contour sample for a long stretch of time, the data may end up
% looking more digital. Also returns some extra information in the form of
% a struct.

if ~exist('degree','var') || isempty(degree)
    degree = 3;
end

if ischar(f)
    f = load(f); % If the file is passed as a string
end

% Load data and crop pre- and post-syllable data out
img = f.data;
img = img(:,f.mapbounds(1)+1:f.mapbounds(2)-1);

[ylen xlen] = size(img);

x = 1:xlen;
y = zeros(size(xlen));

for i=1:xlen
    if sum(img(:,i)) == 0
        y(i) = 0;
    else
%         % Code for generating weighted image
        weights = img(:,i);
        weights = weights(weights~=0); % remove zeros so that weights and val are same size
        vals = find(img(:,i));
%         weights = ones(size(vals)); % uncomment do disable weighting
        y(i) = sum(vals.*weights)/sum(weights); % weighted average
    end
end

% Remove zero-valued points
x = x(y~=0);
y = y(y~=0);

% Convert data to Hz and ms
y = y * ((f.freqrange(2) - f.freqrange(1)) / ylen) + f.freqrange(1);
x = x * ((f.out - f.in) / xlen);

[p,S] = polyfit(x,y,degree);
polydata = polyval(p,x);
% plot(x,polydata,x,y,'.') % For plotting

polynomial = p;
extrainfo = S;
extrainfo.polydata = polydata;
extrainfo.x = x;
extrainfo.y = y;
extrainfo.xlen = xlen;
extrainfo.ylen = ylen;

sstot = sum((y - mean(y)).^2);
ssres = sum((y - polydata).^2);
extrainfo.rsquared = 1 - (ssres / sstot);

end