function [pts_in_order] = collect_points_to_contour(contour, pt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    
% Point nearest "pt"
k = dsearchn(contour,pt);
my_pt = contour(k,:);

% Make temporary array, so that we can keep track of "visited" points
temp = contour;

% Initialize the array that will be returned
pts_in_order = zeros(size(contour));

% Keep track of the distance from point to point
distances = zeros(size(contour(:,1)));

i = 0;

% Do this until we've hit every point
while ~isempty(temp)
    % Find closest point
    [k,d] = dsearchn(temp,my_pt);
    
    dist_no_zeros = distances(distances~=0);
    
    if numel(dist_no_zeros) > 10 && d > mean(dist_no_zeros) + 3 * std(dist_no_zeros)
        % If the distance to this point exceeds a certain threshold,
        % end the loop
        break;
    end
    
    i = i + 1;

    % Store the location of that point
    my_pt = temp(k,:);

    % Store the distance of that point
    distances(i) = d;

    % Add it to the next spot in our returned array
    pts_in_order(i,:) = my_pt;

    % Delete that point
    temp(k,:) = [];
end

% Remove empty points
pts_in_order(~any(pts_in_order,2),:) = [];

end

