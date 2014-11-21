function [pts_in_order] = collect_points_to_contour(contour, pt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    
% Point nearest "pt"
k = dsearchn(contour,pt);
my_pt = contour(k,:);

% Make temporary array, so that we can keep track of "visited" points
temp = contour;

% Initialize the array that will be returned
pts_in_order = [];

% Keep track of the distance from point to point
distances = [];

% Do this until we've hit every point
while ~isempty(temp)
    % Find closest point
    [k,d] = dsearchn(temp,my_pt);

    if numel(distances) > 10 && d > mean(distances) + 3 * std(distances)
        % If the distance to this point exceeds a certain threshold,
        % end the loop
        break;
    end

    % Store the location of that point
    my_pt = temp(k,:);

    % Store the distance of that point
    distances = [distances; d];

    % Add it to the next spot in our returned array
    pts_in_order = [pts_in_order; my_pt];

    % Delete that point
    temp(k,:) = [];
end

% Return our new, ordered array
pts_in_order;

end

