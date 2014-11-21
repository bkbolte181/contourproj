function [r] = cell_to_points(contour)
% CELL_TO_POINTS Convert array of cells to points. An accessory function for
% simple_contours

r = [0,0];
for i=1:numel(contour)
    for j=1:numel(contour{i})
        r = [r; [i, contour{i}(j)]];
    end
end

end

