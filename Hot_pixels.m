% Author: James Dixon, 2021
% Set hot pixels to average of their neighbours

% Image and hot pixel indices taken as input. Output image has hot pixels
% set to average of their neighbours.

function[image_out] = Hot_pixels(image,hot_inds_row,hot_inds_col)

if isempty(hot_inds_row) % If no hot pixels, don't change image.
    image_out = image;
    
else
    % Set image within 2 rows/cols frame of zeros
    image_adj = zeros(size(image,1)+2,size(image,2)+2);
    image_adj(2:size(image,1)+1,2:size(image,2)+1) = image;

    for k=1:length(hot_inds_row)

        den = 8; % denominator for inner pixel averaging

        % Lower den where inds are at the edge of the image
        if hot_inds_row(k)==1 || hot_inds_row(k)==size(image,1)
            den = den - 3;
            if hot_inds_col(k)==1 || hot_inds_col(k)==size(image,2)
                den = den - 2;
            end
        elseif hot_inds_col(k)==1 || hot_inds_col(k)==size(image,2)
            den = den - 3;
        end

        % Increase inds so they correspond to correct pixels in image_adj
        row = hot_inds_row(k)+1;
        col = hot_inds_col(k)+1;

        % Set hot pixel to average value of neighbouring pixels
        image_adj(row,col) = (image_adj(row-1,col) + image_adj(row-1,col-1)...
            + image_adj(row,col-1) + image_adj(row+1,col-1) + image_adj(row+1,col)...
            + image_adj(row+1,col+1) + image_adj(row,col+1) ...
            + image_adj(row-1,col+1))/den;
    end

    % Set image to image with adjusted pixels
    image_out = image_adj(2:size(image,1)+1,2:size(image,2)+1);
end
end