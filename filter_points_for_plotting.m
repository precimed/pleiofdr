function [x, y] = filter_points_for_plotting(x, y, image_size)
%% FILTER_POINTS_FOR_PLOTTING removes points from (x, y) vector
%                             based on the resolution of target plot
%     x                 : vector of x coordinates
%     y                 : vector of y coordinates
%     resolution        : a vector with two values indicating resolution
%                         of the target plot (for example [640, 480] or
%                         [1920, 1080]
    if (isempty(image_size) || any(isnan(image_size)))
        return % skip filtering if resultion is not specific
    end

    if length(x) ~= length(y)
        error('x and y length mismatch')
    end

    % Remove NaNs
    nans = isnan(x) | isnan(y);
    x = x(~nans);
    y = y(~nans);
    index = false(length(x), 1);

    minX = min(x);
    maxX = max(x);
    minY = min(y);
    maxY = max(y);

    bucket_x = (maxX - minX) / image_size(1);
    bucket_y = (maxY - minY) / image_size(2);

    indX = 1 + int32((x - minX) / bucket_x);
    indY = 1 + int32((y - minY) / bucket_y);

    bins = false(image_size(1) + 1, image_size(2) + 1);
    for i = 1:length(indX)
        if (bins(indX(i), indY(i)))
            continue;
        end

        bins(indX(i), indY(i)) = true;
        index(i) = true;
    end

    x = x(index);
    y = y(index);
end
