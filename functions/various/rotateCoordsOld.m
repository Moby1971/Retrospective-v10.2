function [xr, yr] = rotateCoords(x, y, imageSize, k)
    k = mod(k, 4); % Only 0,1,2,3 are valid
    xr = x;
    yr = y;
    w = imageSize(2);
    h = imageSize(1);
    for i = 1:k
        % Rotate 90° counter-clockwise around image center
        x0 = xr - w/2;
        y0 = yr - h/2;
        xr = -y0 + w/2;
        yr =  x0 + h/2;
    end
end
