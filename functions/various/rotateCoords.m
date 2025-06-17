function [xr, yr] = rotateCoords(x, y, k)
    % Rotate (x, y) counterclockwise by 90*k degrees around the origin
    k = mod(k, 4); % Only 0,1,2,3 valid
    switch k
        case 0
            xr = x;
            yr = y;
        case 1
            xr = -y;
            yr = x;
        case 2
            xr = -x;
            yr = -y;
        case 3
            xr = y;
            yr = -x;
    end
end
