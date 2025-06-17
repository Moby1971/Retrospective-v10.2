function [xcList, ycList] = polyclip(x, y, xmin, xmax, ymin, ymax)
%#ok<*AGROW>

xcList = {};
ycList = {};

for k = 1:numel(x)-1
    [accept, x1, y1, x2, y2] = liangBarsky(x(k), y(k), x(k+1), y(k+1), xmin, xmax, ymin, ymax);
    if accept
        xcList{end+1} = [x1, x2];
        ycList{end+1} = [y1, y2];
    end
end

end

function [accept, x1, y1, x2, y2] = liangBarsky(x0, y0, x1, y1, xmin, xmax, ymin, ymax)
dx = x1 - x0; dy = y1 - y0;
p = [-dx dx -dy dy];
q = [x0 - xmin, xmax - x0, y0 - ymin, ymax - y0];
u1 = 0; u2 = 1;
accept = false;
for i = 1:4
    if p(i) == 0
        if q(i) < 0
            return; % Line is parallel and outside
        end
    else
        t = q(i)/p(i);
        if p(i) < 0
            u1 = max(u1, t);
        else
            u2 = min(u2, t);
        end
    end
end

if u1 > u2
    return;
end

accept = true;
x1 = x0 + u1*dx; y1 = y0 + u1*dy;
x2 = x0 + u2*dx; y2 = y0 + u2*dy;

end
