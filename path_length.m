function s = path_length(x, y)
    x = reshape(x,1,[]);
    y = reshape(y,1,[]);
    dx = x(2:end) - x(1:end-1);
    dy = y(2:end) - y(1:end-1);
    ds = sqrt(dx.*dx + dy.*dy);
    s = [0 cumsum(ds)];
end