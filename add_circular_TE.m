function [x, y] = add_circular_TE(xprof, yprof, TEpoints)
    
    if nargin < 3
        TEpoints = 101;
    end
    
    if ~ispolycw([xprof(:); xprof(1)], [yprof(:); yprof(1)])
        xprof = flip(xprof);
        yprof = flip(yprof);
    end
    
    g1 = [(xprof(1)-xprof(2)) (yprof(1)-yprof(2)) 0]';
    g2 = [(xprof(end)-xprof(end-1)) (yprof(end)-yprof(end-1)) 0]';
    
    r1 = [xprof(1); yprof(1); 0];
    r2 = [xprof(end); yprof(end); 0];
    
    d1 = 0.5*dot(g1,(r2-r1))/dot(g1,g1);
    d2 = 0.5*dot(g2,(r1-r2))/dot(g2,g2);
    p1 = r1+d1*g1;
    p2 = r2+d2*g2;
    c = 0.5*(r1+r2);
    r = norm(c-p1);
    
    th1 = angle(p1(1)-c(1)+1i*(p1(2)-c(2)));
    th2 = angle(p2(1)-c(1)+1i*(p2(2)-c(2)));
    
    th = linspace(th2, th1, TEpoints)';
    
    xTE = c(1)+r*cos(th);
    yTE = c(2)+r*sin(th);
    
    [~, iTE] = max(xTE);
    x = [xTE(iTE:end); xprof(:); xTE(1:iTE)];
    y = [yTE(iTE:end); yprof(:); yTE(1:iTE)];

end
