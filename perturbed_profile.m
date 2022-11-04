function [xnew, ynew] = perturbed_profile(xprof,yprof, amp, x1, x2, xtrip, r, theta, tmid)

    switch_back = false;
    if yprof(2) < 0
        xprof = flip(xprof);
        yprof = flip(yprof);
        switch_back = true;
    end

    [~, iLE] = min(xprof);
    x = xprof(1:iLE);
    y = yprof(1:iLE);
    s = curve_length(x,y);

    np = length(x);
    
    n(1,:) = [1 0];
    n(np,:) = [-1 0];
    
    for i=2:np-1
        ds = [xprof(i+1)-xprof(i-1) yprof(i+1)-yprof(i-1)];
        ds = ds/norm(ds);
        n(i,:) = [ds(2) -ds(1)];
    end

    [~,i1] = min(abs(x - x1));
    x1 = x(i1);
    y1 = y(i1);
    [~,i2] = min(abs(x - x2));
    x2 = x(i2);
    y2 = y(i2);
    [~, i3] = min(sqrt((xprof-x2).^2 + (yprof + y2).^2));
    [~, i4] = min(sqrt((xprof-x1).^2 + (yprof + y1).^2));
    [~,iTrip] = min(abs(x - xtrip));
    yTrip = y(iTrip);

    [sTrip, delta] = trip_bump(s(i1), s(iTrip), s(i2), r, theta, tmid, s(end));
    sTripNew = linspace(s(i1),s(i2), 501);
    xTrip = interp1(s(i1:i2),x(i1:i2),sTripNew);
    yTrip = interp1(s(i1:i2),y(i1:i2),sTripNew);
    deltaNew = amp*interp1(sTrip,delta,sTripNew,'spline');

    nTrip = length(xTrip);
    
    nNew(1,:) = n(i1,:);
    nNew(nTrip,:) = n(i2,:);
    
    for i=2:nTrip-1
        ds = [xTrip(i+1)-xTrip(i-1) yTrip(i+1)-yTrip(i-1)];
        ds = ds/norm(ds);
        nNew(i,:) = [ds(2) -ds(1)];
    end

    xTrip = xTrip + deltaNew.*nNew(:,1)';
    yTrip = yTrip + deltaNew.*nNew(:,2)';

    xnew = [xprof(1:i1-1) xTrip xprof(i2+1:i3-1) xTrip(end:-1:1) xprof(i4+1:end)];
    ynew = [yprof(1:i1-1) yTrip yprof(i2+1:i3-1) -yTrip(end:-1:1) yprof(i4+1:end)];
    
    if switch_back
        xnew = flip(xnew);
        ynew = flip(ynew);
        xprof = flip(xprof);
        yprof = flip(yprof);
        
    end

end




