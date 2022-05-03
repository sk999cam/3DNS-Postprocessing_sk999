clear
close all

Rec = 200e3;
Dxiplus = sqrt(80)*Rec^(-13/14)

dy0 = 1e-4;
dy = dy0;
y = dy0;
while y < 0.05
    dy = dy*1.03;
    y(end+1) = y(end) + dy;
end

ny = 96;
Lo = 0.05;
y1 = 0.0001;
r0 = 1.025;
r1 = 1.03;

ys = y1*cumsum((r0*ones(1,ny)).^(0:(ny-1)));
yn = ys(end)