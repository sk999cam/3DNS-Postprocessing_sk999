function [xprof,yprof]=blade_profile(chi,tLE,N)
% create circular arc blade with 'shape-space' thickness distribution

chi_1 = chi(1)*pi/180;
chi_2 = chi(2)*pi/180;

cam = abs(chi_1-chi_2);
stag = (chi_1 + chi_2)*0.5;

c = linspace(0,1,N);
rad = 0.5/sin(cam*0.5);

thet_ex = pi/2 + chi_1;
thet_in = pi/2 + chi_2;

thet = linspace(thet_ex,thet_in,N);
r = linspace(rad,rad,N);

% add thickness
psi = linspace(0,1,N);

% shape space
S = 1 ;
phi(1) = tLE;
phi = sqrt(psi).*(1-psi).*S + psi*phi(1);

t = 5*tLE*phi;
tTE = t(N);

% add TE
rthet = r.*thet;

for nn=1:3
iTE = rthet-rthet(N) < tTE;
tTE = max(t(iTE));
end

a = tTE - abs(rthet(N)-rthet(iTE));
t(iTE) = sqrt(tTE*tTE - a.*a);

rup = r + t*0.5;
rdn = r - t*0.5;

[xc,yc] = pol2cart(thet,r);
[xu,yu] = pol2cart(thet,rup);
[xd,yd] = pol2cart(thet,rdn);

xprof = [xu(1:N-1) xd(N:-1:1)];
yprof = [yu(1:N-1) yd(N:-1:1)];

% remove repeat points and smooth
N = length(xprof);
for mm=1:1
ii = 2:N-1;
xprof(ii) = (xprof(ii-1) + xprof(ii+1))*0.5;
yprof(ii) = (yprof(ii-1) + yprof(ii+1))*0.5;
% xprof(1) = (xprof(N-1) + xprof(2))*0.5;
% yprof(1) = (yprof(N-1) + yprof(2))*0.5;
% xprof(N) = xprof(1);
% yprof(N) = yprof(1);
end


return