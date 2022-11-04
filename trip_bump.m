function [s, deltaHat] = trip_bump(s1,sTrip,s2,r,theta,tmid,sTE)

if nargin < 7
end
%s = linspace(s1,s2,101);

t = 0:0.01:1;
x1 = r*cosd(theta);
y1 = r*sind(theta);

sTrip
sTripHat = (sTrip-s1)/(s2-s1)
sTE

x2 = (1-6*(1-tmid)^2*tmid*x1 - 2*tmid^3)/(6*(1-tmid)*tmid^2);
y2 = (sTripHat-3*(1-tmid)^2*tmid*y1 - tmid^3)/(3*(1-tmid)*tmid^2);


%y1 = max(y1,(8*sTripHat-4)/3);
%y2 = (8*sTripHat-1)/3 - y1;

x = 3*x1.*(1-t).^2.*t + 3*x2.*(1-t).*t.^2 + t.^3;
sHat = 3*y1.*(1-t).^2.*t + 3*y2.*(1-t).*t.^2 + t.^3;
%flip(sHat)
s = s1 + sHat*(s2-s1);
nargin
if nargin == 7
    splot = sTE-s;
else
    splot = s;
end


deltaHat = 64*x.^3.*(1-x).^3;
fontsize = 10;
figure(1)
hold on
tiledlayout(2,2,'TileSpacing','Compact')
nexttile
plot(t,deltaHat)
ylabel('$\hat{\delta}$','Interpreter','latex')
xlabel('$\hat{x}$','Interpreter','latex')
set(gca,'FontSize',fontsize)

nexttile
hold on
plot(x,sHat)
xlabel('$\hat{x}$','Interpreter','latex')
ylabel('$\hat{s}$','Interpreter','latex')
pbaspect([1 1 1])
axis equal
set(gca,'FontSize',fontsize)

nexttile([1 2])
hold on
plot(splot,deltaHat)
xlabel('s/c')
ylabel('$\hat{\delta}$','Interpreter','latex')
pbaspect([4 1 1])
axis([min(splot) max(splot) 0 1])
set(gca,'FontSize',fontsize)


end