function [xi,yi,xi2,yi2]=clip_curve(x,y,i1,i2,n)
% clip out shortest curve between i1 and i2, assuming x,y form a continuous
% loop


N = length(x);
s(1) = 0;
for i=2:N
dx = x(i)-x(i-1);
dy = y(i)-y(i-1);
s(i) = s(i-1) + sqrt(dx*dx + dy*dy);
end
s = s/s(N);


stot = 1;

s1 = s(i1);
s2 = s(i2);

del1 = abs(s1 - s2);
del2 = stot - del1;

if(del1<del2)

flip = 1;
if(i2<i1); flip=-1;end
    
xc = x(i1:flip:i2);
yc = y(i1:flip:i2);
sc = s(i1:flip:i2);

else
  
if(i1<i2)
xc = x([i1:-1:2 N:-1:i2]);
yc = y([i1:-1:2 N:-1:i2]);
sc = s([i1:-1:2 N:-1:i2]);
%sc(sc>0.5)=sc(sc>0.5)-1;
sc(end-(N-i2):end) = sc(end-(N-i2):end)-1;
else
xc = x([i2:-1:2 N:-1:i1]);
yc = y([i2:-1:2 N:-1:i1]);    
sc = s([i2:-1:2 N:-1:i1]); 
%sc(sc>0.5)=sc(sc>0.5)-1;
sc(end-(N-i1):end) = sc(end-(N-i1):end)-1;
end


end

%fc = linspace(0,1,length(sc));
%fi = linspace(0,1,n);
%si = interp1(fc,sc,fi,'linear');
si = linspace(sc(1),sc(end),n);
si2 = 1 - si;

xi = interp1([s-1 s(2:N-1) s+1],[x x(2:N-1) x],si,'spline');
yi = interp1([s-1 s(2:N-1) s+1],[y y(2:N-1) y],si,'spline');

xi2 = interp1([s-1 s(2:N-1) s+1],[x x(2:N-1) x],si2,'spline');
yi2 = interp1([s-1 s(2:N-1) s+1],[y y(2:N-1) y],si2,'spline');


return


