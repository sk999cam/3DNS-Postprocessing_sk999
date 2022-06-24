function s = curve_length(x,y)

N = length(x);
s(1) = 0;

for i=2:N
dx = x(i)-x(i-1);
dy = y(i)-y(i-1);
ds = sqrt(dx*dx + dy*dy);

s(i) = s(i-1) + ds;
    
end



return