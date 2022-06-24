function fx = spacing(n,fx_expan,flag)

if(flag==1) 
    
fx(1) = 0;
df = 1;
for i=2:floor(n/2)
    fx(i) = fx(i-1) + df;
    df = df*fx_expan;
end
for i=floor(n/2):n
    fx(i) = fx(i-1) + df;
    df = df/fx_expan;
end

else
    
fx(1) = 0;
df = 1;
for i=2:n
    fx(i) = fx(i-1) + df;
    df = df*fx_expan;
end

end


fx = fx/fx(n);
fx(1) = 0;
fx(end) = 1;
    
return