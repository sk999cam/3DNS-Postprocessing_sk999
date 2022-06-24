function r = fexpan(L_over_a,N)
% solve (L/a) - r(L/a) + r^N - 1 = 0

if(abs(L_over_a-N) < 1e-12)
r = 1.0
return
elseif( L_over_a>N )
r = 1.00001;
else
r = 0.99999;
end

for iter = 1:50
    
if(r<1.0)
r = 1 - (1 - r^N)/L_over_a;
elseif(r>1)
r = (1 - L_over_a*(1-r))^(1/N);
end

end

if (r<0.1) || (r>1.9)
    r=1.0
end

return