
dely = 1.0;
y1 = 0.01;
N = 52;

fac = 1.0
if(dely/y1>N-1); fac = 1.001; end
if(dely/y1<N-1); fac = 0.999; end
    
    for i=1:10
      if(fac>1.0)
      fac=(1.0 - (dely/y1)*(1.0 - fac))^(1.0/(N-1))
      elseif(fac<1.0) then
      fac =1.0 - (dely/y1)*(1.0 - (fac^float(n-1)))
      end
    end
