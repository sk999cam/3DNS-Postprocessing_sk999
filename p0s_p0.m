function f = p0s_p0(M, gam)

if nargin < 2
    gam = 1.4;
end

f = ( (gam+1)*M^2/2 / (1 + (gam-1)*M^2/2) )^(gam/(gam-1)) * ( 2*gam*M^2/(gam+1) - (gam-1)/(gam+1) )^(1/(1-gam));

end