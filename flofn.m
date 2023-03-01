function f = flofn(M, gam)
if nargin < 2
    gam = 1.4;
end

f = gam * (gam-1)^(-0.5) * M * (1+((gam-1)/2)*M^2)^(-0.5*(gam+1)/(gam-1));

end