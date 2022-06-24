function dy = Re2offset(Re, yplus, c)
if nargin < 3
    c = 1;
end

dy = yplus*sqrt(80)*Re^(-13/14)*c;
end