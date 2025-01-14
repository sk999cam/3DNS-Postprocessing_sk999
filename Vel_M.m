function V = Vel_M(M, T0, cp, gam)

if nargin < 2 || isempty(T0)
    T0 = 300;
end
if nargin < 3 || isempty(cp)
    cp = 1005;
end
if nargin < 4 || isempty(gam)
    gam = 1.4;
end

T = T0*T_T0(M, gam);
V = sqrt(2*cp*(T0-T));