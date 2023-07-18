function M = M_VT0(V,T0,gam,cp)

    if nargin < 3 || isempty(gam)
        gam = 1.4;
    end
    if nargin < 4 || isempty(cp)
        cp = 1005;
    end

    rgas = cp*(gam-1)/gam;

    T = T0 - V.^2/(2*cp);
    M = V./sqrt(gam*rgas*T);
end