function M = flofn_inv(f, gam)

    if nargin < 2
        gam = 1.4;
    end
    
    M0 = 0;
    M1 = 1.0;
    M = 0.5*(M0+M1);
    fnow = flofn(M,gam);
    
    while abs(f/fnow - 1) > 0.0001
        if fnow > f
            M1 = M;
        else
            M0 = M;
        end
        M = 0.5*(M0+M1);
        fnow = flofn(M,gam);
    end
end