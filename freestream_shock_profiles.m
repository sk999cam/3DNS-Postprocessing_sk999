function [x, ro, u, Et, p, M] = freestream_shock_profiles(Min, blk, gas, bcs, xShock, Lshock)

    if nargin < 5
        xShock = 0.5;
        Lshock = 0.005;
    end

    x = blk.x{1}(:,1,1);
    gam = gas.gam;
    cp = gas.cp;
    rgas = cp*(gam-1)/gam;
    
    % Pre-shock conditions
    fM = 1+0.5*(gam-1)*Min^2;
    pin = bcs.Poin*fM^(-gam/(gam-1));
    tin = bcs.Toin/fM;
    roin = pin/(rgas*tin);
    vin = Min*sqrt(gam*rgas*tin);
    Etin = pin/(gam-1) + 0.5*roin*vin^2;
    
    % Post shock conditions
    Ms = sqrt(fM/(gam*Min^2 - 0.5*(gam-1)));
    ps = pin*(1+2*gam*(Min^2-1)/(gam+1));
    ros = 0.5*roin*(gam+1)*Min^2/fM;
    Ts = ps/(ros*rgas);
    vs = Ms*sqrt(gam*rgas*Ts);
    Ets = ps/(gam-1) + 0.5*ros*vs^2;
    
    
    shfn = tanh((x-xShock)/Lshock);
    
    u = 0.5*((vin + vs) - (vin-vs)*shfn);
    ro = 0.5*((roin + ros) - (roin-ros)*shfn);
    Et = 0.5*((Etin+Ets) - (Etin-Ets)*shfn);
    p = (gam-1)*(Et - 0.5*ro.*u.^2);
    T = (Et./ro - u.^2/2)*gam/cp;
    M = u./sqrt(gam*rgas*T);

end