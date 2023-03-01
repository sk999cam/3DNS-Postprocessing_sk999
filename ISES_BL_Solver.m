function [Pr, Ct, CtEQ, pgterm, diffterm] = ISES_CD_outer(s, M, theta, H_ke, H, H_k, delStar, Ue, cf, istart, Kcorr, Kp, Kd)

    if nargin < 10 || isempty(Kcorr)
        Kcorr= 4.2;
    end

    Kcorr

    if nargin < 11 || isempty(Kp)
        Kp = 0;
    end

    if nargin < 12 || isempty(Kd)
        Kd = 0;
    end

    A = 6.7; B = 0.75; % MISES Eqm locus consts
    
    
    dUedx = (Ue(3:end) - Ue(1:end-2))./(s(3:end) - s(1:end-2));
    dUedx = [dUedx(1) dUedx dUedx(end)];
    


    
    Us = 0.5*H_ke.*(1-4*(H_k-1)./(3*H));
    del = theta.*(3.15+1.72./(H_k - 1)) + delStar;
    CtEQ = H_ke.*(0.015./(1-Us)).*(H_k-1).^3./(H_k.^2.*H);
   % CtEQ = 0.02456*((H_k-1)./H_k).^3./(1-Us);
    Ct(1:istart-1) = NaN;
    Ct(istart) = CtEQ(istart);


    for i=istart:length(s)-1
        dCt_ds = Kcorr*(sqrt(CtEQ(i)) - sqrt(Ct(i)));

        diffterm(i) = (2*del(i)*H(i)/(B*theta(i)))*(cf(i)/2 - ((H_k(i) - 1)/(A*H_k(i)))^2);
        dCt_ds = dCt_ds + Kd*diffterm(i);

        pgterm(i) = (2*del(i)/Ue(i))*dUedx(i);
        dCt_ds = dCt_ds - Kp*pgterm(i);
        
        dCt_ds = dCt_ds*Ct(i)/del(i);
        Ct(i+1) = Ct(i) + dCt_ds*(s(i+1)-s(i));
    end

    diffterm(end+1) = NaN;
    pgterm(end+1) = NaN;

    Pr = Ct.*(1-Us);


end

% function val = Hk(H, Me)
%     val = (H - 0.290*Me^2)/(1+0.113*Me^2);
% end

function val = fH(Me, Hk)
    val = 0.290*Me^2 + Hk*(1+0.113*Me^2);
end

function val = fHss(Me, Hk)
    val = (0.064/(Hk-0.8) + 0.251)*Me^2;
end

function val = fUs(Hs, Hk, H)
    val = 0.5*Hs*(1-4*(Hk-1)/(3*H));
end

function val = fCtEQ(Hs, Us, H, Hk)
    val = Hs*(0.015/(1-Us))*(Hk-1)^3/(Hk^2*H);
end

function val = fRet(Me, theta, T0, P0)
    cp = 1005;
    gam = 1.4;
    R = cp*(gam-1)/gam;
    ro0 = P0/(R*T0);
    T = T0/(1+0.5*(gam-1)*Me^2);
    ro = ro0/(1+0.5*(gam-1)*Me^2)^(1/(gam-1));
    Ue = Me*sqrt(gam*R*T);
    mu_ref = 5.83247e-004;
    Tref = 273;
    mu_s = 10.4;
    mue = mu_ref*(T/Tref)^(1.5)*(Tref+mu_s)/(T + mu_s);
    val = ro*Ue*theta/mue;
end

function fCf(Hk, Ret, Me)
    Fc = sqrt(1+0.2*Me^2);
    val = (0.3*exp(-1.33*Hk)*(log10(Ret/Fc))^(-1.74-0.31*Hk) ...
        + 0.00011*(tanh(4-Hk/0.875) - 1))/Fc;
end

function val = fHk(Hk, Me, Ret)  
    if Ret < 400
        H0 = 4;
    else
        H0 = 3+ 400/Ret;
    end
    Hk = 1.505+4/Ret+(0.165-1.6/sqrt(Ret))*
end