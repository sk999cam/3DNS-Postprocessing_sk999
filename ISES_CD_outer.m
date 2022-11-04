function [Pr, Ct, CtEQ] = ISES_CD_outer(mF, xStart, K, pgterm, diffterm)

    if nargin < 3 || isempty(K)
        K = 4.2;
    end

    if nargin < 4 || isempty(pgterm)
        pgterm = false;
    end

    if nargin < 5 || isempty(diffterm)
        diffterm = false;
    end

    istart = mF.x2ind(xStart);
    s = mF.ssurf;
    Us = 0.5*mF.H_ke.*(1-4*(mF.H_k-1)./(3*mF.H));
    del = mF.theta.*(3.15+1.72./(mF.H_k - 1)) + mF.delStar;
    CtEQ = mF.H_ke.*(0.015./(1-Us)).*(mF.H_k-1).^3./(mF.H_k.^2.*mF.H);
    Ct(1:istart-1) = NaN;
    Ct(istart) = CtEQ(istart);

    for i=istart:length(s)-1
        dCt_ds = K*Ct(i)*(sqrt(CtEQ(i)) - sqrt(Ct(i)))/del(i);
        Ct(i+1) = Ct(i) + dCt_ds*(s(i+1)-s(i));
    end

    Pr = Ct.*(1-Us);


end