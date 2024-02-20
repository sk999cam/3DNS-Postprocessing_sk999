function [dsint,kint,alpha] = calcAlpha(ke,dltk,keta,nui,u,tol)

    
    %set wavenumber as centre of bin
    %for m = 1:nmM    
    %    km(m) = kmi+(kmx-kmi)/nmM*(mn(m)-0.5);
    %end

    %dltk = (max(km)-min(km))/(nmM-1);
    k = dltk/2;
    kint = 0;
    dsint = 0;

    while k<keta

        bfk = k/ke;
        expc = exp(-2*(k/keta)^2);
        Ek = u^2/ke*bfk^4/(1+bfk^2)^(17/6)*expc;
        kint = kint+Ek*dltk;
        dsint = dsint+Ek*dltk*k^2;
        k = k+dltk;

    end

    %{
    while 1

        bfk = k/ke;
        expc = exp(-2*(k/keta)^2);
        Ek = u^2/ke*bfk^4/(1+bfk^2)^(17/6)*expc;
        int = int+Ek*dltk;
        dissint = dissint+Ek*dltk*k^2;

        resE = 2/3*k*Ek;

        if resE/int<tol
            break
        end
        k = k+dltk;

    end
    %}

    while 1

            bfk = k/ke;
            expc = exp(-2*(k/keta)^2);
            Ek = 1/ke*bfk^4/(1+bfk^2)^(17/6)*expc;
            kint = kint+Ek*dltk;
            dsint = dsint+Ek*dltk*k^2;

            %resk = keta*ke^(2/3)/4*exp(-2*(k/keta)^2)/k^(8/3);   
            resk = 2/3*k*Ek;
            if resk/kint<tol
                break
            end
            k = k+dltk;

    end

    tke = 3/2*u^2;

    alpha = tke/kint;
    dsint = 2*nui*alpha*dsint;
    kint = alpha*kint;
    

end
