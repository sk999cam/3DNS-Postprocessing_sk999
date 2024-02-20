function preTestSpectrumGrid(paramsName,gridName)
    %function to plot prescribed spectrum, given an inlet box grid as well
    %as a parameter file.
    

    %% READ IN GRID AND BLOCKS

    eval(paramsName)

    fid_g = fopen(gridName,'r'); %
    G = fread(fid_g,inf,'float64');       

    G = reshape(G,3,length(G)/3);  
    xb = reshape(G(1,:),[ni,nj,nk]);
    yb = reshape(G(2,:),[ni,nj,nk]);
    zb = reshape(G(3,:),[ni,nj,nk]);
    fclose(fid_g); 


    %{
    plot3(xb(:,1,1),yb(:,1,1),zb(:,1,1),'k')
    hold on
    plot3(xb(:,end,1),yb(:,end,1),zb(:,end,1),'k')
    plot3(xb(1,:,1),yb(1,:,1),zb(1,:,1),'k')
    plot3(xb(:,1,end),yb(:,1,end),zb(:,1,end),'k')
    hold on
    plot3(xb(:,end,end),yb(:,end,end),zb(:,end,end),'k')
    plot3(xb(1,:,end),yb(1,:,end),zb(1,:,end),'k')
    plot3(squeeze(xb(1,1,:)),squeeze(yb(1,1,:)),squeeze(zb(1,1,:)),'k')
    hold on
    plot3(squeeze(xb(end,end,:)),squeeze(yb(end,end,:)),squeeze(zb(end,end,:)),'k')
    plot3(squeeze(xb(1,end,:)),squeeze(yb(1,end,:)),squeeze(zb(1,end,:)),'k')
    plot3(squeeze(xb(end,1,:)),squeeze(yb(end,1,:)),squeeze(zb(end,1,:)),'k')
    plot3(squeeze(xb(end,:,1)),squeeze(yb(end,:,1)),squeeze(zb(end,:,1)),'k')
    plot3(squeeze(xb(end,:,end)),squeeze(yb(end,:,end)),squeeze(zb(end,:,end)),'k')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %}
   
    lMin = gridSpac(xb,yb,zb)
    
    wrap = @(ke,keta,dltk,b,tol) estInt(ke,keta,dltk,b,tol); 

    tol = 0.00001;

    if L ~= -1 && eps == -1
        eps =  (Tu*vIn)^3/L*sqrt(3/2)^3   
    elseif L ==-1 && eps ~= -1
        L = (Tu*vIn)^3/eps*sqrt(3/2)^3;  
    elseif L ~= -1 && eps ~= -1
        disp('Both eddy lengthscale and dissipation defined. Defaulting to eddy lengthscale.')
        eps = (Tu*vIn)^3/L*sqrt(3/2)^3;  
    else
        disp('RTFD: Wrong usage of L and dissipation')
        return
    end

    %rms of velocity fluctuations
    u = Tu*vIn;
    %ch'ic length scale-based Re
    %ReL = u*Lc/nui;
    ReL = u*L/(nu)
    %Wave length related to location of energy maximum (at sqrt(12/5)ke)
    %relationship to int. length scale fixed for isotropic turbulence
    keSt = 0.746834/(L*0.09);
    %Kolmogorov wave length (= eps^1/4/nu^3/4 = u^(3/4)/(nu^(3/4)*Lint^(1/4))
    %keta = ReL.^(3/4)/(Lint^(1/4)*Lc^(3/4));
    keta = eps^(1/4)/nu^(3/4);

    %Grid spacing
    %dg = maxL/(N-1)
    %minimum resolved wave number (i.e. longest waves - want to take minL
    %here rather than maxL, so that we won't end up with waves too large
    %for one direction.
    kmi = 2*pi/minL;
    %maximum resolved wave number
    kmx = pi/lMin;  


    mn = linspace(1,M,M);  
    km = zeros(M,1);

    ce = 2*nu;
    cu = 2/3;
    b = ce*u^2/(eps*cu);

    %set wavenumber as centre of bin
    for m = 1:M    
        km(m) = kmi+(kmx-kmi)/M*(mn(m)-0.5);
    end

    %width of wave-number bins
    dltk = (max(km)-min(km))/(M-1);

    try
        %Find ke..
        fun = @(x) wrap(x,keta,dltk,b,tol);    % function of x alone
    catch
        
    end
    ke = fzero(fun,[1, keSt]);
    
    disp(['Calculated k_e as ', num2str(ke)])
    [dsint, kint, alpha] = calcAlpha(ke,dltk,keta,nu,u,tol);

    
    disp(['Calculated alpha as ', num2str(alpha)])

    Ek = zeros(M,1);

    for m = 1:M
        bfk = km(m)/ke;
        expc = exp(-2*(km(m)/keta)^2);
        Ek(m) = alpha*u^2/ke*bfk^4/(1+bfk^2)^(17/6)*expc;
    end

    E = Ek/u^2*keta;
    k = km/keta;
    
    disp(1/u^2*keta)
    disp(keta)
    
    figure()
    
    loglog(k,E,'LineWidth',1.5);
    [~,minDex] = min(abs(k-1.5));
    xLimits = [min(k) 1.5];                   %# Limits for the x axis
    yLimits = [E(minDex) max(E)*1.1 ];                      %# Limits for the y axis
    %logScale = diff(yLimits)/diff(xLimits);  %# Scale between the x and y ranges
    %powerScale = diff(log10(yLimits))/...    %# Scale between the x and y powers
    %             diff(log10(xLimits));
    %set(gca,'Xlim',xLimits,'YLim',yLimits,...              %# Set the limits and the
    %        'DataAspectRatio',[1 logScale/powerScale 1]);  %#   data aspect ratio
    set(gca,'Xlim',xLimits,'YLim',yLimits);    
    grid on
    hold on
    i1 = find(E>=0.8*max(E),1,'last');
    i2 = find(E<=0.05*max(E),1,'first');
    Bp = polyfit(log10(k(i1:i2)), log10(E(i1:i2)), 1);
    Yp = polyval(Bp,log10(k(i1:i2)));
    Yp=10.^Yp;
    loglog(k(i1:i2),Yp,'-r');
    
    %find approximate gradient in 
    m = (log10(Yp(end))-log10(Yp(1)))/(log10(k(i2))-log10(k(i1)))
    
    keta
    xlabel('$\kappa/\kappa_\eta$','interpreter','latex','FontSize',15)
    ylabel('$E(\kappa)\kappa_\eta/u''^2$','interpreter','latex','FontSize',15)
   
end
