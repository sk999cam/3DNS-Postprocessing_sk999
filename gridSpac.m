function [lMin, lMax, minL, maxL] = gridSpac(xb,yb,zb)

    [ni,~,~]=size(xb);
    Z=squeeze(zb(1:ni/2,:,:));
    Y=squeeze(yb(1:ni/2,:,:));
    X=squeeze(xb(1:ni/2,:,:));

    

    xyz=cat(4,X,Y,Z);

    dI=xyz(2:end,:,:,:)-xyz(1:end-1,:,:,:);
    lI=vecnorm(dI,2,4);
    lMinI=min(lI(:));
    lMaxI=max(lI(:));
    
    dJ=xyz(:,2:end,:,:)-xyz(:,1:end-1,:,:);
    lJ=vecnorm(dJ,2,4);
    lMinJ=min(lJ(:));
    lMaxJ=max(lJ(:));

    dK=xyz(:,:,2:end,:)-xyz(:,:,1:end-1,:);
    lK=vecnorm(dK,2,4);
    lMinK=min(lK(:));
    lMaxK=max(lK(:));

    lMin=min([lMinI,lMinJ,lMinK]);
    lMax=max([lMaxI,lMaxJ,lMaxK]);


    dI=xyz(end,:,:,:)-xyz(1,:,:,:);
    lI=vecnorm(dI,2,4);
    lMinI=min(lI(:));
    lMaxI=max(lI(:));
    
    dJ=xyz(:,end,:,:)-xyz(:,1,:,:);
    lJ=vecnorm(dJ,2,4);
    lMinJ=min(lJ(:));
    lMaxJ=max(lJ(:));

    dK=xyz(:,:,end,:)-xyz(:,:,1,:);
    lK=vecnorm(dK,2,4);
    lMinK=min(lK(:));
    lMaxK=max(lK(:));

    minL=min([lMinI,lMinJ,lMinK]);
    maxL=max([lMaxI,lMaxJ,lMaxK]);
    
return
    
