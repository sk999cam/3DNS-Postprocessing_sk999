
function [lMin, lMax, minL, maxL] = gridSpac(xb,yb,zb)

    [ni,~,~]=size(xb);
    Z=squeeze(zb(1:ni/2,:,:));
    Y=squeeze(yb(1:ni/2,:,:));
    X=squeeze(xb(1:ni/2,:,:));

    

    xyz=cat(4,X,Y,Z);
    xyz_corner=xyz([1 end],[1 end],[1 end],:);

    dI=xyz(2:end,:,:,:)-xyz(1:end-1,:,:,:);
    lI=vecnorm(dI,2,4);
    lMinI=min(lI(:));
    lMaxI=max(lI(:));

    dI=xyz_corner(2,:,:,:)-xyz_corner(1,:,:,:);
    lI=vecnorm(dI,2,4);
    edgeMinI=min(lI(:));
    edgeMaxI=max(lI(:));

    
    dJ=xyz(:,2:end,:,:)-xyz(:,1:end-1,:,:);
    lJ=vecnorm(dJ,2,4);
    lMinJ=min(lJ(:));
    lMaxJ=max(lJ(:));

    dJ=xyz_corner(:,2,:,:)-xyz_corner(:,1,:,:);
    lJ=vecnorm(dJ,2,4);
    edgeMinJ=min(lJ(:));
    edgeMaxJ=max(lJ(:));

    dK=xyz(:,:,2:end,:)-xyz(:,:,1:end-1,:);
    lK=vecnorm(dK,2,4);
    lMinK=min(lK(:));
    lMaxK=max(lK(:));

    dK=xyz_corner(:,:,2,:)-xyz_corner(:,:,1,:);
    lK=vecnorm(dK,2,4);
    edgeMinK=min(lK(:));
    edgeMaxK=max(lK(:));

    lMin=min([lMinI,lMinJ,lMinK]);
    lMax=max([lMaxI,lMaxJ,lMaxK]);

    minL = min([edgeMinI,edgeMinJ,edgeMinK]);
    maxL=max([edgeMaxI,edgeMaxJ,edgeMaxK]);
    
return
    
