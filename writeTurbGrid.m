function [xb,yb,zb,dx] = writeTurbGrid(inblocks, span, casedir, L)

    if nargin < 4
        L = 1;
    end

    NJ = 0;
    DIM = importdata(fullfile(casedir,'blockdims.txt'));
    grid = cell(2,1);

    for nb = inblocks
        disp(nb);
        ni = DIM(nb,1)
        nj = DIM(nb,2)
        nk = DIM(nb,3)
        
        z = 0;
        G = importdata(fullfile(casedir,['grid_',num2str(nb),'.txt']));
        x = reshape(G(:,1),[ni,nj]);
        y = reshape(G(:,2),[ni,nj]);

        z(1,1,1:nk) = linspace(0,span,nk);
        x = repmat(x,[1,1,nk]);
        y = repmat(y,[1,1,nk]);
        z = repmat(z,[ni,nj,1]);

        grid{nb}.x = x;
        grid{nb}.y = y;
        grid{nb}.z = z;
    
        NJ = NJ+DIM(nb,2);

    end

    NJ = NJ-2;

    xb = cat(2,grid{2}.x,grid{1}.x(:,2:end,:),grid{3}.x(:,2:end,:));
    yb = cat(2,grid{2}.y,grid{1}.y(:,2:end,:),grid{3}.y(:,2:end,:));
    zb = cat(2,grid{2}.z,grid{1}.z(:,2:end,:),grid{3}.z(:,2:end,:));
    
    
    figure()
    plot3(xb(:,1,1),yb(:,1,1),zb(:,1,1))
    hold on
    plot3(xb(:,end,1),yb(:,end,1),zb(:,end,1))
    plot3(xb(1,:,1),yb(1,:,1),zb(1,:,1))
    plot3(xb(:,1,end),yb(:,1,end),zb(:,1,end))
    hold on
    plot3(xb(:,end,end),yb(:,end,end),zb(:,end,end))
    plot3(xb(1,:,end),yb(1,:,end),zb(1,:,end))
    plot3(squeeze(xb(1,1,:)),squeeze(yb(1,1,:)),squeeze(zb(1,1,:)))
    hold on
    plot3(squeeze(xb(end,end,:)),squeeze(yb(end,end,:)),squeeze(zb(end,end,:)))
    plot3(squeeze(xb(1,end,:)),squeeze(yb(1,end,:)),squeeze(zb(1,end,:)))
    plot3(squeeze(xb(end,1,:)),squeeze(yb(end,1,:)),squeeze(zb(end,1,:)))
    plot3(squeeze(xb(end,:,1)),squeeze(yb(end,:,1)),squeeze(zb(end,:,1)))
    plot3(squeeze(xb(end,:,end)),squeeze(yb(end,:,end)),squeeze(zb(end,:,end)))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    
    dx = mean(mean(xb(2,:,:)-xb(1,:,:)));
    
    if L > 0 
        I = int32(round(L/dx));
        iStart = 0;
        iEnd = I;
    end
        
    xb = repmat(xb(1,:,:),I,1,1);
    yb = repmat(yb(1,:,:),I,1,1);
    zb = repmat(zb(1,:,:),I,1,1);
    
    xfac = linspace(double(iStart),double(iEnd-1),double(I))';
    
    xfac = xfac*dx;
    
    xfac = repmat(xfac,1,NJ,nk);
    
    xb = xb+xfac;
    
    clear grid x y z
    
    figure()
    
    plot3(xb(:,1,1),yb(:,1,1),zb(:,1,1))
    hold on
    plot3(xb(:,end,1),yb(:,end,1),zb(:,end,1))
    plot3(xb(1,:,1),yb(1,:,1),zb(1,:,1))
    plot3(xb(:,1,end),yb(:,1,end),zb(:,1,end))
    hold on
    plot3(xb(:,end,end),yb(:,end,end),zb(:,end,end))
    plot3(xb(1,:,end),yb(1,:,end),zb(1,:,end))
    plot3(squeeze(xb(1,1,:)),squeeze(yb(1,1,:)),squeeze(zb(1,1,:)))
    hold on
    plot3(squeeze(xb(end,end,:)),squeeze(yb(end,end,:)),squeeze(zb(end,end,:)))
    plot3(squeeze(xb(1,end,:)),squeeze(yb(1,end,:)),squeeze(zb(1,end,:)))
    plot3(squeeze(xb(end,1,:)),squeeze(yb(end,1,:)),squeeze(zb(end,1,:)))
    plot3(squeeze(xb(end,:,1)),squeeze(yb(end,:,1)),squeeze(zb(end,:,1)))
    plot3(squeeze(xb(end,:,end)),squeeze(yb(end,:,end)),squeeze(zb(end,:,end)))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    G = cat(2,xb(:),yb(:),zb(:));
    G = G';
    turbdir = fullfile(casedir,'turb-gen');
    if ~exist(turbdir,'dir')
        mkdir(turbdir)
    end
    fid_f = fopen(fullfile(turbdir,'turb-grid.dat'),'w'); 
    fwrite(fid_f,G(:),'float64');
    fclose(fid_f);
    
    maxL = dx*double(I);

end


