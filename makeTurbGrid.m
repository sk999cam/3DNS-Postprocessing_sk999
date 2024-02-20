clear;

gridfldr =  '/mnt/Tank/DNSData/Case2a/old/3Re/fineGrid3/';
namePath = 'grid_turb_ax_f';
gridRead = 'ax';

span = 0.009465;

nmBlc = [2,1];%ax  %[1,2];%rad

%% USE THIS IF YOU WANT A SPECIFIC NUMBER OF NODES
iStart = 0; 
iEnd = 715; 
I = iEnd-iStart;

%% USE THIS IF YOU WANT TO HAVE A SPECIFIC LENGTH
L = 0.0408*2;% ax   %0.009533835077005; %rad

NJ = 0;


DIM = importdata([gridfldr,'blockdims.txt']);
grid = cell(2,1);

%% READ IN GRID AND BLOCKS
for nb = nmBlc   
    
    disp(nb)
    
    ni = DIM(nb,1);
    nj = DIM(nb,2);
    nk = DIM(nb,3);    
    
    if strcmp(gridRead,'txt')
        G = importdata([gridfldr,'grid_',num2str(nb),'.txt']);
        x = reshape(G(:,1),[ni,nj,nk]);
        y = reshape(G(:,2),[ni,nj,nk]);
        z = reshape(G(:,3),[ni,nj,nk]);
    elseif strcmp(gridRead,'bin')
        gridName = strcat(gridfldr,'grid_',num2str(nb));
        fid_g = fopen(gridName,'r'); %
        G = fread(fid_g,inf,'float64');       
        
        G = reshape(G,3,length(G)/3);  
        x = reshape(G(1,:),[ni,nj,nk]);
        y = reshape(G(2,:),[ni,nj,nk]);
        z = reshape(G(3,:),[ni,nj,nk]);
        fclose(fid_g); 
    elseif strcmp(gridRead,'ax')
        z = 0;
        G = importdata([gridfldr,'grid_',num2str(nb),'.txt']);
        x = reshape(G(:,1),[ni,nj]);
        y = reshape(G(:,2),[ni,nj]);

        z(1,1,1:nk) = linspace(0,span,nk);
        x = repmat(x,[1,1,nk]);
        y = repmat(y,[1,1,nk]);
        z = repmat(z,[ni,nj,1]);
    else        
        disp('Can\"t read grid')
        break
    end
    clear G
    
    grid{nb}.x = x;
    grid{nb}.y = y;
    grid{nb}.z = z;

    NJ = NJ+DIM(nb,2);

end

NJ = NJ-1;

xb = cat(2,grid{1}.x,grid{2}.x(:,2:end,:));
yb = cat(2,grid{1}.y,grid{2}.y(:,2:end,:));
zb = cat(2,grid{1}.z,grid{2}.z(:,2:end,:));


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


dx = mean(mean(xb(2,:,:)-xb(1,:,:)))

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
fid_f = fopen(namePath,'w'); 
fwrite(fid_f,G(:),'float64');
fclose(fid_f);

maxL = dx*double(I)
