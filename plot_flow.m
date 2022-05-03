close all
clear

nfig = 1;

cp = 1005.0;
gam = 1.4;
cv = cp/gam;
rgas = cp-cv;

spanExtent = 0.1;

temp_slash = '/'; if ispc, temp_slash = '\'; end
dir = 'r150_cwl90_hf2';
 
%load([dir,'span.txt']);
load([dir,temp_slash,'blockdims.txt']);

%%

%nk = length(span);
%nk = 1;
nk = blockdims(1,3);
span = linspace(0,spanExtent,nk);


ivisc =0;

minx = 1e12;
maxx = -1e12;
miny = 1e12;
maxy = -1e12;

for nb=1:size(blockdims,1)

    
    fid = fopen([dir,temp_slash,'grid_',num2str(nb),'.txt'],'r');

   ni = blockdims(nb,1);
   nj = blockdims(nb,2);
      
    x = zeros(ni,nj);
    y = zeros(ni,nj);
    
   
   for j=1:nj
   for i=1:ni
   [A]=fscanf(fid, '%f %f',2);
   x(i,j) = A(1);
   y(i,j) = A(2);
   end
   end
   
   fclose(fid);


    zz(1,1,1:nk) = span; 
    z = repmat(zz,[ni nj 1]);
    x = repmat(x,[1 1 nk]);
    y = repmat(y,[1 1 nk]);
    
    ro = zeros(ni,nj,nk);
    ru = zeros(ni,nj,nk);
    rv = zeros(ni,nj,nk);
    rw = zeros(ni,nj,nk);
    Et = zeros(ni,nj,nk);
    
    fid2 = fopen([dir,temp_slash,'flo2_',num2str(nb)],'r');    
    A = fread(fid2,inf,'float64');
    A = reshape(A,5,length(A)/5);
    
    fid3 = fopen([dir,temp_slash,'nod2_',num2str(nb)],'r');    
    B = fread(fid3,inf,'uint32');
    B = reshape(B,3,length(B)/3);
    
    for n=1:size(A,2)
    i = B(1,n);
    j = B(2,n);
    k = B(3,n);
    ro(i,j,k) = A(1,n);
    ru(i,j,k) = A(2,n);
    rv(i,j,k) = A(3,n);
    rw(i,j,k) = A(4,n);
    Et(i,j,k) = A(5,n);
    end
    
    
    fclose(fid2);
    fclose(fid3);
    
    
    
    
    if(nb==2)
        ylow = min(y(1,:));
        xmin = min(min(min(x)));
    end
    if(nb==3)
        yhigh = max(y(1,:));

    end
    if(nb==10)
        xmax = max(max(max(x)));
    end
    
    %pitch = 0.059045733332343;
    
    nkmid = floor((nk+1)/2)
    %nkmid = 5;
    
    u = ru./ro;
    v = rv./ro;
    w = rw./ro;
    p = (gam-1)*(Et - 0.5*(u.*u + v.*v + w.*w).*ro);
    T = p./(ro*rgas);
    
    s = cp*log(T/300) - rgas*log(p/1e5);
    vel = sqrt(u.*u + v.*v + w.*w);
    mach = vel./sqrt(gam*rgas*T);
    
    To = T.*(1.0 + (gam-1)*0.5*mach.*mach);
    po = p.*((To./T).^(gam/(gam-1)));
    
    prop = mach;
    %prop = w;
    %prop = ro;
    
    figure(nfig)
    pcolor(x(:,:,nkmid),y(:,:,nkmid),prop(:,:,nkmid))
    hold on

end

figure(nfig)
axis([xmin xmax ylow yhigh])
colorbar()
%caxis([-0.5 0.5])
axis equal
shading interp
savefig([dir,temp_slash,'M.png'])