fclose all;

clear


nk = 100;
cp = 1005.0;
gam = 1.4;
cv = cp/gam;
rgas = cp-cv;

pitch =  8.1173923294907724e-01;%0.040578905958682;

 dir = '';% 'rans/';
 
%load([dir,'span.txt']);
load([dir,'blockdims.txt']);

%nk = length(span);
%nk = 1;
%span = linspace(0,1,nk);

ivisc =1;

minx = 1e12;
maxx = -1e12;
miny = 1e12;
maxy = -1e12;

for nb=1:size(blockdims,1)
    
disp(['Reading block ', num2str(nb)]);
    
span = load([dir,'span_',num2str(nb),'.txt']);
    
    fid = fopen([dir,'grid_',num2str(nb),'.txt'],'r');

   ni = blockdims(nb,1);
   nj = blockdims(nb,2);
   %nk = blockdims(nb,3);
   nk = 1;%length(span);
   
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

zz(1,1,1:nk) = span(1:nk); 
z = repmat(zz,[ni nj 1]);
x = repmat(x,[1 1 nk]);
y = repmat(y,[1 1 nk]);
clear zz

ro = zeros(ni,nj,nk);
ru = zeros(ni,nj,nk);
rv = zeros(ni,nj,nk);
rw = zeros(ni,nj,nk);
Et = zeros(ni,nj,nk);

txx = zeros(ni,nj,nk);
tyy = zeros(ni,nj,nk);
tzz = zeros(ni,nj,nk);
txy = zeros(ni,nj,nk);
txz = zeros(ni,nj,nk);
tyz = zeros(ni,nj,nk);
div = zeros(ni,nj,nk);
vort_x = zeros(ni,nj,nk);
vort_y = zeros(ni,nj,nk);
vort_z = zeros(ni,nj,nk);
mut = zeros(ni,nj,nk);


fid2 = fopen([dir,'flo2_',num2str(nb)],'r');    
A = fread(fid2,inf,'float64');
A = reshape(A,5,length(A)/5);

fid3 = fopen([dir,'nod2_',num2str(nb)],'r');    
B = fread(fid3,inf,'int');
B = reshape(B,3,length(B)/3);
% % 
% fid4 = fopen([dir,'visc_',num2str(nb)],'r');    
% C = fread(fid4,inf,'float64');
% C = reshape(C,6,length(C)/6);

 icount(1:ni*nj*nk) = 0;

for n=1:min([size(B,2) size(A,2)])
i = B(1,n);
j = B(2,n);
k = B(3,n);

nid = i + (j-1)*ni + (k-1)*ni*nj;

if(i<=ni & j<=nj & k<=nk & icount(nid)==0)

ro(i,j,k) = A(1,n);
ru(i,j,k) = A(2,n);
rv(i,j,k) = A(3,n);
rw(i,j,k) = A(4,n);
Et(i,j,k) = A(5,n);
% 
% txx(i,j,k) = C(1,n);
% tyy(i,j,k) = C(2,n);
% tzz(i,j,k) = C(3,n);
% txy(i,j,k) = C(4,n);
% txz(i,j,k) = C(5,n);
% tyz(i,j,k) = C(6,n);
% 
% div(i,j,k)    = C(7,n);
% vort_x(i,j,k) = C(8,n);
% vort_y(i,j,k) = C(9,n);
% vort_z(i,j,k) = C(10,n);
% mut(i,j,k)    = C(11,n);

end


%icount(nid) = icount(nid)+1;
icount(nid) = 1;

end


fclose(fid2);
fclose(fid3);
%fclose(fid4);

%stop


if min(min(min(x))) < minx
    minx = min(min(min(x)));
end

if max(max(max(x))) > maxx
    maxx = max(max(max(x)));
end

if min(min(min(y))) < miny
    miny = min(min(min(y)));
end

if max(max(max(y))) > maxy
    maxy = max(max(max(y)));
end

nkmid=floor((nk+1)*0.5);

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
% 
for k=1:nk
[dudx,dudy] = gradHO(squeeze(x(:,:,nkmid)),squeeze(y(:,:,nkmid)),squeeze(u(:,:,k)));
[dvdx,dvdy] = gradHO(squeeze(x(:,:,nkmid)),squeeze(y(:,:,nkmid)),squeeze(v(:,:,k)));
[drdx,drdy] = gradHO(squeeze(x(:,:,nkmid)),squeeze(y(:,:,nkmid)),squeeze(ro(:,:,k)));
%div = drudx + drvdy;
f =sqrt( drdx.*drdx + drdy.*drdy );
strain(:,:,k) = dudy + dvdx;
vortz(:,:,k) = dudy - dvdx;
schli(:,:,k) = sqrt(drdx.*drdx + drdy.*drdy);
end

alpha = atan(v./u)*180/pi;

figure(1)

urms = repmat(std(u,0,3),[1 1 nk]);
vrms = repmat(std(v,0,3),[1 1 nk]);
wrms = repmat(std(w,0,3),[1 1 nk]);


prop = p;%mut;%vortz;%p;%mut;%vort_z;%mach;%vortz;%strain;%vortz;%100*sqrt((urms.*urms + vrms.*vrms + wrms.*wrms)/3.0)/13.8;%vortz - repmat(mean(vortz,3),[1 1 nk]);%vortz;%schli;%vortz;%w;%v - vref;%mach;%p-1e5;%w;%vortz;

minp(nb) = min(min(min(prop)));
maxp(nb) = max(max(max(prop)));

% 
pcolor(x(:,:,nkmid),y(:,:,nkmid),prop(:,:,nkmid))
hold on
pcolor(x(:,:,nkmid),y(:,:,nkmid)+pitch,prop(:,:,nkmid))
pcolor(x(:,:,nkmid),y(:,:,nkmid)-pitch,prop(:,:,nkmid))


if(nb == 3 | nb ==4 | nb ==5 | nb ==7)
figure(2)
plot(squeeze(x(:,end,:)),p(:,end,:),'b')
if(nb==5)
plot(squeeze(x(:,end,:)),p(:,end,:),'bx')
end
hold on
end

if(nb==1 || nb==2 || nb==8 || nb==9)
figure(2)
plot(x(:,end,:),p(:,end,:),'r'),hold on
plot(x(:,end-1,:),p(:,end-1,:),'b')
plot(x(:,floor(end/2),:),p(:,floor(end/2),:),'g')
end
% 
if nk>1
figure(6)
surf(squeeze(x(:,end,:)),squeeze(y(:,end,:)),squeeze(z(:,end,:)),squeeze(prop(:,end,:)))
hold on

end

% 
clear vortz schli strain

% 
%        tdata=[];
%        tdata.Nvar=10;
%        tdata.varnames={'x','y','z','p','T','ro','u','v','w','mut'};
%        tdata.cubes(1).zonename=['block ',num2str(nb)];
%        tdata.cubes(1).x=x;
%        tdata.cubes(1).y=y;  
%        tdata.cubes(1).z=z;
%        tdata.cubes(1).v(1,:,:,:)=p;
%        tdata.cubes(1).v(2,:,:,:)=T;
%        tdata.cubes(1).v(3,:,:,:)=ro;
%        tdata.cubes(1).v(4,:,:,:)=u;
%        tdata.cubes(1).v(5,:,:,:)=v;
%        tdata.cubes(1).v(6,:,:,:)=w;
%        tdata.cubes(1).v(7,:,:,:)=mut;
%        tdata.vformat(1:10) = 1; 
%        tdata.cubes.solutiontime=1;
%        
%        mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])



end

% % 
% % figure(4)
% % camlight; lighting phong
% % axis equal
% % axis off
% % 
figure(6)
shading interp
axis equal



figure(1)
axis([minx maxx miny-pitch*0.5 maxy+pitch*0.5])

%caxis([-2e6 2e6])
axis equal
shading interp