function [u,v,w] = turbid()

Ly = 0.667;
Lz = 0.1;
Lx = Ly*4;

ni = 444;
nj = 111;
nk = 64;

dx = Lx/ni;
dy = Ly/nj;
dz = Lz/nk;

[del,i] = max([dx,dy,dz]);

% set grid size
turb.nx = ni;
turb.ny = nj;
turb.nz = nk;

% set min length-scale
len = 6*dy;
% points per min length-scale
turb.lsx = floor(len/dx);
turb.lsy = floor(len/dy);
turb.lsz = floor(len/dz);

% outer length-scale set to 2 x min
turb.nfx = 2*turb.lsx;
turb.nfy = 2*turb.lsy;
turb.nfz = 2*turb.lsz;

% seed for random numbers
turb.iseed = 10023459;

% initialize flow
turb.velvecs = [];
u = zeros(ni,nj,nk);
v = zeros(ni,nj,nk);
w = zeros(ni,nj,nk);
turb = turin(turb);   

for i=1:turb.nx   
disp(['Processing i-slice ',num2str(i)])
turb = turin(turb);    
u(i,:,:) = turb.velvecs.velx;
v(i,:,:) = turb.velvecs.vely;
w(i,:,:) = turb.velvecs.velz;
end

% make periodic in i,j,k
u=periodic_ijk(u);
v=periodic_ijk(v);
w=periodic_ijk(w);

% remove mean and normalize by std
u = u - repmat(mean(u),[ni,1,1]);
v = v - repmat(mean(v),[ni,1,1]);
w = w - repmat(mean(w),[ni,1,1]);
% 
up = repmat(std(u),[ni,1,1]);
vp = repmat(std(v),[ni,1,1]);
wp = repmat(std(w),[ni,1,1]);

u = u./up;
v = v./vp;
w = w./wp;

% compute and plot spectra
N = ni*nj*nk;

v1 = fftn(u);
v2 = fftn(v);
v3 = fftn(w);

Ek = (v1.*conj(v1) + v2.*conj(v2) + v3.*conj(v3))/N;
Ek_1 = (v1.*conj(v1))/N;
Ek_2 = (v2.*conj(v2))/N;
Ek_3 = (v3.*conj(v3))/N;

Ek_x = squeeze(sum(sum(Ek,2),3))/(nj*nk);
Ek_y = squeeze(sum(sum(Ek,1),3))/(ni*nk);
Ek_z = squeeze(sum(sum(Ek,1),2))/(ni*nj);

Ek_1_x = squeeze(sum(sum(Ek_1,2),3));
Ek_1_y = squeeze(sum(sum(Ek_1,1),3));
Ek_1_z = squeeze(sum(sum(Ek_1,1),2));

Ek_2_x = squeeze(sum(sum(Ek_2,2),3));
Ek_2_y = squeeze(sum(sum(Ek_2,1),3));
Ek_2_z = squeeze(sum(sum(Ek_2,1),2));

Ek_3_x = squeeze(sum(sum(Ek_3,2),3));
Ek_3_y = squeeze(sum(sum(Ek_3,1),3));
Ek_3_z = squeeze(sum(sum(Ek_3,1),2));

k_x = linspace(0,ni/Lx,ni);
k_y = linspace(0,nj/Ly,nj);
k_z = linspace(0,nk/Lz,nk);

Ekm = mean(mean(mean(Ek)));
k_kol =  linspace(0,max([ni,nj,nk])/min([Lx,Ly,Lz]),max([ni,nj,nk]));
Ek_kol = 100*Ekm*(k_kol.^(-5/3));


figure(1)

subplot(2,2,1)
loglog(k_x,Ek_x,k_y,Ek_y,k_z,Ek_z);
hold on
loglog(k_kol,Ek_kol,'r--');

subplot(2,2,2)
loglog(k_x,Ek_1_x,k_x,Ek_2_x,k_x,Ek_3_x);
hold on
loglog(k_kol,Ek_kol,'r--');

subplot(2,2,3)
loglog(k_y,Ek_1_y,k_y,Ek_2_y,k_y,Ek_3_y);
hold on
loglog(k_kol,Ek_kol,'r--');

subplot(2,2,4)
loglog(k_z,Ek_1_z,k_z,Ek_2_z,k_z,Ek_3_z);
hold on
loglog(k_kol,Ek_kol,'r--');

% write to file
fid = fopen('inflow_turb.dat','w');
A = zeros(3,ni,nj,nk);
A(1,:,:,:) = u;
A(2,:,:,:) = v;
A(3,:,:,:) = w;
A = reshape(A,3*ni*nj*nk,1,1,1); 
fwrite(fid,A,'float64');
fclose(fid);

return
end
	   
       
function turb = turin(turb)
% Method of Philips and Fyfe implemented in Matlab
% ! If startup then initialize full random number fields:
% ! we are starting a new problem; if not startup then
% ! shift random fields and calculate temporally and
% ! spatially correlated data.
% ! lsx and lsy are the length scales as multiples
% ! of gridcells.
% ! The filter widths nfx, nfy, and nfz should be at least
% ! twice the length scales. Normally nfy = nfz.

ny = turb.ny;
nz = turb.nz;
iseed = turb.iseed;
lsx = turb.lsx;
lsy = turb.lsy;
lsz = turb.lsz;
nfx = turb.nfx;
nfy = turb.nfy;
nfz = turb.nfz;
velvecs = turb.velvecs;
startup = isempty(velvecs);

if(startup)
disp('initializing random numbers');
% initialize random number generator
rng(iseed,'twister' );
end

I = -nfx : nfx;
J = -nfy+1 : ny + nfy;
K = -nfz+1 : nz + nfz;

ni = length(I);
nj = length(J);
nk = length(K);


if (startup)
    rx = rand([ni,nj,nk]);
    ry = rand([ni,nj,nk]);
    rz = rand([ni,nj,nk]);
    rfill = zeros([nj,nk]);
    
    rx = 2. * rx - 1.0;
    ry = 2. * ry - 1.0;
    rz = 2. * rz - 1.0;
else
    % Shift random fields:
    rx = turb.rx;
    ry = turb.ry;
    rz = turb.rz;
    rfill = rand([nj,nk]);
    rfill = 2. * rfill - 1.0;
    rx(2:end,:,:) = rx(1:end-1,:,:);
    rx(1,:,:) = rfill;
    
    rfill = rand([nj,nk]);
    rfill = 2. * rfill - 1.0;
    ry(2:end,:,:) = ry(1:end-1,:,:);
    ry(1,:,:) = rfill;
    
    rfill = rand([nj,nk]);
    rfill = 2. * rfill - 1.0;
    rz(2:end,:,:) = rz(1:end-1,:,:);
    rz(1,:,:) = rfill;
end
velvecs = vfromr(rx, ry, rz, nfx, nfy, nfz, ny, nz, lsx,lsy, lsz);

turb.rx = rx;
turb.ry = ry;
turb.rz = rz;
turb.velvecs = velvecs;

return
end
%	  
function velvecs = vfromr(rx, ry, rz, nfx, nfy, nfz, ny, nz, lsx, lsy,lsz)

bijk = filco(lsx, lsy, lsz, nfx, nfy, nfz);
velvecs.velx = zeros(ny,nz);
velvecs.vely = zeros(ny,nz);
velvecs.velz = zeros(ny,nz);

kp = (-nfz:nfz)+nfz+1;
jp = (-nfy:nfy)+nfy+1;
ip = (-nfx:nfx)+nfx+1;
        
%pool = parpool;

for k = 1:nz
    for j = 1:ny

        velvecs.velx(j,k) =  ...
            sum(sum(sum(bijk(ip,jp,kp) .* rx(ip, j+jp-1, k+kp-1))));
        velvecs.vely(j,k) = ...
            sum(sum(sum(bijk(ip,jp,kp) .* ry(ip, j+jp-1, k+kp-1))));
        velvecs.velz(j,k) = ...
            sum(sum(sum(bijk(ip,jp,kp) .* rz(ip, j+jp-1, k+kp-1))));
%         dumx(j) =  ...
%             sum(sum(sum(bijk(ip,jp,kp) .* rx(ip, j+jp-1, k+kp-1))));
%         dumy(j) = ...
%             sum(sum(sum(bijk(ip,jp,kp) .* ry(ip, j+jp-1, k+kp-1))));
%         dumz(j) = ...
%             sum(sum(sum(bijk(ip,jp,kp) .* rz(ip, j+jp-1, k+kp-1))));
    end
end
return
end  
      
function f = ffunc(k, n)
%! The actual filtering function. Notice that it's a
%! Gaussian.
f = exp(-(pi*k*k / (2.*n*n)));
return
end

function f = filco(lsx, lsy, lsz, nfx, nfy, nfz)
%!Returns the array of filter coefficients.
%! lsx, etc. are the length scales;
%! nfy, etc. are the filter widths.
%!The 1d coefficients that are multiplied together to form
%!the filco array:
I = -nfx : nfx;
J = -nfy : nfy;
K = -nfz : nfz;

ni = length(I);
nj = length(J);
nk = length(K);

bx(1:ni) = 0.0;
by(1:nj) = 0.0;
bz(1:nk) = 0.0;

f = zeros(ni,nj,nk);

s = 0.0;
for k = I
    s = s + ffunc(k, lsx);
end
s = sqrt(s);
for k = I
    bx(k+nfx+1) = ffunc(k, lsx) / s ;
end

s = 0.0;
for k = J
    s = s + ffunc(k, lsy);
end
s = sqrt(s);
for k = J
by(k+nfy+1) = ffunc(k, lsy) / s;
end

s = 0.0;
for k = K
    s = s + ffunc(k, lsz);
end
s = sqrt(s);
for k = K
    bz(k+nfz+1) = ffunc(k, lsz) / s;
end

for k = K+nfz+1
for j = J+nfy+1
for i = I+nfx+1
f(i,j,k) = bx(i)*by(j)*bz(k);
end
end
end

return
end

function f=periodic_ijk(f)

f(1,:,:) = (f(2,:,:) + f(end-1,:,:))*0.5;
f(end,:,:) = f(1,:,:);

f(:,1,:) = (f(:,2,:) + f(:,end-1,:))*0.5;
f(:,end,:) = f(:,1,:);

f(:,:,1) = (f(:,:,2) + f(:,:,end-1))*0.5;
f(:,:,end) = f(:,:,1);

return
end
       
