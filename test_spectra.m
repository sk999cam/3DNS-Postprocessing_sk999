
ni=64*2;
nj=64*4;
nk=64*5;

Lx = 1;
Ly = 1;
Lz = 1;

x = linspace(0,Lx,ni);
y = linspace(0,Ly,nj);
z = linspace(0,Lz,nk);



N = (ni*nj*nk);

[Y,X,Z]=meshgrid(y,x,z);

nmodes = 10;
f = zeros(ni,nj,nk);
for n=1:nmodes
f = f + (sin(Y*n*pi) + sin(X*n*pi) + sin(Z*n*pi))/n;
end

v = fftn(f);

Ek = v.*conj(v)/N; 

% sum of Ek = sum of f^2;
Total_Energy_time = sum(sum(sum(f.*f)))
Total_Energy_freq = sum(sum(sum(Ek)))

Ek_x = squeeze(sum(sum(Ek,2),3));
Ek_y = squeeze(sum(sum(Ek,1),3));
Ek_z = squeeze(sum(sum(Ek,1),2));
% 
% Ek_x = squeeze(mean(mean(Ek,2),3));
% Ek_y = squeeze(mean(mean(Ek,1),3));
% Ek_z = squeeze(mean(mean(Ek,1),2));

k_x = linspace(0,(ni-1)/Lx,ni);
k_y = linspace(0,(nj-1)/Ly,nj);
k_z = linspace(0,(nk-1)/Lz,nk);


loglog(k_x,Ek_x,k_y,Ek_y,k_z,Ek_z);