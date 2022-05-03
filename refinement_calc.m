clear
close all

base_folder = cd;
temp_slash = '/'; if ispc, temp_slash = '\'; end

case_name = 'r150_cwl90_3d';
blk = read_grid(case_name);

RANS_folder = [base_folder temp_slash 'RANS' temp_slash 'r150_cwl90' temp_slash 'Re200k_data'];

% Y_shear_data = readmatrix([RANS_folder temp_slash 'y-shear']);
% N = size(Y_shear_data,1)-1;
% points = [N:-1:135 1:134];
% Y_shear_data = Y_shear_data(points,:);
% [~,iLE] = min(Y_shear_data(:,1));
% Y_shear_data = circshift(Y_shear_data,(1-iLE),1);
% points = circshift(points,(1-iLE));
% [~,iTE] = max(Y_shear_data(:,1));
% points = points(N:-1:iTE);
% points2 = points(1:iTE);
% x = Y_shear_data(N:-1:iTE,1);
% x2 = Y_shear_data(1:iTE);
Y_shear_data = readmatrix([RANS_folder temp_slash 'y-shear']);
shear_data = readmatrix([RANS_folder temp_slash 'wall-shear']);
p_surf_data = readmatrix([RANS_folder temp_slash 'p_surf']);
density_data = readmatrix([RANS_folder temp_slash 'density']);
y_plus_data = readmatrix([RANS_folder temp_slash 'yplus']);
visc_data = readmatrix([RANS_folder temp_slash 'mu']);

%%
% Y_shear = Y_shear_data(points,2);
% shear = shear_data(points,2);
% density = density_data(points,2);
Y_shear_data = Y_shear_data(1:end-1,:);
shear_data = shear_data(1:end-1,:);
p_surf_data = p_surf_data(1:end-1,:);
density_data = density_data(1:end-1,:);
y_plus_data = y_plus_data(1:end-1,:);
visc_data = visc_data(1:end-1,:);

x_grid = [];
y_grid = [];
o_blocks = [4 6 9 5];
flip = [0, 0, 1, 1];
for i=1:length(o_blocks)
    x_tmp = blk{o_blocks(i)}.x(:,end:-1:end-1);
    y_tmp = blk{o_blocks(i)}.y(:,end:-1:end-1);
    if flip(i) == 1
        x_tmp = x_tmp(end:-1:1,:);
        y_tmp = y_tmp(end:-1:1,:);
    end
    x_grid = [x_grid x_tmp'];
    y_grid = [y_grid y_tmp'];
end
x = x_grid(1,:);
[~,iLE] = min(x);
[~,iTE] = max(x);
%x = x(iLE:iTE);
x_grid = x_grid(:,iLE:iTE);
y_grid = y_grid(:,iLE:iTE);
figure(2)
plot(x_grid(1,:),y_grid(1,:))
hold on
plot(x_grid(2,:),y_grid(2,:))
axis equal

for i = 1:size(x_grid,2)
    offset(i) = sqrt((x_grid(2,i)-x_grid(1,i))^2+(y_grid(2,i)-y_grid(1,i))^2);
end


for i=1:size(x_grid,2)
    [~,j] = min(abs(Y_shear_data(:,1)/0.02 - x_grid(1,i)));
    tau_w(i) = shear_data(j,2);
    rho(i) = density_data(j,2);
    Ps(i) = p_surf_data(j,2);
    yplus(i) = y_plus_data(j,2);
    mu(i) = visc_data(j,2);
end

nu = mu./rho;

Ufric = sqrt(tau_w./rho);
Yplus = 0.5*0.02*offset.*Ufric./nu; % Cell centroid at half offset

figure(4)
plot(x_grid(1,:),Yplus)
hold on
%plot(x,yplus)

[~,iLE] = min(x);
[~,iTE] = max(x);

P0 = 74575.2;
ga = 1.4;

Misen = real(sqrt((2/(ga-1))*((P0./Ps).^((ga-1)/ga) - 1)));
figure(5)
plot(x_grid(1,:),Misen)
write = false;
if write
    fid = fopen([base_folder temp_slash 'RANS_loading.txt'],'w+');
    fprintf(fid,'%7.5f\t%7.5f\n',[x(iLE:iTE); Misen]);
    fclose(fid);
end

%%

fid = fopen('RANS_Ufric_data_Re200k','w+');
fprintf(fid,'%f %f %f %e\n',[x_grid(1,:); tau_w; rho; mu]);
