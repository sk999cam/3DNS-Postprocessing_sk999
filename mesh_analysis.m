function mesh_analysis(blk, skip)

if nargin<2, skip=8; end

temp_slash = '/'; if ispc, temp_slash = '\'; end

C = colororder;

NB = length(blk.x);
try
    npp = blk.npp;
catch
    npp = 1;
end

nij_pts = 0;
figure(1)
hold on
for i=1:NB
    x = blk.x{i};
    y = blk.y{i};
    x = x(:,:,1);
    y = y(:,:,1);
    procdims(i,:) = size(x)/npp;
    ni = size(x,1);
    nj = size(x,2);
    nij_pts = nij_pts + ni*nj;
    nk = blk.nk;
    for j=[1:skip:ni ni]
        if (j==1) || (j==ni)
            plot(x(j,:),y(j,:),'r','LineWidth',1)
        else
            plot(x(j,:),y(j,:),'k')
        end
    end
    for j=[1:skip:nj nj]
        if (j==1) || (j==nj)
            plot(x(:,j),y(:,j),'r','LineWidth',1)
        else
            plot(x(:,j),y(:,j),'k')
        end
    end
%     plot(x(1,:),y(1,:),'r','LineWidth',1)
%     plot(x(end,:),y(end,:),'r','LineWidth',1)
%     plot(x(:,1),y(:,1),'r','LineWidth',1)
%     plot(x(:,end),y(:,end),'r','LineWidth',1)
end

axis equal
%%

% data = readmatrix([case_name temp_slash 'probe.txt']);
% for i=1:size(data,1)
%     probe{i}.nb = data(i,1);
%     probe{i}.ni = data(i,2);
%     probe{i}.nj = data(i,3);
%     probe{i}.nk = data(i,4);
% end
% 
% 
% for i=1:length(probe)
%     probe_plot(i) = scatter(blk{probe{i}.nb}.x(probe{i}.ni, probe{i}.nj),blk{probe{i}.nb}.y(probe{i}.ni, probe{i}.nj));
%     probe_legend(i) = string(sprintf('Probe %d',i));
% end
% 
% trip.nb = 6;
% trip.i = 20;
% trip.j = size(blk{trip.nb}.x,2);
% trip.x = blk{probe{i}.nb}.x(probe{i}.ni, probe{i}.nj);
% trip.y = blk{probe{i}.nb}.y(probe{i}.ni, probe{i}.nj);
% trip
% 
% probe_plot(length(probe)+1) = scatter(blk{trip.nb}.x(trip.i,trip.j),blk{trip.nb}.y(trip.i,trip.j));
% probe_legend(length(probe)+1) = string(sprintf('Trip location'));
% 
% 
% legend(probe_plot,probe_legend)
% 
% set(gca,'FontSize',14)
% 
% axis equal

x = blk.x{6}(:,1,1);
y = blk.y{6}(:,1,1);
dx_surf = (x(end)-x(1))/length(x);


%%
figure(2)
hold on
i=3;
skip = 8;
x = blk.x{i};
y = blk.y{i};
x = x(:,:,1);
y = y(:,:,1);
ni = size(x,1);
nj = size(x,2);



for j=[1 skip:skip:ni]
    if (j==1) || (j==ni)
        plot(x(j,:),y(j,:),'r','LineWidth',1)
    else
        plot(x(j,:),y(j,:),'k')
    end
end
for j=[1 skip:skip:nj]
    if (j==1) || (j==nj)
        plot(x(:,j),y(:,j),'b','LineWidth',1)
    else
        plot(x(:,j),y(:,j),'k')
    end
end

%% Cell aspect ratio

figure(3)
hold on
for i=1:NB
    x = blk.x{i};
    y = blk.y{i};
    x = x(:,:,1);
    y = y(:,:,1);

    dxi = x(2:end,:)-x(1:end-1,:);
    dxi(end+1,:) = dxi(end,:);
    dxj = x(:,2:end)-x(:,1:end-1);
    dxj(:,end+1) = dxj(:,end);
    
    dyi = y(2:end,:)-y(1:end-1,:);
    dyi(end+1,:) = dyi(end,:);
    dyj = y(:,2:end)-y(:,1:end-1);
    dyj(:,end+1) = dyj(:,end);
    
    dsi = sqrt(dxi.^2+dyi.^2);
    dsj = sqrt(dxj.^2+dyj.^2);
    
    aspect_ratio{i} = max(dsi./dsj, dsj./dsi);
    pcolor(blk.x{i},blk.y{i},aspect_ratio{i})
    contour(blk.x{i},blk.y{i},aspect_ratio{i},[11 13 15],'-k')
end
shading('interp')

axis equal
cb = colorbar;

%% O grid analysis

% A = readmatrix('RANS_Ufric_data_Re200k');
% x_RANS = A(:,1);
% [x_RANS, inds, ~] = unique(x_RANS);
% c_RANS = 0.02;
% tau_w = A(inds,2);
% rho = A(inds,3);
% mu = A(inds,4);
% nu = mu./rho;

% Ufric = sqrt(tau_w./rho);

x_grid = [];
y_grid = [];
if NB == 12
    o_blocks = [4 6 9 5];
    flip = [0, 0, 1, 1];
elseif NB == 9
    o_blocks = [3 5 7 4];
    flip = [0, 0, 1, 1];
end
    
for i=1:length(o_blocks)
    x_tmp = blk.x{o_blocks(i)}(:,end:-1:1);
    y_tmp = blk.y{o_blocks(i)}(:,end:-1:1);
    if flip(i) == 1
        x_tmp = x_tmp(end:-1:1,:);
        y_tmp = y_tmp(end:-1:1,:);
    end
    x_grid = [x_grid; x_tmp];
    y_grid = [y_grid; y_tmp];
end
x = x_grid(:,1);
y = y_grid(:,1);
[~,iLE] = min(x);
[~,iTE] = max(x);
iMid = round(0.5*(iLE+iTE));
x = x(iLE:iTE);
tmp = iLE:iTE;
[x, inds, ~] = unique(x);
inds = tmp(inds);
dy = sqrt((x_grid(inds,2)-x_grid(inds,1)).^2 + (y_grid(inds,2) - y_grid(inds,1)).^2);
dx = sqrt((x_grid(inds(2:end),1) - x_grid(inds(1:end-1),1)).^2 + (y_grid(inds(2:end),1) - y_grid(inds(1:end-1),1)).^2);
span = 0.1;
dz = span/(nk-1);

x_cut = x_grid(iMid,:);
y_cut = y_grid(iMid,:);
dy_bl = sqrt((x_cut(2:end)-x_cut(1:end-1)).^2 + (y_cut(2:end) - y_cut(1:end-1)).^2);
expansion_r = dy_bl(2:end)./dy_bl(1:end-1);


% Ufric = interp1(x_RANS,Ufric,x,'pchip','extrap');
% 
% nu = interp1(x_RANS,nu,x,'pchip','extrap');
% nu = nu/c_RANS;
% Xplus = dx.*Ufric(1:end-1)./nu(1:end-1);
% Yplus = dy.*Ufric./nu;
% Zplus = dz.*Ufric./nu;
% 
% figure(3)
% plot(x,Yplus)
% ylabel('First point y^{+}')
% xlabel('x/c')
% xlim([0 1])
% set(gca,'FontSize',16)
% grid on
% 
% figure(4)
% plot(x(1:end-1),Xplus)
% hold on
% plot(x,Zplus)
% xlabel('x/c')
% xlim([0 1])
% legend('\Delta x^{+}','\Delta z^{+}')
% set(gca,'FontSize',16)
% grid on

figure(5)
plot(2:length(expansion_r)+1,expansion_r)
xlabel('y point')
ylabel('O grid expansion ratio')
set(gca,'FontSize',16)

%%

nij_pts
Nprocs = nk*sum(prod(procdims,2))/npp
Npts = nij_pts*nk
end