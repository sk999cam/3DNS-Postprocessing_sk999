clear
close all
case_name = 'r150_cwl';
g = ts_read_hdf5('ts_results/turbostream_r150_hst_1Re_trans_66_0_w_patches.hdf5');
gam=1.4;
cp = g{1,1}.av.cp;
R = cp*(gam-1)/gam;

base_folder = cd;
temp_slash = '/'; if ispc, temp_slash = '\'; end

coords = readmatrix('r150_aerofoil.txt');
x = coords(:,1);
y = coords(:,2);
%plot(x,y)
%axis equal

c = max(x) - min(x);
x = x/c;
y = y/c;
thickness = max(y)-min(y);
p0 = 74575.203;
T0 = 300;
Minf = 0.84;
Tinf = T0/(1+((gam-1)/2)*Minf^2);
mu0 = 1.8e-5;
mu_inf = mu0*(Tinf/288)^0.62;
Vinf = Minf*sqrt(gam*R*Tinf);
ro0 = p0/(R*T0);
ro_inf = ro0*((1+((gam-1)/2)*Minf^2)^(-1/(gam-1)));
Re_c = ro_inf*Vinf*c/mu_inf;

Re_init = 3.2e4;
visc_mult = Re_c/(Re_init*c);



g = ts_secondary(g);
h.fig = 2;
h.sub = 111;
type.view = 'B2B';
%type.mesh = 'on';
%type.res = 'fine';
ts_plot_structured(g,'M',[-1,3,-1],type,[],[0:4],h)
ts_P_surf = g{1,1}.P(:,3,1)';
ts_x_surf = g{1,1}.x(:,3,1)';
ts_y_surf = g{1,1}.rt(:,3,1)';
[~,iLE] = min(ts_x_surf);
x_c = (ts_x_surf-ts_x_surf(iLE))/(ts_x_surf(1) - ts_x_surf(iLE));
x_c = x_c(iLE:end);
s_surf = 0;
for i=2:length(ts_x_surf)
    s_surf(i) = s_surf(i-1)+sqrt((ts_x_surf(i) - ts_x_surf(i-1))^2+(ts_y_surf(i) - ts_y_surf(i-1))^2);
end
Misen = sqrt((2/(gam-1))*((ts_P_surf/p0).^(-(gam-1)/gam)-1));
s_surf = s_surf(iLE:end);
Misen = Misen(iLE:end);

figure(3)
loading = readmatrix([case_name '/aerofoil_loading.txt']);
fluent_loading = readmatrix([base_folder temp_slash 'RANS' temp_slash 'r150_cwl90' temp_slash 'RANS_loading.txt']);

hold on
plot(loading(:,1),loading(:,2))
plot(x_c,Misen)
plot(fluent_loading(:,1),fluent_loading(:,2))
grid on
ylabel('M_{isen}')
xlabel('x/c')
legend('3DNS: Re = 32k','RANS - Turbostream','RANS - Fluent: SA, Re = 200k','location','southeast')
set(gca,'FontSize',16)
figure(4)

ts_P_inlet = squeeze(g{4,1}.P(1,3,:))';
ts_x_inlet = squeeze(g{4,1}.x(1,3,:))';
ts_y_inlet = squeeze(g{4,1}.rt(1,3,:))';
ts_M_inlet = squeeze(g{4,1}.M(1,3,:))';
ts_V_inlet = squeeze(g{4,1}.V(1,3,:))';

ts_P_upstream = squeeze(g{4,1}.P(8,3,:))';
ts_x_upstream = squeeze(g{4,1}.x(8,3,:))';
ts_y_upstream = squeeze(g{4,1}.rt(8,3,:))';
ts_M_upstream = squeeze(g{4,1}.M(8,3,:))';
ts_V_upstream = squeeze(g{4,1}.V(8,3,:))';
hold on
plot(ts_V_inlet/272.87,ts_y_inlet,'b:')
plot(ts_V_upstream/272.87,ts_y_inlet,'b-')

plot(ts_M_inlet,ts_y_inlet,'r:')
plot(ts_M_upstream,ts_y_inlet,'r-')


