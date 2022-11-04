clear
close all

load DNS_loading.mat

run = 18;
run_fluent = true;

params.amp = 0.0045;
params.x1 = 0.4300;
params.x2 = 0.0650;
params.xTrip = 0.110;
params.r = 1.1000;
params.theta = 80;
params.tmid = 0.5000;

cp = 1005;
gam = 1.4;
p0 = 74575.2;

folder = ['E:\MATLAB\RANS_optimisation\run' num2str(run)];

p = mesh_and_run_Fluent(params,folder,run_fluent);


%%

run = run;
folder = ['E:\MATLAB\RANS_optimisation\run' num2str(run)];
fid = fopen(fullfile(folder,'loading.txt'),'r');

%load(fullfile(folder,'params_loading.mat'))

fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);

data = fscanf(fid,'%f %f',[2 Inf]);

x = data(1,:);
p = data(2,:);
M = sqrt(((p0./p).^((gam-1)/gam) - 1)*(2/(gam-1)));

figure
plot(x_DNS,M_DNS)
hold on
plot(x,M)
set(gca,'FontSize',12)
xlabel('x/c')
ylabel('M_{isen}')
grid on
legend('DNS','RANS')
save(fullfile(folder,'params_loading.mat'),'x','M','p','params')


%%
run2=5;
folder2 = ['E:\MATLAB\RANS_optimisation\run' num2str(run2)];
load(fullfile(folder2,'params_loading.mat'))
plot(x,M)

%%


