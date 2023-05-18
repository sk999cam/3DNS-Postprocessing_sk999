clear
close all

%%

g2 = DNS_case('300k_grid2_tripped',[8 9 10 11 12]);
g2.readMeanFlows

%%

g300M = DNS_channel('M1-2_g300M',1);
g300M.readMeanFlow
g300M.meanFlow.iSmoothBL = true;

%%
x_foil = g2.meanFlow.xSurf;
inds = x_foil > 0.25 & x_foil < 0.95;
pr_foil = g2.meanFlow.blPr;
th_foil = g2.meanFlow.theta;
del_foil = g2.meanFlow.delStar;
[~, ps_ind_foil] = max(1./pr_foil'.*(x_foil>0.4 & x_foil<0.6));
theta_ps_foil = th_foil(ps_ind_foil);
del_ps_foil = del_foil(ps_ind_foil);
M_foil = g2.meanFlow.Msurf;
i = g2.meanFlow.x2ind(0.4);
while M_foil(i)>1
    i = i+1;
end
xShock_foil = x_foil(i);
x_foil = x_foil(inds);
M_foil = M_foil(inds);
th_foil = th_foil(inds);
pr_foil = pr_foil(inds);
Re_foil = g2.meanFlow.Re_theta(inds);
del_foil = del_foil(inds);

%%

x_ch = g300M.meanFlow.xSurf;
pr_ch = g300M.meanFlow.blPr;
th_ch = g300M.meanFlow.theta;
M_ch = g300M.meanFlow.Msurf;
Re_ch = g300M.meanFlow.Re_theta;
del_ch = g300M.meanFlow.delStar;

[~, ps_ind_ch] = max(1./pr_ch'.*(x_ch>0.35 & x_ch<0.6));
theta_ps_ch = th_ch(ps_ind_ch);
del_ps_ch = del_ch(ps_ind_ch);
i = g300M.meanFlow.x2ind(0.2);
while M_ch(i)>1
    i = i+1;
end
xShock_ch = x_ch(i);

%%


f1 = figure;
ax = axes(f1);
plot(ax,(x_foil-xShock_foil)/del_ps_foil,pr_foil, 'k','LineWidth',1);
hold on
plot(ax,(x_ch-xShock_ch)/del_ps_ch,pr_ch, 'b','LineWidth',1);
xlabel('$(x-x_{shock})/\delta^*_{pre-shock}$','Interpreter','latex');
ylabel('Turbulence productiom', 'Interpreter','latex')
pbaspect([3 1 1]);
grid on
set(gca, 'FontSize',12)
legend('Isolated aerofoil','Channel', 'Location','northwest');


%%

f2 = figure;
ax2 = axes(f2);
plot(ax2,(x_foil-xShock_foil)/del_ps_foil,Re_foil, 'k','LineWidth',1);
hold on
plot(ax2,(x_ch-xShock_ch)/del_ps_ch,Re_ch, 'b', 'LineWidth',1);

%%

figure();
a = g300M.meanFlow.blDevPlot('H_k',[],[1 3],[],'k')
a.LineWidth=2;
grid on
xlabel('$x/L$','Interpreter','latex')
ylabel('$H_k$','Interpreter','latex')
set(gca, 'FontSize',22)
pbaspect([5 1 1])
figure()
b = g300M.meanFlow.blDevPlot('cf',[],[0 0.005],[],'k');
b.LineWidth=2;
xlabel('$x/L$','Interpreter','latex')
pbaspect([5 1 1])
ylabel('$c_f$','Interpreter','latex')
grid on
set(gca, 'FontSize',22)


%%
figure()
ax = gca;
a = g300M.meanFlow.blDevPlot('blPr',ax,[],[],'k');
a.LineWidth = 2;
hold on
a = g300M.meanFlow.blDevPlot('blPr_eq',ax,[],[],'k:');
a.LineWidth = 2;
legend('DNS','Equilibrium');
set(gca, 'FontSize',26)
pbaspect([3 1 1]);
grid on
xlabel('$x/L$','Interpreter','latex')
ylabel('Turbulence productiom', 'Interpreter','latex')
ylim([0 0.004])
