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
[~, ps_ind_foil] = max(1./pr_foil'.*(x_foil>0.4 & x_foil<0.6));
theta_ps_foil = th_foil(ps_ind_foil);
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

%%

x_ch = g300M.meanFlow.xSurf;
pr_ch = g300M.meanFlow.blPr;
th_ch = g300M.meanFlow.theta;
M_ch = g300M.meanFlow.Msurf;
Re_ch = g300M.meanFlow.Re_theta;

[~, ps_ind_ch] = max(1./pr_ch'.*(x_ch>0.35 & x_ch<0.6));
theta_ps_ch = th_ch(ps_ind_ch);
i = g300M.meanFlow.x2ind(0.2);
while M_ch(i)>1
    i = i+1;
end
xShock_ch = x_ch(i);

%%


f1 = figure;
ax = axes(f1);
plot(ax,(x_foil-xShock_foil)/theta_ps_foil,pr_foil);
hold on
plot(ax,(x_ch-xShock_ch)/theta_ps_ch,pr_ch);

%%

f2 = figure;
ax2 = axes(f2);
plot(ax2,(x_foil-xShock_foil)/theta_ps_foil,Re_foil);
hold on
plot(ax2,(x_ch-xShock_ch)/theta_ps_ch,Re_ch);