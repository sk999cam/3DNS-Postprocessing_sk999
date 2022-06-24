function blade_mesher()
%clear
% profile inputs
% N = 5000;
% chi(1) = 70.0;
% chi(2) = 50.0;
% tLE = 0.03;

% geometry inputs
% pitch = 0.75;
% stag = (chi(1) + chi(2))*0.5;

% initial mesh inputs
npp(1) = 4;
Lup = 0.75;
Ldn = 1.0;
Lo = 0.05;
ywall(1) = 0.005;
msmooths(1) = 500;

% final mesh inputs
npp(2) = 28;
refine_fac = 4;
ywall(2) = 0.0005;
msmooths(2) = 50;


% solver inputs
% boundary conditions
bcs.Toin = 300.0;
bcs.Poin = 100000.0;
bcs.pexit = 0.99625e5;
bcs.vin = 50.0;
bcs.alpha = 40.0;
bcs.gamma = 0.0;
bcs.aturb = -1.0;
bcs.lturb = 10*ywall(2);
bcs.ilength = 500;
bcs.radprof = 0;
bcs.g_z = 0.0;
bcs.twall = -1.0;
bcs.cax = 1.0;


% gas properties
gas.gamma=1.4;
gas.cp = 1005.0;
gas.mu_ref = 100*18.2e-6;
gas.mu_tref = 273.0;
gas.mu_cref = 110.4;
gas.pr = 0.72;


% create cicular arc profile
%[xprof,yprof]=blade_profile(chi,tLE,N);

% read in profile
[xprof,yprof,pitch,stag]=read_profile('harrison.txt');  

% get surface distance
sprof = curve_length(xprof,yprof);

% create 9-block topology and initial coarse mesh
[blk,next_block,next_patch,corner] = blade_topology(xprof,yprof,pitch,Lup,Ldn,Lo,stag,npp(1),ywall(1),msmooths(1));

% create final refined mesh
[blk] = mesh_refinement(blk,refine_fac,npp(2));
[blk] = mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths(2),xprof,yprof,ywall(2));

% write 3DNS case files
write_case('harrison',blk,next_block,next_patch,corner,bcs,gas)

% plot profile, topology and mesh

figure(1)
plot(xprof,yprof,'-r.'), axis equal
hold on

NB = length(blk);
for ib=1:NB
xnew=blk{ib}.x;
ynew=blk{ib}.y;
%pcolor(xnew,ynew,xnew), hold on
%pcolor(xnew,ynew+pitch,xnew)
plot(xnew,ynew+pitch,'k'), hold on
plot(xnew',ynew'+pitch,'k')
axis equal
end
   
% plot topology
for ib=1:NB
xnew=blk{ib}.x;
ynew=blk{ib}.y;
plot(xnew(1,:),ynew(1,:),'k'), hold on
axis equal
plot(xnew(end,:),ynew(end,:),'k')
plot(xnew(:,end),ynew(:,end),'k'), hold on
plot(xnew(:,1),ynew(:,1),'k'), hold on
end
% 
% 
return