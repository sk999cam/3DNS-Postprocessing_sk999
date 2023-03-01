clear
close all
proffile = '../isolated_aerofoils/generic/hst_aerofoil.txt';
% proffile = '../isolated_aerofoils/r150_cwl_90/r150_aerofoil.txt';

[xprof,yprof,pitch,stag] = read_profile(proffile,true);
msmooths=50;
npp=6;
pitch=2;
bound_angle = 5.0;
Lup=2.0;
Ldn=1.0;
Lo=0.05;
ywall=Lo/npp;

[blk,next_block,next_patch,corner]=single_blade_topology(xprof,yprof,pitch,bound_angle,Lup,Ldn,Lo,stag,npp,ywall,msmooths,0.1);

NB = length(blk);

% Intermediate mesh imputs
npp = 12; % points in i,j direction for each block are set to a multiple of npp- this ensures load balancing when running on many cores
refine_fac = 2.0; % the grid dimensions are scaled by this number in both directions
ywall = 0.0008; % the final near wall cell height
msmooths = 100; % number of iterations for the poisson solver

blk=mesh_refinement(blk,refine_fac,npp);

blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall);



for i=1:NB
    int_blk_dims(i,:) = size(blk{i}.x);
end
disp('Intermediate gris block dims:')
int_blk_dims

npp = 48;
refine_fac = 8.0;
ywall = 0.0001;
msmooths = 20;

blk=mesh_refinement(blk,refine_fac,npp);
sprintf("BL points: %d", length(blk{6}.x(1,:)));
blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall,1);

blk = thin_blocks(blk, [1,2,3], [2,2,2], 2.8, 48);
blk = thin_blocks(blk, [10,11,12], [1,1,1], 2.0, 48);
blk = thin_blocks(blk, [3,8,12], [3,3,3], 2.0, 48);
blk = thin_blocks(blk, [2,7,11], [4,4,4], 2.0, 48);



fig = figure;
hold on
for i=1:NB
    xnew = blk{i}.x;
    ynew = blk{i}.y;
    plot(xnew,ynew,'k');
    plot(xnew',ynew','k');
end
axis equal
blktmp = blk;
clear blk
for i=1:NB
    blk_dims(i,:) = size(blktmp{i}.x);
    blk.x{i} = blktmp{i}.x;
    blk.y{i} = blktmp{i}.y;
end
disp('Final grid block dims:')
blk_dims
disp('Final grid proc dims:')
blk_dims/npp
sprintf('Total points: %d', sum(prod(blk_dims,2)))
sprintf('Total ij cores: %d', sum(prod(blk_dims/npp,2)))