clear
close all
profile = fullfile('r150_aerofoil.txt');



[xprof,yprof,pitch,stag] = read_profile(profile,true);
msmooths=50;
npp=6;
pitch=2;
bound_angle = 5.0;
Lup=2.0;
Ldn=1.0;
Lo=0.05;
ywall=Lo/npp;

[blk,next_block,next_patch,corner] = single_blade_topology(xprof,yprof,pitch,bound_angle,Lup,Ldn,Lo,stag,npp,ywall,msmooths,0.1);

NB = length(blk);

% Intermediate mesh imputs
npp = 36; % points in i,j direction for each block are set to a multiple of npp- this ensures load balancing when running on many cores
refine_fac = 6.0; % the grid dimensions are scaled by this number in both directions
ywall = 0.0002; % the final near wall cell height
msmooths = 100; % number of iterations for the poisson solver

blk=mesh_refinement(blk,refine_fac,npp);

blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall);



for i=1:NB
    int_blk_dims(i,:) = size(blk{i}.x);
end
disp('Intermediate grid block dims:')
int_blk_dims

npp = 28;
refine_fac = 4.0;
ywall = 0.000075;
msmooths = 20;

blk=mesh_refinement(blk,refine_fac,npp);
sprintf("BL points: %d", length(blk{6}.x(1,:)));
blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall,1);


blk = thin_blocks(blk, [1,2,3], [2,2,2], 2, npp);
blk = thin_blocks(blk, [10,11,12], [1,1,1], 2, npp);
blk = thin_blocks(blk, [3,8,12], [3,3,3], 2, npp);
blk = thin_blocks(blk, [2,7,11], [4,4,4], 2, npp);



fig = figure;
hold on
for i=1:NB
    xnew = blk{i}.x;
    ynew = blk{i}.y;
    plot(xnew,ynew,'k');
    plot(xnew',ynew','k');
end
axis equal

for i=1:NB
    blk_dims(i,:) = size(blk{i}.x);
end
disp('Final grid block dims:')
blk_dims
disp('Final grid proc dims:')
blk_dims/npp
sprintf('Total points: %d', sum(prod(blk_dims,2)))
sprintf('Total ij cores: %d', sum(prod(blk_dims/npp,2)))

for i=1:NB
    blk2.x{i}=blk{i}.x;
    blk2.y{i}=blk{i}.y;
    blk2.nk{i}=1;
end
blk2.oblocks = [4 6 9 5];
blk2.oblocks_flip = [0 0 1 1];
blk2.span = 0.1;
blk2.io_surfaces.blks = [1 2 2 7 11 11 10 12 12 8 3 3];
blk2.io_surfaces.types = [1 1 3 3 3 2 2 2 4 4 4 1];
blk2.blockdims = blk_dims;

blk = blk2;

