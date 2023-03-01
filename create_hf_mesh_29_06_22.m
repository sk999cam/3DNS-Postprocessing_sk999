clear
close all
proffile = '../isolated_aerofoils/r150_cwl_90/r150_aerofoil.txt';

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
npp = 24; % points in i,j direction for each block are set to a multiple of npp- this ensures load balancing when running on many cores
refine_fac = 4.0; % the grid dimensions are scaled by this number in both directions
ywall = 0.0006; % the final near wall cell height
msmooths = 10; % number of iterations for the poisson solver

blk=mesh_refinement(blk,refine_fac,npp);

blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall);



for i=1:NB
    int_blk_dims(i,:) = size(blk{i}.x);
end
disp('Intermediate grid block dims:')
int_blk_dims

npp = 56;
refine_fac = 8.0;
ywall = 0.000075;
msmooths = 1;
nkprocs = 3;

blk=mesh_refinement(blk,refine_fac,npp);
sprintf("BL points: %d", length(blk{6}.x(1,:)));
blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall,1);

for i=1:NB
    fine_blk_dims(i,:) = size(blk{i}.x);
end
disp('Fine grid block dims:')
fine_blk_dims

%% 


blk = thin_blocks(blk, [1,2,3], [2,2,2], 5, npp);
blk = thin_blocks(blk, [10,11,12], [1,1,1], 4, npp);
blk = thin_blocks(blk, [3,8,12], [3,3,3], 5, npp);
blk = thin_blocks(blk, [2,7,11], [4,4,4], 5, npp);

%%



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
for i=1:NB
    blktmp.x{i}=blk{i}.x;
    blktmp.y{i}=blk{i}.y;
    blktmp.nk{i} = nkprocs*npp;
end
blk = blktmp;
blk.npp = npp;

fprintf('Total ij points: %d\n', sum(prod(blk_dims,2)))
fprintf('Total ij cores: %d\n', sum(prod(blk_dims/npp,2)))
fprintf('Total points: %dM\n', round(sum(prod(blk_dims,2))*blk.nk{1}/1e6))
fprintf('Total cores: %d\n', nkprocs*sum(prod(blk_dims/npp,2)))



