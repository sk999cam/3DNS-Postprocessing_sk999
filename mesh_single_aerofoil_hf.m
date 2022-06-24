clear
close all
proffile = '../isolated_aerofoils/r150_cwl_90/r150_aerofoil/r150_aerofoil.txt';

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
figure
hold on
for i=1:NB
    xnew = blk{i}.x;
    ynew = blk{i}.y;
    plot(xnew,ynew,'k');
    plot(xnew',ynew','k');
end
axis equal

