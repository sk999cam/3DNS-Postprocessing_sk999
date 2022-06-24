function [blk,next_block,next_patch,corner] = blade_topology(xprof,yprof,pitch,Lup,Ldn,Lo,stag,npp,ywall,msmooths)

% create 9-block topology with o-grid
% returns a coarse mesh defining an initial topology

% get LE and TE points

[xLE,iLE] = min(xprof);
yLE = yprof(iLE);

[xTE,iTE] = max(xprof);
yTE = yprof(iTE);

% O-grid boundary
N = length(xprof);
for i=1:N
im1=i-1;
ip1=i+1;

if(i==1)
    im1 = N-1;
end
if(i==N)
    ip1 = 2;
end

dx = xprof(ip1)-xprof(im1);
dy = yprof(ip1)-yprof(im1);
ds = sqrt(dx*dx + dy*dy);   

xo(i) = xprof(i) - Lo*dy/ds;
yo(i) = yprof(i) + Lo*dx/ds;
    
end

xo(N) = (xo(1)+xo(N))*0.5;
yo(N) = (yo(1)+yo(N))*0.5;
xo(1) = xo(N);
yo(1) = yo(N);

%plot bounding box
%rotate by stagger about LE to find bounding box;
[~,ile]=min(xo);
xdat = xLE;
ydat = yLE;
xr = xo-xdat;
yr = yo-ydat;
if(stag>30*pi/180)
    stag=30*pi/180;
end
stag = -stag;
xr_now = xr.*cos(stag) - yr.*sin(stag);
yr_now = xr.*sin(stag) + yr.*cos(stag);

xr = xr_now;
yr = yr_now;

minx = min(xr);
maxx = max(xr);
miny = min(yr);
maxy = max(yr);


% rotate back
xb1 = minx*cos(stag) - miny*sin(-stag) + xdat;
xb2 = minx*cos(stag) - maxy*sin(-stag) + xdat;
xb3 = maxx*cos(stag) - maxy*sin(-stag) + xdat;
xb4 = maxx*cos(stag) - miny*sin(-stag) + xdat;

yb1 = minx*sin(-stag) + miny*cos(stag) + ydat;
yb2 = minx*sin(-stag) + maxy*cos(stag) + ydat;
yb3 = maxx*sin(-stag) + maxy*cos(stag) + ydat;
yb4 = maxx*sin(-stag) + miny*cos(stag) + ydat;


% plot([xb1 xb2 xb3 xb4 xb1],[yb1 yb2 yb3 yb4 yb1],'k')
% text(xb1,yb1,'1b')
% text(xb2,yb2,'2b')
% text(xb3,yb3,'3b')
% text(xb4,yb4,'4b')

x1 = xLE - Lup;
y1 = yb1;

x2 = x1;
y2 = y1 + pitch/2;

x3 = x1;
y3 = y1 - pitch/2;

x4 = xTE + Ldn;
y4 = yb4;

x5 = x4;
y5 = y4 + pitch/2;

x6 = x4;
y6 = y4 - pitch/2;


%find point on o-bound nearest to xb2,yb2
d = (xo-xb2).^2 + (yo-yb2).^2;
[dmin2,imin2] = min(sqrt(d));

x7 = xo(imin2);
y7 = yo(imin2);

x8 = x7;
y8 = y7 - pitch;


%find point on o-bound nearest to xb1,yb1
d = (xo-xb1).^2 + (yo-yb1).^2;
[dmin3,imin3] = min(sqrt(d));

x9 = xo(imin3);
y9 = yo(imin3);


%find point on o-bound nearest to xb3,yb3
d = (xo-xb3).^2 + (yo-yb3).^2;
[dmin5,imin5] = min(sqrt(d));

x10 = xo(imin5);
y10 = yo(imin5);

x11 = x10;
y11 = y10 - pitch;

% %find point on o-bound nearest to xb4,yb4
d = (xo-xb4).^2 + (yo-yb4).^2;
[dmin6,imin6] = min(sqrt(d));


x12 = xo(imin6);
y12 = yo(imin6);

x13 = x12;
y13 = y12+pitch;

% % plot topology vertices
% plot([x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13],[y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13],'rs');
% text(x1,y1,'1')
% text(x2,y2,'2')
% text(x3,y3,'3')
% text(x4,y4,'4')
% text(x5,y5,'5')
% text(x6,y6,'6')
% text(x7,y7,'7')
% text(x8,y8,'8')
% text(x9,y9,'9')
% text(x10,y10,'10')
% text(x11,y11,'11')
% text(x12,y12,'12')
% text(x13,y13,'13')

%
% get curve lengths to set mesh points based on smallest arc length (block
% 3 LE)
[xp,yp]=clip_curve(xprof,yprof,imin3,imin2,200);
% get curve length
s = curve_length(xp,yp);
sle = s(end);

% estimate profile surface length
sprof = curve_length(xprof,yprof);
ssurf = sprof(end)/2;



%
% Creat initial mesh to solve for block boundaries
%
% Block 1
%
ni=npp*ceil(Lup/sle);
nj=npp*ceil((y2-y1)/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

% extract part of o-grid
[xp,yp]=clip_curve(xo,yo,imin3,imin2,nj);

nic = 2;
njc = nj;

if(yp(njc)<yp(1)) %flip if needed
yp = yp(njc:-1:1);
xp = xp(njc:-1:1);
end

xp = repmat(xp,[nic 1]);
yp = repmat(yp,[nic 1]);

xp(1,:) = linspace(x1,x2,njc);
yp(1,:) = linspace(y1,y2,njc);

% create arrays for interpolation
i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);

% use 2-d interpolation to generate base grid
x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{1}.x = x;
blk{1}.y = y;

clear xp yp

%
% Block 2
%
ni=npp*ceil(Lup/sle);
nj=npp*ceil((y9-y8)/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

nic = 2;
njc = nj;

xp(1,:) = linspace(x3,x1,njc);
yp(1,:) = linspace(y3,y1,njc);

xp(nic,:) = linspace(x8,x9,njc);
yp(nic,:) = linspace(y8,y9,njc);

i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);

x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{2}.x = x;
blk{2}.y = y;


clear xp yp

%
% Block 3
%
ni=npp*ceil((y2-y1)/sle);
nj=npp*ceil(Lo/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp,yp]=clip_curve(xo,yo,imin3,imin2,ni);
[xp2,yp2]=clip_curve(xprof,yprof,imin3,imin2,ni);

nic = ni;
njc = 2;

if(yp(nic)<yp(1)) %flip
yp = yp(nic:-1:1);
xp = xp(nic:-1:1);
yp2 = yp2(nic:-1:1);
xp2 = xp2(nic:-1:1);
end

xp = repmat(xp',[1 njc]);
yp = repmat(yp',[1 njc]);

xp(:,2) = xp2';
yp(:,2) = yp2';

i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);


x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);


blk{3}.x = x;
blk{3}.y = y;



clear xp yp

%
% Block 4
%
ni=npp*ceil(ssurf/sle);
nj=npp*ceil(Lo/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp,yp]=clip_curve(xo,yo,imin3,imin6,ni);
[xp2,yp2]=clip_curve(xprof,yprof,imin3,imin6,ni);

nic = ni;
njc = 2;

if(xp(nic)<xp(1)) %flip
yp = yp(nic:-1:1);
xp = xp(nic:-1:1);
yp2 = yp2(nic:-1:1);
xp2 = xp2(nic:-1:1);
end

xp = repmat(xp',[1 njc]);
yp = repmat(yp',[1 njc]);

xp(:,njc) = xp2;
yp(:,njc) = yp2;


i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);

x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{4}.x = x;
blk{4}.y = y;


clear xp yp

%
% Block 5
%
ni=npp*ceil(ssurf/sle);
nj=npp*ceil(Lo/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp,yp]=clip_curve(xo,yo,imin2,imin5,ni);
[xp2,yp2,xp3,yp3]=clip_curve(xprof,yprof,imin2,imin5,ni);

if(mean(yp2)<mean(yp3))
xp2 = xp3;
yp2 = yp3;
end

nic = ni;
njc = 2;

if(xp(nic)<xp(1)) %flip
yp = yp(nic:-1:1);
xp = xp(nic:-1:1);
yp2 = yp2(nic:-1:1);
xp2 = xp2(nic:-1:1);
end

xp = repmat(xp',[1 njc]);
yp = repmat(yp',[1 njc]);

xp(:,njc) = xp2;
yp(:,njc) = yp2;


i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);


x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{5}.x = x;
blk{5}.y = y;


clear xp yp

%
% Block 6
%
ni=npp*ceil(ssurf/sle);
nj=npp*ceil((y9-y8)/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp2,yp2]=clip_curve(xo,yo,imin3,imin6,ni);
[xp,yp]=clip_curve(xo,yo,imin2,imin5,ni);
nic = ni;
njc = 2;

if(xp(nic)<xp(1)) %flip
yp = yp(nic:-1:1);
xp = xp(nic:-1:1);
end

if(xp2(nic)<xp2(1)) %flip
yp2 = yp2(nic:-1:1);
xp2 = xp2(nic:-1:1);
end

xp = repmat(xp',[1 njc]);
yp = repmat(yp'-pitch,[1 njc]);

xp(:,njc) = xp2;
yp(:,njc) = yp2;


i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);

x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{6}.x = x;
blk{6}.y = y;

clear xp yp

%
% Block 7
%
ni=npp*ceil((y5-y4)/sle);
nj=npp*ceil(Lo/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp,yp]=clip_curve(xo,yo,imin5,imin6,ni);
[xp2,yp2]=clip_curve(xprof,yprof,imin5,imin6,ni);

nic = ni;
njc = 2;

if(yp(nic)<yp(1)) %flip
yp = yp(nic:-1:1);
xp = xp(nic:-1:1);
yp2 = yp2(nic:-1:1);
xp2 = xp2(nic:-1:1);
end


xp = repmat(xp',[1 njc]);
yp = repmat(yp',[1 njc]);

xp(:,2) = xp2';
yp(:,2) = yp2';

i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);


x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{7}.x = x;
blk{7}.y = y;


clear xp yp

%
% Block 8
%
ni=npp*ceil(Ldn/sle);
nj=npp*ceil((y5-y4)/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

[xp,yp]=clip_curve(xo,yo,imin5,imin6,nj);
nic = 2;
njc = length(xp);

if(yp(njc)<yp(1)) %flip
yp = yp(njc:-1:1);
xp = xp(njc:-1:1);
end

xp = repmat(xp,[nic 1]);
yp = repmat(yp,[nic 1]);

xp(nic,:) = linspace(x4,x5,njc);
yp(nic,:) = linspace(y4,y5,njc);


i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);


x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);


blk{8}.x = x;
blk{8}.y = y;


clear xp yp

%
% Block 9
%
ni=npp*ceil(Ldn/sle);
nj=npp*ceil((y9-y8)/sle);

i=1:ni;
j=1:nj;
fi = (i-1)/(ni-1);
fj = (j-1)/(nj-1);

[J,I]=meshgrid(fj,fi);

nic = 2;
njc = nj;

xp(1,:) = linspace(x11,x12,njc);
yp(1,:) = linspace(y11,y12,njc);

xp(nic,:) = linspace(x6,x4,njc);
yp(nic,:) = linspace(y6,y4,njc);


i=1:nic;
j=1:njc;

fic = (i-1)/(nic-1);
fjc = (j-1)/(njc-1);

[Jc,Ic]=meshgrid(fjc,fic);


x=interp2(Jc,Ic,xp,J,I);
y=interp2(Jc,Ic,yp,J,I);

blk{9}.x = x;
blk{9}.y = y;



clear xc yc

NB = length(blk);
for ib=1:NB
[nib{ib},njb{ib}] = size(blk{ib}.x);    
NI{ib} = nib{ib};
NJ{ib} = njb{ib};
end

% patch interfaces

% Block 1
next_block{1}.im = 0;
next_block{1}.ip = 3;
next_block{1}.jm = 2;
next_block{1}.jp = 2;

next_patch{1}.im = 1; % inlet
next_patch{1}.ip = 3;
next_patch{1}.jm = 4;
next_patch{1}.jp = 3;

% Block 2
next_block{2}.im = 0;
next_block{2}.ip = 6;
next_block{2}.jm = 1;
next_block{2}.jp = 1;

next_patch{2}.im = 1; % inlet
next_patch{2}.ip = 1;
next_patch{2}.jm = 4;
next_patch{2}.jp = 3;

next_block{3}.im = 4;
next_block{3}.ip = 5;
next_block{3}.jm = 1;
next_block{3}.jp = 0;

next_patch{3}.im = 1;
next_patch{3}.ip = 1;
next_patch{3}.jm = 2;
next_patch{3}.jp = 3;  % wall

next_block{4}.im = 3;
next_block{4}.ip = 7;
next_block{4}.jm = 6;
next_block{4}.jp = 0;

next_patch{4}.im = 1;
next_patch{4}.ip = 1;
next_patch{4}.jm = 4;
next_patch{4}.jp = 3;  % wall


next_block{5}.im = 3;
next_block{5}.ip = 7;
next_block{5}.jm = 6;
next_block{5}.jp = 0;

next_patch{5}.im = 2;
next_patch{5}.ip = 2;
next_patch{5}.jm = 3;
next_patch{5}.jp = 3;  % wall


next_block{6}.im = 2;
next_block{6}.ip = 9;
next_block{6}.jm = 5;
next_block{6}.jp = 4;

next_patch{6}.im = 2;
next_patch{6}.ip = 1;
next_patch{6}.jm = 3;
next_patch{6}.jp = 3;

next_block{7}.im = 4;
next_block{7}.ip = 5;
next_block{7}.jm = 8;
next_block{7}.jp = 0;

next_patch{7}.im = 2;
next_patch{7}.ip = 2;
next_patch{7}.jm = 1;
next_patch{7}.jp = 3;  % wall

next_block{8}.im = 7;
next_block{8}.ip = 0; 
next_block{8}.jm = 9;
next_block{8}.jp = 9;

next_patch{8}.im = 3;
next_patch{8}.ip = 2; % exit
next_patch{8}.jm = 4;
next_patch{8}.jp = 3;

next_block{9}.im = 6;
next_block{9}.ip = 0;
next_block{9}.jm = 8;
next_block{9}.jp = 8;

next_patch{9}.im = 2;
next_patch{9}.ip = 2; % exit
next_patch{9}.jm = 4;
next_patch{9}.jp = 3;

% corners
corner{1}.block = {1 2 3 4 6};
corner{1}.i = {nib{1} nib{2} 1 1 1};
corner{1}.j = {1 njb{2} 1 1 njb{6}};

corner{2}.block = {1 2 3 5 6};
corner{2}.i = {nib{1} nib{2} nib{3} 1 1};
corner{2}.j = {njb{1} 1 1 1 1};

corner{3}.block = {8 9 7 4 6};
corner{3}.i = {1 1 1 nib{4} nib{6}};
corner{3}.j = {1 njb{9} 1 1 njb{6}};

corner{4}.block = {8 9 7 5 6};
corner{4}.i = {1 1 nib{7} nib{5} nib{6}};
corner{4}.j = {njb{8} 1 1 1 1};
for n=1:length(corner)
corner{n}.Nb = length(corner{n}.block);    
end

% call mesh smooth to solve for block boundaries

[blk] = mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall);

return
