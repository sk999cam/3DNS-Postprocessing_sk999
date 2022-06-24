clear all

dir = '';

% read in patch interfaces from datum input file.
blockdims = load([dir,'blockdims.txt']);
fidin = fopen([dir,'input_gpu.txt'],'r');
A = fscanf(fidin,'%d %d',2);
NB=A(1);
for ib=1:NB
    
   nijk = fscanf(fidin,'%d %d %d',3);
  
   %nprocs = fscanf(fidin,'%d %d %d', 3);    
      
   patches{ib} = fscanf(fidin,'%d %d %d %d',4);
   
   im = patches{ib}(1);
   ip = patches{ib}(2);
   jm = patches{ib}(3);
   jp = patches{ib}(4);
   
   im_next_block{ib}= 0;
   im_next_patch{ib}= im;
   
   ip_next_block{ib}= 0;
   ip_next_patch{ib}= ip;
   
   jm_next_block{ib}= 0;
   jm_next_patch{ib}= jm;
   
   jp_next_block{ib}= 0;
   jp_next_patch{ib}= jp;
   
   if(im==0)
   A=fscanf(fidin,'%d %d',2);
   im_next_block{ib}=A(1);
   im_next_patch{ib}=A(2);
   end
   
   if(ip==0)
   A=fscanf(fidin,'%d %d',2);
   ip_next_block{ib}=A(1);
   ip_next_patch{ib}=A(2);
   end
   
   if(jm==0)
   A=fscanf(fidin,'%d %d',2);
   jm_next_block{ib}=A(1);
   jm_next_patch{ib}=A(2);
   end
   
   if(jp==0)
   A=fscanf(fidin,'%d %d',2);
   jp_next_block{ib}=A(1);
   jp_next_patch{ib}=A(2);
   end
   
end


% now read in corners information
ncorner=fscanf(fidin,'%d',1);
for n=1:ncorner
A=fscanf(fidin,'%d',1);    
corner{n}.Nb = A(1);
for m=1:A(1)
A=fscanf(fidin,'%d %d %d',3);        
corner{n}.block{m} = A(1);
corner{n}.i{m} = A(2);
corner{n}.j{m} = A(3);
end
end

fclose(fidin);


% now read in datum mesh
for ib=1:NB
    
   fid = fopen([dir,'grid_',num2str(ib),'.txt'],'r');
   ni = blockdims(ib,1);
   nj = blockdims(ib,2);
   nk = blockdims(ib,3);

   x = zeros(ni,nj);
   y = zeros(ni,nj);
   
   for j=1:nj
   for i=1:ni  
   [A]=fscanf(fid, '%f %f',2);
   x(i,j) = A(1);
   y(i,j) = A(2);
   end
   end
   fclose(fid);
   
   Ni(ib) = ni;
   Nj(ib) = nj;
   
blk{ib}.x = x;
blk{ib}.y = y;
   
end

for ib=1:NB
next_block{ib}.im = im_next_block{ib};
next_block{ib}.ip = ip_next_block{ib};
next_block{ib}.jm = jm_next_block{ib};
next_block{ib}.jp = jp_next_block{ib};
next_patch{ib}.im = im_next_patch{ib};
next_patch{ib}.ip = ip_next_patch{ib};
next_patch{ib}.jm = jm_next_patch{ib};
next_patch{ib}.jp = jp_next_patch{ib};
end

% get pitch
ymin = 1e12;
ymax = -1e12;
for ib=1:NB

if patches{ib}(1)==1

miny = min(blk{ib}.y(1,:));
maxy = max(blk{ib}.y(1,:));

if miny < ymin
    ymin = miny;
end

if maxy > ymax
    ymax = maxy;
end

end

end

pitch = ymax - ymin;

% 
% % fix periodics
% % now check periodic interfaces
% for ib=1:NB 
% [blk]=make_periodic(blk,next_block,next_patch,ib,pitch);
% end


ng = 2; % ghost points for block boundary smoothing
% call multiblock to find the block boundary interfaces
for ib=1:NB 
[up{ib},dn{ib}] = multiblock(blk,next_block,next_patch,ib,pitch,ng);
end

% check patches here
for ib=1:NB 

im = patches{ib}(1);
ip = patches{ib}(2);
jm = patches{ib}(3);
jp = patches{ib}(4);
 

if(im == 0)
x = blk{ib}.x(1,:); 
y = blk{ib}.y(1,:);
xp = up{ib}.xi(1,:);
yp = up{ib}.yi(1,:);
dx = x - xp;
dy = y - yp;
d = sqrt(dx.*dx + dy.*dy);
err = sum(d);

if(err>1e-24)
    'im'
    ib
    err
    [maxd,imax]=max(d)
end
plot(x,y,'ro',xp,yp,'kx');
hold on
end

if(ip == 0)
x = blk{ib}.x(end,:); 
y = blk{ib}.y(end,:);
xp = dn{ib}.xi(1,:);
yp = dn{ib}.yi(1,:);
dx = x - xp;
dy = y - yp;
d = sqrt(dx.*dx + dy.*dy);
err = sum(d);

if(err>1e-24)
    'ip'
    ib
    err
    [maxd,imax]=max(d)
end
plot(x,y,'ro',xp,yp,'kx');
hold on
end



if(jm == 0)
x = blk{ib}.x(:,1); 
y = blk{ib}.y(:,1);
xp = up{ib}.xj(:,1);
yp = up{ib}.yj(:,1);
dx = x - xp;
dy = y - yp;
d = sqrt(dx.*dx + dy.*dy);
err = sum(d);

if(err>1e-24)
    'jm'
    ib
    err
    [maxd,imax]=max(d)
end
plot(x,y,'ro',xp,yp,'kx');
hold on
end

if(jp == 0)
x = blk{ib}.x(:,end); 
y = blk{ib}.y(:,end);
xp = dn{ib}.xj(:,1);
yp = dn{ib}.yj(:,1);
dx = x - xp;
dy = y - yp;
d = sqrt(dx.*dx + dy.*dy);
err = sum(d);

if(err>1e-24)
    'jp'
    ib
    err
    [maxd,imax]=max(d)
end
plot(x,y,'ro',xp,yp,'kx');
hold on
end

end

% 
% 
% % write grid
% for ii=1:NB
% x = blk{ii}.x;
% y = blk{ii}.y;
% [ni,nj,nk]=size(x);
%     
% fid = fopen([dir,'grid_',num2str(ii),'.txt'],'w')
% for k=1:nk
% for j=1:nj
% for i=1:ni
% fprintf(fid,'%20.16e %20.16e\n',x(i,j,k),y(i,j,k));
% end
% end
% end
% fclose(fid);
% 
% figure(10)
% pcolor(x(:,:,1),y(:,:,1),x(:,:,1)), hold on
% 
% end
