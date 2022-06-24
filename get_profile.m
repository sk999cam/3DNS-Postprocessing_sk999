clear all

%dir = 'original_grid/';
dir = 'harrison/'
blockdims = load([dir,'blockdims.txt']);

Nb = size(blockdims,1);

minx=1e12;
miny=1e12;

maxx=-1e12;
maxy=-1e12;

for nb=1:Nb
    
ni = blockdims(nb,1);
nj = blockdims(nb,2);
nk = blockdims(nb,3);

   fid = fopen([dir,'grid_',num2str(nb),'.txt'],'r');      
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
   
   xb{nb} =x;
   yb{nb} =y;
   
end


pitch = yb{2}(1,end) - yb{1}(1,1);


xprof = [];
yprof = [];

for nb = [3 5 7 4]
xnow = xb{nb}(:,end);
ynow = yb{nb}(:,end);

if(nb==3 || nb==5)
xprof = [xprof; xnow];
yprof = [yprof; ynow];
else
xprof = [xprof; xnow(end:-1:1)];
yprof = [yprof; ynow(end:-1:1)];    
end

end


% get surface distance
sprof = curve_length(xprof,yprof);
slen = sprof(end);
sprof = sprof/slen;
[snow,I,J] = unique(sprof);

% interpolate to remove unique points and distribute uniformly
si = linspace(0,1,5000);
xi = interp1(sprof(I),xprof(I),si,'spline');
yi = interp1(sprof(I),yprof(I),si,'spline');


% write profile
N = length(xi);
N = N-1; % remove last point to avoid repeat

% normalize
[xLE,iLE] = min(xi);
yLE = yi(iLE);
xi = (xi-xLE)/(slen*0.5);
yi = (yi-yLE)/(slen*0.5);
pitch = pitch/(slen*0.5);

fid = fopen(['harrison.txt'],'w');
fprintf(fid,'%d\n',[N]);
fprintf(fid,'%20.16e\n',[pitch]);
for ii=1:N
fprintf(fid,'%20.16e %20.16e\n',xi(ii),yi(ii));
end
fclose(fid);




