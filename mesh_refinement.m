function [blk] = mesh_refinement(blk,refine_fac,npp)

NB = length(blk);

for ib=1:NB
    
x = blk{ib}.x;
y = blk{ib}.y;

[ni,nj]=size(x);

ni_new = npp*ceil(refine_fac*ni/npp);
nj_new = npp*ceil(refine_fac*nj/npp);          
        
I = linspace(0,1,ni);
J = linspace(0,1,nj);

[X,Y] = meshgrid(J,I);

II = linspace(0,1,ni_new);
JI = linspace(0,1,nj_new);

[XI,YI] = meshgrid(JI,II);

%xnew = interp2(X,Y,x,XI,YI,'spline');
%ynew = interp2(X,Y,y,XI,YI,'spline');
xnew = interp2(X,Y,x,XI,YI,'cubic');
ynew = interp2(X,Y,y,XI,YI,'cubic');

blk{ib}.x = xnew;
blk{ib}.y = ynew;

end

return