function plot_mesh(blk, skip)
if nargin < 2
    skip=1;
end
NB = length(blk.x);
figure;
hold on
for i=1:NB
    
    xnew = blk.x{i};
    ynew = blk.y{i};
    [ni,nj] = size(xnew);
    is = 1:skip:ni;
    is = [is(1:end-1) ni];
    js = 1:skip:nj;
    js = [js(1:end-1) nj];
    xnew = xnew(is, js);
    ynew = ynew(is,js);
    plot(xnew,ynew,'k');
    plot(xnew',ynew','k');
    plot(blk.x{i}(1,:),blk.y{i}(1,:),'r');
    plot(blk.x{i}(end,:),blk.y{i}(end,:),'r');
    plot(blk.x{i}(:,1),blk.y{i}(:,1),'r');
    plot(blk.x{i}(:,end),blk.y{i}(:,end),'r');

end
axis equal
end