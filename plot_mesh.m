function plot_mesh(blk)

NB = length(blk)
figure;
hold on
for i=1:NB
    xnew = blk{i}.x;
    ynew = blk{i}.y;
    plot(xnew,ynew,'k');
    plot(xnew',ynew','k');
end