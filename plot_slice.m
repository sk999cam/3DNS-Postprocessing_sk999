function plot_slice(slice,blk,prop,ax)
    NB = length(slice);
    hold on
    for i=1:NB
        i
        q = slice(1,i).(prop);
        size(blk{i}.x)
        size(blk{i}.y)
        size(q)
        pcolor(ax,blk{i}.x,blk{i}.y,q);
    end
    shading('interp')
    axis equal
end