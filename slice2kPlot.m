function slice2kPlot(slice, prop, fpath, lims, label)
    h = figure('Visible','off');
    ax = axes(h);
  
    q = slice.(prop);

    hold on
    for i=1:slice.NB
        pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
    end
    shading('interp')
    axis([-0.6 2 -0.5 0.5])
    axis equal
    cb = colorbar(ax);
    if nargin > 3 && ~isempty(lims)
        caxis(lims);
    end
    if nargin > 4 && ~isempty(label)
        cb.Label.String = label;
    end
    set(ax,'FontSize',16);
    exportgraphics(h, fpath, 'Resolution', 600);
end