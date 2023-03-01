function slice2kPlot(slice, blk, prop, fpath, lims, label, area, aspect)
    h = figure('Visible','off');
    ax = axes(h);
  
    q = slice.(prop);

    hold on
    for i=1:slice.NB
        pcolor(ax, blk.x{i}, blk.y{i}, q{i});
    end
    shading('interp')
    axis equal
    pbaspect(aspect)
    axis(area)
    axis off
    cb = colorbar(ax);
    if nargin > 4 && ~isempty(lims)
        caxis(lims);
    end
    if nargin > 5 && ~isempty(label)
        cb.Label.String = label;
    end
    set(ax,'FontSize',16);
    exportgraphics(h, fpath, 'Resolution', 600);
end
