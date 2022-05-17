function slice2jPlot(slice, prop, fpath, lims, label)
    h = figure('Visible','off');
    ax = axes(h);

    q = slice.(prop);
    dx = abs(slice.X(end,1)-obj.X(1,1));
    dy = abs(slice.Z(1,end)-obj.Z(1,1));
    
    pcolor(ax, slice.X, slice.Z, q);
    xlabel('Surface distance')
    shading('interp')
    axis equal
    
    cb = colorbar('southoutside');
    if nargin > 3 && ~isempty(lims)
        caxis(lims)
    end
    if nargin > 4 && ~isempty(label)
        cb.Label.String = label;
    end
    pbaspect(ax, [dx dy dx]);

    if nargin > 3 && ~isempty(lims)
        caxis(lims);
    end
    if nargin > 4 && ~isempty(label)
        cb.Label.String = label;
    end

    set(ax,'FontSize',16);
    exportgraphics(h, fpath, 'Resolution', 600);
end