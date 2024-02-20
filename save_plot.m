function save_plot(dir, fname, save, fmt)

    if nargin < 3 || isempty(save)
        save = true;
    end
    if nargin < 4
        fmt = 'png';
    end
    if save
        f = gcf;
        savefig(f, fullfile(dir, [fname '.fig']))
        switch fmt
            case 'png'
                exportgraphics(f, fullfile(dir, [fname '.' fmt]), 'Resolution',450)
            case 'eps'
                exportgraphics(f, fullfile(dir, [fname '.' fmt]))
            case 'svg'
                set(f, 'Color', 'None');
                figure(f)
                set(gca, 'Color', 'None');
                plot2svg(fullfile(dir, [fname '.' fmt]), f, 'png')
        end
    end

end