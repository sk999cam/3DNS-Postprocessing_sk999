function nodes = mesh2nodes(blk,points)
    nblocks = length(blk.x);
    node_xs = points.x_coordinate;
    [node_xs, I] = sort(node_xs);

    node_ys = points.y_coordinate;
    node_ys = node_ys(I);
    nodenums = points.nodenumber;
    nodenums = nodenums(I);

    for nb=1:nblocks
        
        nib = size(blk.x{nb},1);
        njb = size(blk.x{nb},2);
        nodestmp = zeros(size(blk.x{nb}));
        xtmp = blk.x{nb};
        ytmp = blk.y{nb};
        indsnow = node_xs <= max(xtmp,[],'all') & node_xs >= min(xtmp,[],'all') ...
            & node_ys <= max(ytmp,[],'all') & node_ys >= min(ytmp,[],'all');
        nodx = node_xs(indsnow);
        nody = node_ys(indsnow);
        nodnum = nodenums(indsnow);

        for i=1:nib
            fprintf('Block: %d, i= %d / %d\n',nb,i,nib)
            parfor j=1:njb

                dx = nodx - xtmp(i,j);
                shortlist = abs(dx) < 0.0002;
                %dx = dx(shortlist);
                dy = nody - ytmp(i,j);
                dy = dy(shortlist);
                nodesnow = nodenums(shortlist);

                [~,ind] = min(dx.^2+dy.^2);
                nodestmp(i,j) = nodnum(ind);
            end
        end
        nodes{nb} = nodestmp;
    end
end