function nodes = mesh2nodes(blk,points)
    nblocks = length(blk.x);
    for nb=1:nblocks
        node_xs = points.x_coordinate;
        node_ys = points.y_coordinate;
        nodenums = points.nodenumber;
        nib = size(blk.x{nb},1);
        njb = size(blk.x{nb},2);
        nodestmp = zeros(size(blk.x{nb}));
        xtmp = blk.x{nb};
        ytmp = blk.y{nb};
        for i=1:nib
            fprintf('Block: %d, i= %d / %d\n',nb,i,nib)
            parfor j=1:njb
                dx = node_xs - xtmp(i,j);
                dy = node_ys - ytmp(i,j);

                [~,ind] = min(dx.^2+dy.^2);
                nodestmp(i,j) = nodenums(ind);
            end
        end
        nodes{nb} = nodestmp;
    end
end