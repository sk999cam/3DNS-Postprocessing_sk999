function err = write_plot3d_extruded(blk,filename,nk,span)
    NB = length(blk.x);
    blkdims = zeros(2,NB);
    for i=1:NB
        blkdims(:,i) = size(blk.x{i})';
    end
    blkdims(3,:) = nk*ones(1, NB);
    zv(1,1,:) = linspace(0,span,nk);
    f = fopen(filename,'w');
    fprintf(f,'%d\n',NB);
    fprintf(f,'%d %d %d\n',blkdims);
    for i=1:NB
        x = repmat(blk.x{i}, 1, 1, nk);
        y = repmat(blk.y{i}, 1, 1, nk);
        z = repmat(zv, blkdims(1,i), blkdims(2,i), 1);
        fprintf(f,'%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',x);
        fprintf(f,'%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',y);
        fprintf(f,'%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',z);
    end
    err = 0;
    fclose(f);
end