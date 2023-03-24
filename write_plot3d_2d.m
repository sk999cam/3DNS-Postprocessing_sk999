function err = write_plot3d_2d(blk,filename)
    NB = length(blk.x);
    blkdims = zeros(2,NB);
    for i=1:NB
        blkdims(:,i) = size(blk.x{i})';
    end
    f = fopen(filename,'w');
    fprintf(f,'%d\n',NB);
    fprintf(f,'%d %d\n',blkdims);
    for i=1:NB
        fprintf(f,'%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',blk.x{i});
        fprintf(f,'%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',blk.y{i});
    end
    err = 0;
    fclose(f);
end