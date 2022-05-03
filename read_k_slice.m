function slice = read_k_slice(casedir, nSlice, blockdims, gas)
    gam = gas.gam;
    cp = gas.cp;
    cv = cp/gam;
    rgas = cp-cv;
    temp_slash = '/'; if ispc, temp_slash = '\'; end
    slice = struct([]);
    for nb = 1:size(blockdims, 1)
        flofile = fopen([casedir temp_slash 'kcu2_' num2str(nb) '_' num2str(nSlice)],'r');
        nodfile = fopen([casedir temp_slash 'knd2_' num2str(nb) '_' num2str(nSlice)],'r');
        A = fread(flofile,inf,'float64');
        A = reshape(A,5,length(A)/5);
        
        B = fread(nodfile,inf,'uint32');
        B = reshape(B,3,length(B)/3);

        fclose(flofile);
        fclose(nodfile);

        ro = zeros(blockdims(nb,1),blockdims(nb,2));
        ru = zeros(blockdims(nb,1),blockdims(nb,2));
        rv = zeros(blockdims(nb,1),blockdims(nb,2));
        rw = zeros(blockdims(nb,1),blockdims(nb,2));
        Et = zeros(blockdims(nb,1),blockdims(nb,2));

        slice(nb).ro = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).u = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).v = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).w = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).Et = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).p = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).T = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).M = zeros(blockdims(nb,1),blockdims(nb,2)); 
        slice(nb).s = zeros(blockdims(nb,1),blockdims(nb,2));
        slice(nb).vel = zeros(blockdims(nb,1),blockdims(nb,2));
        nb
        size(slice(nb).ro)

        for n=1:size(A,2)
            i = B(1,n);
            j = B(2,n);
            k = B(3,n);
            ro(i,j) = A(1,n);
            ru(i,j) = A(2,n);
            rv(i,j) = A(3,n);
            rw(i,j) = A(4,n);
            Et(i,j) = A(5,n);
        end
        size(A)
        size(B)
        

        size(ro)
        blockdims(nb,:)


        slice(nb).ro = ro;
        u = ru./ro;
        v = rv./ro;
        w = rw./ro;
        slice(nb).u = u;
        slice(nb).v = v;
        slice(nb).w = w;
        slice(nb).Et = Et;
        p = (gam-1)*(Et - 0.5*(u.*u + v.*v + w.*w).*ro);
        slice(nb).p = p;
        T = p./(ro*rgas);
        slice(nb).T = T;
        
        slice(nb).s = cp*log(T/300) - rgas*log(p/1e5);
        slice(nb).vel = sqrt(u.*u + v.*v + w.*w);
        slice(nb).M = slice(nb).vel./sqrt(gam*rgas*T);

    end
end