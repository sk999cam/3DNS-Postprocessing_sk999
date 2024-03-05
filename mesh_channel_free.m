function [blk] = mesh_channel_free(y_in, y_out, ni, nj, nk, wall_aspect)

    if nargin < 5 || isempty(nk)
        nk = 1;
    end

    if nargin < 6 || isempty(wall_aspect)
        wall_aspect = 15;
    end
    
    nib = length(ni);
    njb = length(nj);
    nis = ni;
    njs = nj;
    ni = sum(nis) - (nib-1);
    nj = sum(njs) - (njb-1);
    
    
    fi = linspace(0, 1, ni);
    L = 1;
    
    ywall = (L/(ni-1))/wall_aspect;
    msmooths = 200;
    fx_expan = 1;
    
    for i=1:ni
        xnow = fi(i)*L;
        yprof = y_in + fi(i)*(y_out - y_in);
        fy_expan = fexpan(yprof/ywall,nj);
        fyexs(i) = fy_expan;
        fj = spacing(nj,fy_expan,0);
        for j=1:nj
            x(i,j) = xnow;
            y(i,j) = yprof*fj(j);
        end
    end
    
    if nib == 1 && njb == 1
        blk.x{1} = x;
        blk.y{1} = y;
        blk.blockdims = [ni, nj, nk];
        nb = 1;
    else
        ibounds = [0 cumsum(nis-1)]+1;
        jbounds = [0 cumsum(njs-1)]+1;
        nb = 0;
        for jb = 1:njb
            for ib = 1:nib
                nb = nb+1;
                blk.x{nb} = x(ibounds(ib):ibounds(ib+1),jbounds(jb):jbounds(jb+1));
                blk.y{nb} = y(ibounds(ib):ibounds(ib+1),jbounds(jb):jbounds(jb+1));
                blk.blockdims(nb,:) = [size(blk.x{nb}) nk];
            end
        end
    end
    
    
    blk.nk = nk;
    
    block_grid = zeros(nib+2,njb+2);
    block_grid(2:nib+1, 2:njb+1) = reshape(1:nb,nib,njb);
    
    nb = 0;
    for jb = 1:njb
        for ib = 1:nib
            nb = nb+1;
            next_block{nb}.im = block_grid(ib, jb+1);
            next_block{nb}.ip = block_grid(ib+2, jb+1);
            next_block{nb}.jm = block_grid(ib+1, jb);
            next_block{nb}.jp = block_grid(ib+1, jb+2);
            if ib == 1
                next_patch{nb}.im = 1;
            else
                next_patch{nb}.im = 2;
            end
            if ib == nib
                next_patch{nb}.ip = 2;
            else
                next_patch{nb}.ip = 1;
            end
            if jb == 1
                next_patch{nb}.jm = 3;
            else
                next_patch{nb}.jm = 4;
            end
            if jb == njb
                next_patch{nb}.jp = 6;
            else
                next_patch{nb}.jp = 3;
            end
        end
    end
    
    if nib > 1 && njb > 1
        nc = 0;
        for jb = 1:njb-1
            for ib = 1:nib-1
                nc = nc+1;
                corner{nc}.Nb = 4;
                b1 = block_grid(ib+1, jb+1);
                b2 = block_grid(ib+2, jb+1);
                b3 = block_grid(ib+1, jb+2);
                b4 = block_grid(ib+2, jb+2);
                corner{nc}.block = {b1 b2 b3 b4};
                corner{nc}.i = {nis(ib) 1 nis(ib) 1};
                corner{nc}.j = {njs(ib) njs(ib) 1 1};
            end
        end
    else
        corner = [];
    end

    blk.nbg = nib*njb;
    blk.block_groups = {};
    for ibg = 1:blk.nbg
        blk.block_groups{ibg} = [ibg];
    end
    
    blk.next_block = next_block;
    blk.next_patch = next_patch;
    blk.corner = corner;

end