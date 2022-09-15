function blk = thin_blocks(blk, ibs, dirs, factor, npp)

% Dir = 1: Expand from i=1 edge
% Dir = 2: Expand from i=ni edge
% Dir = 3: Expand from j=1 edge
% Dir = 4: Expand from j=nj edge

    for ii = 1:length(ibs)
        sprintf('Thinning block %d/%d\n', ii, length(ibs));
        ib = ibs(ii);
        dir = dirs(ii);
    
        x_old = blk{ib}.x;
        y_old = blk{ib}.y;
        [ni, nj] = size(x_old);
    
        if dir == 1
            ni_new = npp*ceil(ni/npp/factor);
            fi = get_new_spacing(ni, ni_new);
            fj = linspace(0,1,nj);

        elseif dir == 2
            ni_new = npp*ceil(ni/npp/factor);
            fi = get_new_spacing(ni, ni_new, 1);
            fj = linspace(0,1,nj);
            fi;

        elseif dir == 3
            nj_new = npp*ceil(nj/npp/factor);
            fj = get_new_spacing(nj, nj_new);
            fi = linspace(0,1,ni);

        elseif dir == 4
            nj_new = npp*ceil(nj/npp/factor);
            fj = get_new_spacing(nj, nj_new, 1);
            fi = linspace(0,1,ni);
        end

        
        
        fic=linspace(0,1,ni);
        fjc=linspace(0,1,nj);
    
        [J,I]=meshgrid(fj,fi);
        [Jc,Ic]=meshgrid(fjc,fic);
    
        blk{ib}.x=interp2(Jc,Ic,x_old,J,I);
        blk{ib}.y=interp2(Jc,Ic,y_old,J,I);
    end
end

function f = get_new_spacing(n1, n2, iflip)

    if nargin < 3
        iflip = 0;
    end

    r0 = 1;
    r = 1.1;
    r1 = 1.2;

    while r - r0 > 1e-6
        ds = 1.0/(n1-1);
        dist = cumsum(ds*r.^(0:n2-2));
        dist = dist(end);
        if dist > 1.0
            r1 = r;
        else
            r0 = r;
            
        end
        r = (r0+r1)/2;
    end

    if iflip
        f = [0.0 cumsum(flip(ds*r.^(0:n2-2)))];
        f = f/f(end);
        
        
    else
        f = [0.0 cumsum(ds*r.^(0:n2-2))];
        f = f/f(end);
    end

end