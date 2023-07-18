function blk = thin_blocks(blk, ibs, dirs, factor, npp)

if length(factor) == 1
    factor(1:length(ibs)) = factor;
end

% Dir = 1: Expand from i=1 edge
% Dir = 2: Expand from i=ni edge
% Dir = 3: Expand from j=1 edge
% Dir = 4: Expand from j=nj edge

    for ii = 1:length(ibs)
        fprintf('Thinning block %d/%d\n', ii, length(ibs));
        ib = ibs(ii);
        dir = dirs(ii);
    
        x_old = blk.x{ib};
        y_old = blk.y{ib};
        [ni, nj] = size(x_old);
    
        if dir == 1
            ni_new = npp*ceil(ni/npp/factor(ii));
            fi = get_new_spacing(ni, ni_new);
            fj = linspace(0,1,nj);

        elseif dir == 2
            ni_new = npp*ceil(ni/npp/factor(ii));
            fi = get_new_spacing(ni, ni_new, 1);
            fj = linspace(0,1,nj);
            fi;

        elseif dir == 3
            nj_new = npp*ceil(nj/npp/factor(ii));
            fj = get_new_spacing(nj, nj_new);
            fi = linspace(0,1,ni);

        elseif dir == 4
            nj_new = npp*ceil(nj/npp/factor(ii));
            fj = get_new_spacing(nj, nj_new, 1);
            fi = linspace(0,1,ni);

        elseif dir == 5
            ni_new = npp*ciel(ni/npp/factor(ii));
            fi = thin_middle(ni, ni_new);
            fj = linspace(0, 1, nj);
        elseif dir == 6
            nj_new = npp*ceil(nj/npp/factor(ii));
            fj = thin_middle(nj, nj_new);
            fi = linspace(0, 1, ni);
        end

        
        
        fic=linspace(0,1,ni);
        fjc=linspace(0,1,nj);
    
        [J,I]=meshgrid(fj,fi);
        [Jc,Ic]=meshgrid(fjc,fic);
    
        blk.x{ib}=interp2(Jc,Ic,x_old,J,I);
        blk.y{ib}=interp2(Jc,Ic,y_old,J,I);
    end
end

function f = get_new_spacing(n1, n2, iflip)

    if nargin < 3
        iflip = 0;
    end

    r0 = 1;
    r = 1.1;
    r1 = 1.2;

    while r - r0 > 1e-10
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

function f = thin_middle(n1, n2)

    ds = 1/(n1-1);
    r0 = 1;
    r = 1.1;
    r2 = 1.2;
    n = floor(n2/2);

    while r - r0 > 1e-10
        dist = cumsum(ds*r.^([0:(n-1) (n2-n-2):-1:0]));
        tot = dist(end);
        if tot > 1.0
            r1 = r;
        else
            r0 = r;
            
        end
        r = (r0+r1)/2;
    end

    f = [0 dist]/dist(end);
end