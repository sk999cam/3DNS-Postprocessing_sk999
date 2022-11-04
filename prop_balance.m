function [vol_int, surf_int] = prop_balance(prop,blk,cell_areas,regions, plot)

    if nargin < 5 || isempty(plot)
        plot = false;
    end

    NB = length(blk.x);

    celldims = blk.blockdims(:,1:2) - 1;
%     cells.x = zeros([sum(prod(celldims,2)) 1]);
%     cells.y = zeros([sum(prod(celldims,2)) 1]);
    cells.area = zeros([sum(prod(celldims,2)) 1]);
    cells.prop = zeros([sum(prod(celldims,2)) 1]);
    cells.i = zeros([sum(prod(celldims,2)) 1]);
    cells.j = zeros([sum(prod(celldims,2)) 1]);
    cells.ib = zeros([sum(prod(celldims,2)) 1]);


    
    for ir = 1:length(regions)
        nb = regions{ir}.nb;
        irange = regions{ir}.is:regions{ir}.ie;
        jrange = regions{ir}.js:regions{ir}.je;
        xline = [blk.x{nb}(irange,jrange(1))' blk.x{nb}(irange(end),jrange(2:end)) ...
            blk.x{nb}(irange(end-1:-1:1),jrange(end))' blk.x{nb}(irange(1), jrange(end-1:-1:1))];
        yline = [blk.y{nb}(irange,jrange(1))' blk.y{nb}(irange(end),jrange(2:end)) ...
            blk.y{nb}(irange(end-1:-1:1),jrange(end))' blk.y{nb}(irange(1), jrange(end-1:-1:1))];
        surf_area(ir) = polyarea(xline,yline);


        if size(prop{1},3) > 1
            dA1 = yline(2:end) - yline(1:end-1);
            dA2 = xline(1:end-1) - xline(2:end);

            if ispolycw(xline, yline)
                dA1 = -dA1;
                dA2 = -dA2;
            end

            u1 = [prop{nb}(irange,jrange(1),1)' prop{nb}(irange(end),jrange(2:end),1) ...
            prop{nb}(irange(end-1:-1:1),jrange(end),1)' prop{nb}(irange(1), jrange(end-1:-1:1),1)];

            u2 = [prop{nb}(irange,jrange(1),2)' prop{nb}(irange(end),jrange(2:end),2) ...
            prop{nb}(irange(end-1:-1:1),jrange(end),2)' prop{nb}(irange(1), jrange(end-1:-1:1),2)];

            dF = 0.5*( dA1.*(u1(1:end-1) + u1(2:end)) + dA2.*(u2(1:end-1) + u2(2:end)) );
            surf_int(ir) = sum(dF);
        else
            surf_int(ir) = 0;
        end
    end

%     for ib = 1:NB
% 
%         x11 = blk.x{ib}(1:end-1,1:end-1);
%         x12 = blk.x{ib}(1:end-1,2:end);
%         x21 = blk.x{ib}(2:end,1:end-1);
%         x22 = blk.x{ib}(2:end,2:end);
% 
%         y11 = blk.y{ib}(1:end-1,1:end-1);
%         y12 = blk.y{ib}(1:end-1,2:end);
%         y21 = blk.y{ib}(2:end,1:end-1);
%         y22 = blk.y{ib}(2:end,2:end);
% 
%         xmid{ib} = 0.25*(x11+x12+x21+x22);
%         ymid{ib} = 0.25*(y11+y12+y21+y22);
    

    parfor ib=1:NB
        nentries = celldims(ib,1)*celldims(ib,2);
    
        areatmp{ib} = zeros(nentries, 1);
        proptmp{ib} = zeros(nentries, 1);
        itmp{ib} = zeros(nentries, 1);
        jtmp{ib} = zeros(nentries, 1);
        ibtmp{ib} = ib*ones(nentries, 1);

        if size(prop{1},3) > 1
            [dudx, ~] = gradHO(blk.x{ib}, blk.y{ib}, prop{ib}(:,:,1));
            [~, dvdy] = gradHO(blk.x{ib}, blk.y{ib}, prop{ib}(:,:,2));
            propnow = dudx + dvdy;
        else
            propnow = prop{ib};
        end
        
        for i=1:celldims(ib,1)
            if mod(i, 100) == 0
%                  fprintf('Block %d, i=%d/%d\n', ib, i, celldims(ib,1))
            end
            for j=1:celldims(ib,2)
                pos = celldims(ib,2)*(i-1)+j;
                dA = cell_areas{ib}(i,j);
                areatmp{ib}(pos) = dA
                itmp{ib}(pos) = i;
                jtmp{ib}(pos) = j;

                proptmp{ib}(pos) = 0.25*(propnow(i,j)+propnow(i+1,j)+propnow(i+1,j+1)+propnow(i,j+1))*dA;

            end
        end
    end

    for ib = 1:NB
    
        offset = sum(prod(celldims(1:ib,:),2));
        nentries = celldims(ib,1)*celldims(ib,2);
   
        cells.area(offset-nentries+1:offset) = areatmp{ib};
        cells.prop(offset-nentries+1:offset) = proptmp{ib};
        cells.i(offset-nentries+1:offset) = itmp{ib};
        cells.j(offset-nentries+1:offset) = jtmp{ib};
        cells.ib(offset-nentries+1:offset) = ibtmp{ib};

        blk.region{ib} = zeros(size(blk.x{ib}));
    end

    for ir = 1:length(regions)
        inds = cells.ib == regions{ir}.nb & ...
            cells.i >= regions{ir}.is & cells.i < regions{ir}.ie & ...
            cells.j >= regions{ir}.js & cells.j < regions{ir}.je;
        vol_int(ir) = sum(cells.prop(inds));
        vol_area(ir) = sum(cells.area(inds));
        cells.region(inds) = ir/(length(regions));
        if surf_int(ir) == 0
            fprintf('Region %d: Total volume = %s, \x222b dV  = %s, \x222b\x03d5 dV = %s\n', ir, surf_area(ir), vol_area(ir), vol_int(ir));
        else
            fprintf('Region %d: Total volume = %s, \x222b dV  = %s, \x222bu.dA = %s, \x222b\x2207.u dV = %s\n', ir, surf_area(ir), vol_area(ir), surf_int(ir), vol_int(ir));
        end
        points = find(inds);
        for ip = 1:length(points)
            ind = points(ip);
            blk.region{cells.ib(ind)}(cells.i(ind),cells.j(ind)) = ir/length(regions);
        end
    end

    if plot
        ax = gca;
        hold on
        for ib=1:NB
            pcolor(ax,blk.x{ib},blk.y{ib},blk.region{ib})
        end
        shading flat
    end

end