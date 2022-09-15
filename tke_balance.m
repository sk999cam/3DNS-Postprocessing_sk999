function [e_Pr, e_diss] = tke_balance(mF, regions)

    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
    
    %%
    
    celldims = mF.blk.blockdims(:,1:2) - 1;
    cells.x = zeros([sum(prod(celldims,2)) 1]);
    cells.y = zeros([sum(prod(celldims,2)) 1]);
    cells.area = zeros([sum(prod(celldims,2)) 1]);
    cells.prod = zeros([sum(prod(celldims,2)) 1]);
    cells.diss = zeros([sum(prod(celldims,2)) 1]);
    cells.i = zeros([sum(prod(celldims,2)) 1]);
    cells.j = zeros([sum(prod(celldims,2)) 1]);
    cells.ib = zeros([sum(prod(celldims,2)) 1]);
    
    
    %%
    
    parfor ib=1:mF.NB
        nentries = celldims(ib,1)*celldims(ib,2);
    
        Pr_prop = mF.Pr{ib};
        diss_prop = mF.diss{ib};
    
        xtmp{ib} = zeros(nentries, 1);
        ytmp{ib} = zeros(nentries, 1);
        areatmp{ib} = zeros(nentries, 1);
        Prtmp{ib} = zeros(nentries, 1);
        disstmp{ib} = zeros(nentries, 1);
        itmp{ib} = zeros(nentries, 1);
        jtmp{ib} = zeros(nentries, 1);
        ibtmp{ib} = zeros(nentries, 1);
        
        for i=1:celldims(ib,1)
            if mod(i, 100) == 0
                 fprintf('Block %d, i=%d/%d\n', ib, i, celldims(ib,1))
            end
            for j=1:celldims(ib,2)
                pos = celldims(ib,2)*(i-1)+j;
                xnow = [mF.blk.x{ib}(i,j) mF.blk.x{ib}(i+1,j) ...
                    mF.blk.x{ib}(i+1,j+1) mF.blk.x{ib}(i,j+1)];
                ynow = [mF.blk.y{ib}(i,j) mF.blk.y{ib}(i+1,j) ...
                    mF.blk.y{ib}(i+1,j+1) mF.blk.y{ib}(i,j+1)];
    
                xtmp{ib}(pos) = mean(xnow);
                ytmp{ib}(pos) = mean(ynow);
    
                area = polyarea(xnow,ynow);
                areatmp{ib}(pos) = abs(area);
    
                disstmp{ib}(pos) = 0.25*(diss_prop(i,j)+diss_prop(i+1,j)+diss_prop(i+1,j+1)+diss_prop(i,j+1))*area;
                Prtmp{ib}(pos) = 0.25*(Pr_prop(i,j)+Pr_prop(i+1,j)+Pr_prop(i+1,j+1)+Pr_prop(i,j+1))*area;
                itmp{ib}(pos) = i;
                jtmp{ib}(pos) = j;
                ibtmp{ib}(pos) = ib;
            end
        end
    end
    
    for ib = 1:mF.NB
    
        offset = sum(prod(celldims(1:ib,:),2));
        nentries = celldims(ib,1)*celldims(ib,2);
    
        cells.x(offset-nentries+1:offset) = xtmp{ib};
        cells.y(offset-nentries+1:offset) = ytmp{ib};
        cells.area(offset-nentries+1:offset) = areatmp{ib};
        cells.diss(offset-nentries+1:offset) = disstmp{ib};
        cells.Pr(offset-nentries+1:offset) = Prtmp{ib};
        cells.i(offset-nentries+1:offset) = itmp{ib};
        cells.j(offset-nentries+1:offset) = jtmp{ib};
        cells.ib(offset-nentries+1:offset) = ibtmp{ib};
    end
    
    
    for ir = 1:length(regions)
        inds = cells.ib == regions{ir}.nb & ...
            cells.i >= regions{ir}.is & cells.i <= regions{ir}.ie & ...
            cells.j >= regions{ir}.js & cells.j <= regions{ir}.je;
%         if ir == 1
%             inds = cells.ib == 6 & cells.i>= 100 & cells.i <= 400;
%         end
        e_Pr(ir) = sum(cells.Pr(inds));
        e_diss(ir) = sum(cells.diss(inds));
    end
end
