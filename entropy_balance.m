function [e_unst, e_s, e_phi, e_irrev, e_rev] = entropy_balance(mF, regions)

pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

%%

celldims = mF.blk.blockdims(:,1:2) - 1;
cells.x = zeros([sum(prod(celldims,2)) 1]);
cells.y = zeros([sum(prod(celldims,2)) 1]);
cells.area = zeros([sum(prod(celldims,2)) 1]);
cells.unst = zeros([sum(prod(celldims,2)) 1]);
cells.conv = zeros([sum(prod(celldims,2)) 1]);
cells.diss = zeros([sum(prod(celldims,2)) 1]);
cells.irrev = zeros([sum(prod(celldims,2)) 1]);
cells.rev = zeros([sum(prod(celldims,2)) 1]);
cells.i = zeros([sum(prod(celldims,2)) 1]);
cells.j = zeros([sum(prod(celldims,2)) 1]);
cells.ib = zeros([sum(prod(celldims,2)) 1]);


%%

parfor ib=1:mF.NB
    nentries = celldims(ib,1)*celldims(ib,2);

   % dsdt = (s2{ib} - s1{ib})/hfcase.meanFlow.meanTime;

    [dsdx,dsdy] = gradHO(mF.blk.x{ib},mF.blk.y{ib},mF.s{ib});
    [dqx_Tdx,~] = gradHO(mF.blk.x{ib},mF.blk.y{ib},mF.rev_gen_x{ib});
    [~,dqy_Tdy] = gradHO(mF.blk.x{ib},mF.blk.y{ib},mF.rev_gen_y{ib});
    [drousdx, ~] = gradHO(mF.blk.x{ib},mF.blk.y{ib},mF.rous{ib});
    [~, drovsdy] = gradHO(mF.blk.x{ib},mF.blk.y{ib},mF.rovs{ib});

    unst_prop = mF.ro{ib}.*mF.s{ib};

%     conv_prop = mF.ro{ib}.*(mF.u{ib}.*dsdx + ...
%         mF.v{ib}.*dsdy);
    
    conv_prop = drousdx + drovsdy;

    diss_prop = mF.diss_T{ib};
    irrev_prop = mF.irrev_gen{ib};
    rev_prop = dqx_Tdx + dqy_Tdy;

    xtmp{ib} = zeros(nentries, 1);
    ytmp{ib} = zeros(nentries, 1);
    areatmp{ib} = zeros(nentries, 1);
    unsttmp{ib} = zeros(nentries,1);
    convtmp{ib} = zeros(nentries, 1);
    disstmp{ib} = zeros(nentries, 1);
    irrevtmp{ib} = zeros(nentries, 1);
    revtmp{ib} = zeros(nentries, 1);
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

            unsttmp{ib}(pos) = 0.25*(unst_prop(i,j)+unst_prop(i+1,j)+unst_prop(i+1,j+1)+unst_prop(i,j+1))*area;
            convtmp{ib}(pos) = 0.25*(conv_prop(i,j)+conv_prop(i+1,j)+conv_prop(i+1,j+1)+conv_prop(i,j+1))*area;
            disstmp{ib}(pos) = 0.25*(diss_prop(i,j)+diss_prop(i+1,j)+diss_prop(i+1,j+1)+diss_prop(i,j+1))*area;
            irrevtmp{ib}(pos) = 0.25*(irrev_prop(i,j)+irrev_prop(i+1,j)+irrev_prop(i+1,j+1)+irrev_prop(i,j+1))*area;
            revtmp{ib}(pos) = 0.25*(rev_prop(i,j)+rev_prop(i+1,j)+rev_prop(i+1,j+1)+rev_prop(i,j+1))*area;
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
    cells.unst(offset-nentries+1:offset) = unsttmp{ib};
    cells.conv(offset-nentries+1:offset) = convtmp{ib};
    cells.diss(offset-nentries+1:offset) = disstmp{ib};
    cells.irrev(offset-nentries+1:offset) = irrevtmp{ib};
    cells.rev(offset-nentries+1:offset) = revtmp{ib};
    cells.i(offset-nentries+1:offset) = itmp{ib};
    cells.j(offset-nentries+1:offset) = jtmp{ib};
    cells.ib(offset-nentries+1:offset) = ibtmp{ib};
end


for ir = 1:length(regions)
    inds = cells.ib == regions{ir}.nb & ...
        cells.i >= regions{ir}.is & cells.i <= regions{ir}.ie & ...
        cells.j >= regions{ir}.js & cells.j <= regions{ir}.je;
%     if ir == 1
%         inds = cells.ib == 6 & cells.i>= 100 & cells.i <= 400;
%     end
    e_unst(ir) = sum(cells.unst(inds));
    e_s(ir) = sum(cells.conv(inds));
    e_phi(ir) = sum(cells.diss(inds));
    e_irrev(ir) = sum(cells.irrev(inds));
    e_rev(ir) = sum(cells.rev(inds));
end
end
