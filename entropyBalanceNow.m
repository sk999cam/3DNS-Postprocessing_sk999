hfcase = tripcase;

pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

%%

celldims = hfcase.blk.blockdims(:,1:2) - 1;
cells.x = zeros([sum(prod(celldims,2)) 1]);
cells.y = zeros([sum(prod(celldims,2)) 1]);
cells.area = zeros([sum(prod(celldims,2)) 1]);
%cells.unst = zeros([sum(prod(celldims,2)) 1]);
cells.conv = zeros([sum(prod(celldims,2)) 1]);
cells.diss = zeros([sum(prod(celldims,2)) 1]);
cells.irrev = zeros([sum(prod(celldims,2)) 1]);
cells.rev = zeros([sum(prod(celldims,2)) 1]);
cells.i = zeros([sum(prod(celldims,2)) 1]);
cells.j = zeros([sum(prod(celldims,2)) 1]);
cells.ib = zeros([sum(prod(celldims,2)) 1]);


%%

parfor ib=1:hfcase.NB
    nentries = celldims(ib,1)*celldims(ib,2);

   % dsdt = (s2{ib} - s1{ib})/hfcase.meanFlow.meanTime;

    [dsdx,dsdy] = gradHO(hfcase.blk.x{ib},hfcase.blk.y{ib},hfcase.meanFlow.s{ib});
    [dqx_Tdx,~] = gradHO(hfcase.blk.x{ib},hfcase.blk.y{ib},hfcase.meanFlow.rev_gen_x{ib});
    [~,dqy_Tdy] = gradHO(hfcase.blk.x{ib},hfcase.blk.y{ib},hfcase.meanFlow.rev_gen_y{ib});

    %unst_prop = hfcase.meanFlow.ro{ib}.*dsdt;

    conv_prop = hfcase.meanFlow.ro{ib}.*(hfcase.meanFlow.u{ib}.*dsdx + ...
        hfcase.meanFlow.v{ib}.*dsdy);

    diss_prop = hfcase.meanFlow.diss_T{ib};
    irrev_prop = hfcase.meanFlow.irrev_gen{ib};
    rev_prop = dqx_Tdx + dqy_Tdy;

    xtmp{ib} = zeros(nentries, 1);
    ytmp{ib} = zeros(nentries, 1);
    areatmp{ib} = zeros(nentries, 1);
    %unsttmp{ib} = zeros(nentries,1);
    convtmp{ib} = zeros(nentries, 1);
    disstmp{ib} = zeros(nentries, 1);
    irrevtmp{ib} = zeros(nentries, 1);
    revtmp{ib} = zeros(nentries, 1);
    itmp{ib} = zeros(nentries, 1);
    jtmp{ib} = zeros(nentries, 1);
    ibtmp{ib} = zeros(nentries, 1);
    
    for i=1:celldims(ib,1)
        if mod(i, 50) == 0
             sprintf('Block %d, i=%d/%d', ib, i, celldims(ib,1))
        end
        for j=1:celldims(ib,2)
            pos = celldims(ib,2)*(i-1)+j;
            xnow = [hfcase.blk.x{ib}(i,j) hfcase.blk.x{ib}(i+1,j) ...
                hfcase.blk.x{ib}(i+1,j+1) hfcase.blk.x{ib}(i,j+1)];
            ynow = [hfcase.blk.y{ib}(i,j) hfcase.blk.y{ib}(i+1,j) ...
                hfcase.blk.y{ib}(i+1,j+1) hfcase.blk.y{ib}(i,j+1)];

            xtmp{ib}(pos) = mean(xnow);
            ytmp{ib}(pos) = mean(ynow);

            area = polyarea(xnow,ynow);
            areatmp{ib}(pos) = area;

          %  unsttmp{ib}(pos) = 0.25*(unst_prop(i,j)+unst_prop(i+1,j)+unst_prop(i+1,j+1)+unst_prop(i,j+1))*area;
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

for ib = 1:hfcase.NB

    offset = sum(prod(celldims(1:ib,:),2));
    nentries = celldims(ib,1)*celldims(ib,2);

    cells.x(offset-nentries+1:offset) = xtmp{ib};
    cells.y(offset-nentries+1:offset) = ytmp{ib};
    cells.area(offset-nentries+1:offset) = areatmp{ib};
   % cells.unst(offset-nentries+1:offset) = unsttmp{ib};
    cells.conv(offset-nentries+1:offset) = convtmp{ib};
    cells.diss(offset-nentries+1:offset) = disstmp{ib};
    cells.irrev(offset-nentries+1:offset) = irrevtmp{ib};
    cells.rev(offset-nentries+1:offset) = revtmp{ib};
    cells.i(offset-nentries+1:offset) = itmp{ib};
    cells.j(offset-nentries+1:offset) = jtmp{ib};
    cells.ib(offset-nentries+1:offset) = ibtmp{ib};
end


% for ir = 1:length(regions)
%     inds = cells.ib == regions{ir}.nb & ...
%         cells.i >= regions{ir}.is & cells.i <= regions{ir}.ie & ...
%         cells.j >= regions{ir}.js & cells.j <= regions{ir}.je;
%     e_unst(ir) = 0;%sum(cells.unst(inds));
%     e_s(ir) = sum(cells.conv(inds));
%     e_phi(ir) = sum(cells.diss(inds));
%     e_irrev(ir) = sum(cells.irrev(inds));
%     e_rev(ir) = sum(cells.rev(inds));
% end


%%

indices = cells.x > 0.75 & cells.x < 1.2 & cells.y > -0.1 & cells.y < 0.1;

pre_shock_inds = cells.ib == 6 & cells.i>= 100 & cells.i <= 400;
post_shock_inds = cells.ib == 6 & cells.i >= 475 & cells.i <= 672;
te_inds = cells.ib == 9;
wake_inds = cells.ib == 10  & cells.i <= 150;

%e_unst = sum(cells.unst(indices));
%%





e_s(1) = sum(cells.conv(pre_shock_inds));
e_phi(1) = sum(cells.diss(pre_shock_inds));
e_irrev(1) = sum(cells.irrev(pre_shock_inds));
e_rev(1) = sum(cells.rev(pre_shock_inds));

e_s(2) = sum(cells.conv(post_shock_inds));
e_phi(2) = sum(cells.diss(post_shock_inds));
e_irrev(2) = sum(cells.irrev(post_shock_inds));
e_rev(2) = sum(cells.rev(post_shock_inds));

e_s(3) = sum(cells.conv(te_inds));
e_phi(3) = sum(cells.diss(te_inds));
e_irrev(3) = sum(cells.irrev(te_inds));
e_rev(3) = sum(cells.rev(te_inds));

e_s(4) = sum(cells.conv(wake_inds));
e_phi(4) = sum(cells.diss(wake_inds));
e_irrev(4) = sum(cells.irrev(wake_inds));
e_rev(4) = sum(cells.rev(wake_inds));

e_unst(1:4) = 0;

for ng = 1:length(e_s)
    stackData(ng, 1, 1) = e_s(ng);
    stackData(ng, 1, 2) = e_unst(ng);
    stackData(ng, 1, 3:5) = 0;

    stackData(ng, 2, 1:2) = 0;
    stackData(ng, 2, 3) = e_phi(ng);
    stackData(ng, 2, 4) = e_irrev(ng);
    stackData(ng, 2, 5) = e_rev(ng);
end
% 
% 
% %%
% 
% groupLabels = {"Pre-shock", "Post-shock", "Triling edge","Wake"};
% plotBarStackGroups(stackData, groupLabels)