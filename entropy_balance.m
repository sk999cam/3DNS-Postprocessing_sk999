clear
close all




% initcase = DNS_case('cwl90_300k_tripped',2);
% initcase.readInstFlow;
% s1 = initcase.instFlow.k_ave('s');
% clear initcase
% 
% 
% hfcase = DNS_case('cwl90_300k_tripped',3);
% hfcase.readInstFlow;
% disp('read inst flows');
% s2 = hfcase.instFlow.k_ave('s');
% clear hfcase

hfcase = DNS_case('cwl90_300k_tripped',3);
hfcase.readMeanFlow;

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
end

%%

indices = cells.x > 0.75 & cells.x < 1.2 & cells.y > -0.1 & cells.y < 0.1;

%e_unst = sum(cells.unst(indices));
e_s = sum(cells.conv(indices));
e_phi = sum(cells.diss(indices));
e_irrev = sum(cells.irrev(indices));
e_rev = sum(cells.rev(indices));

