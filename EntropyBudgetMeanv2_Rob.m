% Calculate the entropy budget for a flow from a mean slice.

% v2 is for use on restructured mesh where y increases with j. All co-ords
% now increase with i,j and k respectively.
close all
clear all
fclose all

%% Set Directory and key variables
             % Location of slice
isst = 1;    % set = 1 if no slice in that direction (i.e we have all points)
jsst = 228;
ksst = 1;
             % Set slice size
nis  = 0;    % If 0, then nis=ni etc.
njs  = 25;
nks  = 0;

nibuf   = 10; %10; % seems to be a huge number at i = 30 for the irev heat transfer term
nibufex = 30;

             % Stats
nstats = 23;
nmean_start = 2;
nmean_end   = 2;
             % Gas properties
cp   = 1005.0;
gam  = 1.4;
cv   = cp/gam;
rgas = cp-cv;
Pr   = 0.72;
mu_c = 3.42E-5 ;

%dir = '/home/rgp32/rds/rds-rr-3dtbl-LlPtkm5bqO8/3DNS_DuctExp/2022-09-07-UnderstandingParrelel/'
dir = '/home/rgp32/rds/rds-rr-3dtbl-LlPtkm5bqO8/3DNS_DuctExp/2022-09-30-Turbulent_finer400k_root3/'
dir = '/home/rgp32/rds/rds-rr-3dtbl-LlPtkm5bqO8/3DNS_DuctExp/2022-10-04-Turbulent_finer400k_restructured/'

% Figure formatting
fsfrac = 24;
fs     = 20; % Font size
lw     = 1.5; % Line wdith
figures_dir = [dir, 'Figures_mean/'];

%% Load Geometry
span_1 = load([dir,'span_1.txt']);
span = span_1(:,1);
load([dir,'blockdims.txt']);
Nb = size(blockdims,1)
for nb=1:Nb
    fid = fopen([dir,'grid_',num2str(nb),'.txt'],'r');
    ni = blockdims(nb,1);
    nj = blockdims(nb,2);
    nk = blockdims(nb,3);
    x = zeros(ni,nj);
    y = zeros(ni,nj);
    for j=1:nj
        for i=1:ni
            [A]=fscanf(fid, '%f %f',2);
            x(i,j) = A(1);
            y(i,j) = A(2);
        end
    end
    fclose(fid);
    zz(1,1,1:nk) = span;
    z = repmat(zz,[ni nj 1]);
    x = repmat(x,[1 1 nk]);
    y = repmat(y,[1 1 nk]);
end
njmid = floor(njs./2);

% Set slice size
if nis == 0; nis=ni; end
if njs == 0; njs=nj; end
if nks == 0; nks=nk; end

% Get slice coords
xc = x(isst:isst+nis-1, jsst:jsst+njs-1, ksst:ksst+nks-1);
yc = y(isst:isst+nis-1, jsst:jsst+njs-1, ksst:ksst+nks-1);
zc = z(isst:isst+nis-1, jsst:jsst+njs-1, ksst:ksst+nks-1);
sl = cumsum(sqrt( diff(xc,1,1).^2 + diff(yc,1,1).^2 ),1);
sl = vertcat(zeros(1,njs,nks),sl);
%% mean
load([dir,'mean_time.txt']) 
tot_time = sum(mean_time(nmean_start:nmean_end,3));


%% Read mean file in
for nmean=nmean_start:nmean_end
    delt_now = mean_time(nmean,3);

    % Create data arrays
    if(nmean==nmean_start)
        ro = zeros(nis,njs,nks);
        ru = zeros(nis,njs,nks);
        rv = zeros(nis,njs,nks);
        rw = zeros(nis,njs,nks);
        Et = zeros(nis,njs,nks);
        ro2 = zeros(nis,njs,nks);

        ruu =  zeros(nis,njs,nks);
        rvv =  zeros(nis,njs,nks);
        rww =  zeros(nis,njs,nks);
        ruv =  zeros(nis,njs,nks);
        ruw =  zeros(nis,njs,nks);
        rvw =  zeros(nis,njs,nks);

        rus =  zeros(nis,njs,nks);
        rvs =  zeros(nis,njs,nks);
        rws =  zeros(nis,njs,nks);

        diss = zeros(nis,njs,nks);
        dissT = zeros(nis,njs,nks);
        
        qxT = zeros(nis,njs,nks);
        qyT = zeros(nis,njs,nks);
        qzT = zeros(nis,njs,nks);

        qx2T2 = zeros(nis,njs,nks);
        qy2T2 = zeros(nis,njs,nks);
        qz2T2 = zeros(nis,njs,nks);
        %snow = zeros(nis,njs,nks);
        %pnow = zeros(nis,njs,nks);
        %tnow = zeros(nis,njs,nks);
        %ponow = zeros(nis,njs,nks);
        %tonow = zeros(nis,njs,nks);
    end

    % Read in mean file
    % Check the nmean is correct (the write out of mean_time was initially
    % messed up - mean_count was not being updated in the time.txt file)
    fid2 = fopen([dir,'meanslice2_',num2str(nb),'_',num2str(nmean)],'r');
    A = fread(fid2,inf,'float64');
    A = reshape(A,nstats,length(A)/nstats);
    fclose(fid2)
% mnods     _1_1

    fid3 = fopen([dir,'mnodslice2_',num2str(nb),'_',num2str(nmean)],'r');
    B = fread(fid3,inf,'int');
    disp(length(B))
    B = reshape(B,3,[]);
    fclose(fid3);

    NN = size(B,2);
    if size(A,2) < NN
        NN = size(A,2);
    end
    icount(1:nis*njs*nks) = 0;
    for n=1:NN;
        i = B(1,n) -isst+1;
        j = B(2,n) -jsst+1;
        k = B(3,n) -ksst+1;
        nid = i + (j-1)*nis + (k-1)*nis*njs;
        if(i<=nis & j<=njs & k<=nks & icount(nid)==0)
            ro(i,j,k) = ro(i,j,k) + A(1,n)./tot_time;
            ru(i,j,k) = ru(i,j,k) + A(2,n)./tot_time;
            rv(i,j,k) = rv(i,j,k) + A(3,n)./tot_time;
            rw(i,j,k) = rw(i,j,k) + A(4,n)./tot_time;
            Et(i,j,k) = Et(i,j,k) + A(5,n)./tot_time;

            ro2(i,j,k) = ro2(i,j,k) + A(6,n)./tot_time;
            ruu(i,j,k) = ruu(i,j,k) + A(7,n)./tot_time;
            rvv(i,j,k) = rvv(i,j,k) + A(8,n)./tot_time;
            rww(i,j,k) = rww(i,j,k) + A(9,n)./tot_time;
            ruv(i,j,k) = ruv(i,j,k) + A(10,n)./tot_time;
            ruw(i,j,k) = ruw(i,j,k) + A(11,n)./tot_time;
            rvw(i,j,k) = rvw(i,j,k) + A(12,n)./tot_time;

            rus(i,j,k) = rus(i,j,k) + A(13,n)./tot_time;
            rvs(i,j,k) = rvs(i,j,k) + A(14,n)./tot_time;
            rws(i,j,k) = rws(i,j,k) + A(15,n)./tot_time;

            diss(i,j,k) = diss(i,j,k)  + A(16,n)./tot_time; % Viscous dissipation (total)
            dissT(i,j,k)= dissT(i,j,k) + A(17,n)./tot_time; % Viscous dissipation (total) / T

            % Reversible heat transfer terms
            qxT(i,j,k) = qxT(i,j,k) + A(18,n)./tot_time; % qx/T
            qyT(i,j,k) = qyT(i,j,k) + A(19,n)./tot_time; % qy/T
            qzT(i,j,k) = qzT(i,j,k) + A(20,n)./tot_time; % qz/T

            % Irreversible heat transfer terms
            qx2T2(i,j,k) = qx2T2(i,j,k) + A(21,n)./tot_time; % qx^2/T^2
            qy2T2(i,j,k) = qy2T2(i,j,k) + A(22,n)./tot_time; % qy^2/T^2
            qz2T2(i,j,k) = qz2T2(i,j,k) + A(23,n)./tot_time; % qz^2/T^2
        end
        icount(nid) = icount(nid)+1;
        %icount(nid) = 1;
    end
end

[drsdt, ~] = calc_drsdt(dir, nb, tot_time);
drsdt = drsdt(isst:isst+nis-1, jsst:jsst+njs-1, ksst:ksst+nks-1);
drsdt = drsdt./tot_time;

%% Calculate Entropy budget

% % Use fliplr to flip grid so that y increases with j
% % At the moment the bend curves upwards (increasing y co-ord) but the inner
% % radius (greater y) is the j=1 line. Therefore as we increase j, y
% % actually decreases. Need to use fliplr to correct this.
% xc = fliplr(xc);
% yc = fliplr(yc);
% zc = fliplr(zc);
% drsdt = fliplr(drsdt);
% dissT = fliplr(dissT);
% qx2T2 = fliplr(qx2T2);
% qy2T2 = fliplr(qy2T2);
% qz2T2 = fliplr(qz2T2);
% rus = fliplr(rus);
% rvs = fliplr(rvs);
% rws = fliplr(rws);
% qxT = fliplr(qxT);
% qyT = fliplr(qyT);
% qzT = fliplr(qzT);


% Now proceed
% First calculate the areas of the faces require for finite vol derivatives
gridVar = finVolDeriv3D_grid(xc,yc,zc);

% Volume integrals
dissTV = trapz(squeeze(zc(1,1,:)), cellAvg(dissT).*gridVar{5}, 3);
qx2T2V = trapz(squeeze(zc(1,1,:)), cellAvg(qx2T2).*gridVar{5}, 3);
qy2T2V = trapz(squeeze(zc(1,1,:)), cellAvg(qy2T2).*gridVar{5}, 3);
qz2T2V = trapz(squeeze(zc(1,1,:)), cellAvg(qz2T2).*gridVar{5}, 3);
rsV    = trapz(squeeze(zc(1,1,:)), cellAvg(drsdt).*gridVar{5}, 3);

% Entropy Fluxes
[drusdx, ~, ~] = finVolDeriv3D_nG(gridVar,rus(:,:,:));
[~, drvsdy, ~] = finVolDeriv3D_nG(gridVar,rvs(:,:,:));
[~, ~, drwsdz] = finVolDeriv3D_nG(gridVar,rws(:,:,:));
[dqxTdx, ~, ~] = finVolDeriv3D_nG(gridVar,qxT(:,:,:));
[~, dqyTdy, ~] = finVolDeriv3D_nG(gridVar,qyT(:,:,:));
[~, ~, dqzTdz] = finVolDeriv3D_nG(gridVar,qzT(:,:,:));

% Integrate fluxes through BL
drusdxV  = trapz(squeeze(zc(1,1,:)), drusdx.*gridVar{5}, 3);
drvsdyV  = trapz(squeeze(zc(1,1,:)), drvsdy.*gridVar{5}, 3);
drwsdzV  = trapz(squeeze(zc(1,1,:)), drwsdz.*gridVar{5}, 3);
dqxTdxV  = trapz(squeeze(zc(1,1,:)), dqxTdx.*gridVar{5}, 3);
dqyTdyV  = trapz(squeeze(zc(1,1,:)), dqyTdy.*gridVar{5}, 3);
dqzTdzV  = trapz(squeeze(zc(1,1,:)), dqzTdz.*gridVar{5}, 3);

% Group Terms
Ent_conv = drusdxV+drvsdyV+drwsdzV;
Ent_qrev = dqxTdxV+dqyTdyV+dqzTdzV;
Ent_qirr = qx2T2V+qy2T2V+qz2T2V;

linIdx = find(abs(Ent_qirr) > 1000);
[rows,columns] = ind2sub(size(Ent_qirr), linIdx);
Ent_qirr(rows,columns) = 0;

numerics = rsV + Ent_conv - dissTV - Ent_qirr + Ent_qrev;

% Calculate over whole volume
Ent_convVol = sum(Ent_conv(nibuf:end-nibufex,:), 'all');
Ent_rsVol   = sum(rsV(nibuf:end-nibufex,:), 'all');
Ent_lhsVol  = Ent_convVol + Ent_rsVol;
Ent_qrevVol = sum(Ent_qrev(nibuf:end-nibufex,:), 'all');
Ent_qirrVol = sum(Ent_qirr(nibuf:end-nibufex,:), 'all');
Ent_dissVol = sum(dissTV(nibuf:end-nibufex,:), 'all');
NumericsVol = sum(numerics(nibuf:end-nibufex,:), 'all');

% Plot Entropy along duct
figure()
plot(sl(1:end-1,njmid,1), Ent_conv(:,1)+rsV(:,1), 'k-'), hold on
plot(sl(1:end-1,njmid,1), dissTV(:,1), 'r--')
plot(sl(1:end-1,njmid,1), dissTV(:,1) + Ent_qirr(:,1), 'b--')
plot(sl(1:end-1,njmid,1), dissTV(:,1) + Ent_qirr(:,1) - Ent_qrev(:,1), 'k--')

% Plot local entropy budget
figure()
b = bar([Ent_dissVol./Ent_lhsVol,...
    Ent_qirrVol./Ent_lhsVol,...
    -Ent_qrevVol./Ent_lhsVol,...
    NumericsVol./Ent_lhsVol;...
    0,0,0,0], 'stacked');
xlim([0,2])
set(gca,'XTick',[])
set(b, {'DisplayName'}, {'Viscous Dissipation', 'Irreversible Heat Transfer', 'Reversible Heat Transfer', 'Numerics'}')
legend('Location','eastoutside')
clear b
saveas(gcf,([figures_dir,'EntropyBudget.pdf']))
% 
% % Integrate entropy budget through BL
% %edge = nk;
% % Find edege of BL
% velCellAvg2 = cellAvg(cellAvg(vel));
% UinfCellAvg2 = mean(velCellAvg2(:,njmidcell2,nkmid:end));
% % Integrate entropy through BL
% for i = 1:ni
%     edge = k_fs999(i);
%     zvec = squeeze(z(1,njmid,1:edge));
%     sconv_tot(i) = trapz(zvec, squeeze(sconv(1,njmidcell2,1:edge)));
%     sdiss_tot(i) = trapz(zvec, squeeze(sdiss(1,njmidcell2,1:edge)));
%     sirr_tot(i)  = trapz(zvec, squeeze(sirr(1,njmidcell2,1:edge)));
%     srev_tot(i)  = trapz(zvec, squeeze(srev(1,njmidcell2,1:edge)));
%     numerics_tot(i) = trapz(zvec, squeeze(numerics(1,njmidcell2,1:edge)));
% end
% % Plot total entropy budget
% % figure()
% % b = bar([sdiss_tot./sconv_tot, sirr_tot./sconv_tot, -srev_tot./sconv_tot, numerics_tot./sconv_tot; 0,0,0,0], 'stacked');
% % xlim([0,2])
% % set(gca,'XTick',[])
% % set(b, {'DisplayName'}, {'Viscous Dissipation', 'Irreversible Heat Transfer', 'Reversible Heat Transfer', 'Numerics'}')
% % legend('Location','eastoutside')
% figure()
% plot(Ret, sdiss_tot./sconv_tot, 'b'), hold on
% plot(Ret, (sdiss_tot+sirr_tot)./sconv_tot, 'r')
% plot(Ret, (sdiss_tot+sirr_tot-srev_tot)./sconv_tot, 'o')
% plot(Ret, (sdiss_tot+sirr_tot-srev_tot+numerics_tot)./sconv_tot, 'k')

return

%% debugging when Ent_qirr is stupidly large
stupidcells = Ent_qirr(abs(Ent_qirr) > 1000);
linIdx = find(abs(Ent_qirr) > 1000);
[rows,columns] = ind2sub(size(Ent_qirr), linIdx);
