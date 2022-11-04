fclose all
close all
clear all

% Script to read in probe data and calculate the autocorrelation and TKE
% spectra

% Figures are saved to a file called 'ProbeFigs' within the flow directory

%% Setup
% Hardcoded flow variables
cp = 1005.0;
gam = 1.4;
cv = cp/gam;
rgas = cp-cv;
mu_ref =  3.42E-5; 
nu   = mu_ref./1.225;
% Approximate Lturb 
Lturb = 0.1;     % Set Lturb equal to span
Lturb = Lturb./6;   % Set Lturb equal to BL thickness

% Plot formatting
lw = 1.5; % Linewidth

% Flow directory where probe files are stored
dir = '/run/media/ojj23/Crucial X6/isolated_aerofoils/r150_cwl_90/g1_sigma4/run3/';

% code will read probe data between tstart and tend
tstart = 0; 
tend   = 'end'; % Set to 'end' to read to the end of probe file

% Save figures to a folder called ProbeFigs within the flow directory
figures_dir = [dir, 'ProbeFigs/'];

%% Read probe data
% Load number of probes
fid = fopen([dir,'probe.txt'],'r');
[A] = fscanf(fid, '%d %d',2);
nprobe = A(1);
fclose(fid)

% Loop through probes
for np = 1:nprobe
    % Load Probe Data
    probe = load([dir,'probe_',num2str(np)],'r');
    [~, idx_tstart] = min( abs(probe(:,1)-tstart) );
    if tend ~= 'end'
        [~, idx_tend] = min( abs(probe(:,1)-tend) );
    else
        idx_tend = length(probe(:,1));
    end
    time{np} = probe(idx_tstart:idx_tend,1);
    time{np} = time{np} - time{np}(1);
    T{np}    = mean(diff(time{np}));
    Fs{np}   = 1./T{np};
    ro{np}   = probe(idx_tstart:idx_tend,2);
    ru{np}   = probe(idx_tstart:idx_tend,3);
    rv{np}   = probe(idx_tstart:idx_tend,4);
    rw{np}   = probe(idx_tstart:idx_tend,5);
    Et{np}   = probe(idx_tstart:idx_tend,6);

    % Calculate Properties
    u{np}    = ru{np}./ro{np};
    v{np}    = rv{np}./ro{np};
    w{np}    = rw{np}./ro{np};
    p{np}    = (gam-1)*(Et{np} - 0.5*(u{np}.^2 + v{np}.^2 + w{np}.^2).*ro{np});
    T{np}    = p{np}./(ro{np}*rgas);

    % Find mean and fluctuating velocities
    umean{np}= mean(u{np});
    vmean{np}= mean(v{np});
    wmean{np}= mean(w{np});
    vel{np}  = sqrt( umean{np}.^2 + vmean{np}.^2 + wmean{np}.^2 );
    ufluc{np}= u{np} - umean{np};
    vfluc{np}= v{np} - vmean{np};
    wfluc{np}= w{np} - wmean{np};
    
    % Calculate TKE
    TKE(np)    = 0.5.* sqrt(mean(ufluc{np}.^2) +  mean(vfluc{np}.^2) + mean(wfluc{np}.^2));
    uturb(np)  = sqrt(2.*TKE(np)./3);
    diss(np)   = uturb(np).^3./Lturb;
    Lk(np)     = (nu.^3./diss(np)).^0.25;
    tk(np)     = (nu./diss(np)).^0.5;
    uk(np)     = (nu.*diss(np)).^0.25;

    % Temporal Autocorrelation
    R_uu{np} = (cumsum((ufluc{np}(1).*ufluc{np}))./((1:length(ufluc{np})).') ./ sqrt(ufluc{np}(1).^2 .* cumsum(ufluc{np}.^2)./((1:length(ufluc{np})).')));
    R_vv{np} = (cumsum((vfluc{np}(1).*vfluc{np}))./((1:length(vfluc{np})).') ./ sqrt(vfluc{np}(1).^2 .* cumsum(vfluc{np}.^2)./((1:length(vfluc{np})).')));
    R_ww{np} = (cumsum((wfluc{np}(1).*wfluc{np}))./((1:length(wfluc{np})).') ./ sqrt(wfluc{np}(1).^2 .* cumsum(wfluc{np}.^2)./((1:length(wfluc{np})).')));

    R_uu_S{np} = cumsum((ufluc{np}(1).*ufluc{np}))./((1:length(ufluc{np})).');
    L{np}      = length(R_uu_S{np});

    % Find integral timescale:
    T_int(np,1) = trapz(time{np}, abs(R_uu{np}));
    T_int(np,2) = trapz(time{np}, abs(R_vv{np}));
    T_int(np,3) = trapz(time{np}, abs(R_ww{np}));

    % kolmogorov and turb frequencies
    fk(np)     = 1./tk(np);
    fL(np)     = 1./max(T_int(np,:));

    % Plot autocorrelation
    figure()
    hold on
    plot(time{np}, R_uu{np}, 'k-')
    plot(time{np}, R_vv{np}, 'b-')
    plot(time{np}, R_ww{np}, 'r-')
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize',20)
    ylabel('$R_{uu}$', 'Interpreter','latex','FontSize',20)
    title(['Probe ',num2str(np)], 'Interpreter','latex','FontSize',20)
    ylim([-1 1])
    legend('$R_{uu}$', '$R_{vv}$','$R_{ww}$', 'Interpreter','latex', 'Fontsize',18, 'location', 'eastoutside')
    grid on
    txt = ['$T_{int}$= ', num2str(round(T_int(np,1),2,"significant")), 's'];
    text(time{np}(floor(length(time{np})/2)), -0.8, txt, 'Interpreter','latex','FontSize',16 )
%    saveas(gcf,([figures_dir, 'R_uu_Probe', num2str(np),'.pdf']))
end

%return
%% FFT to find energy spectrum
for np = 1:nprobe


    % Pwelch method
    %[pu{np}, fu{np}] = pwelch(ufluc{np}, length(ufluc{np}), 0, [], 'psd', Fs{np}); % Pwelch with no windowing
    
    windowL = max(T_int(np,:)) .* Fs{np};
    [pu{np}, fu{np}] = pwelch(ufluc{np}, windowL, windowL./2, [], Fs{np}, 'psd'); % pwelch with windowing equal to the integral time
    %[pu{np}, fu{np}] = pwelch(ufluc{np}, [], [], [], Fs{np}, 'psd'); %pwelch with default windowing

    % Convert from frequency domain to wavenumber space using Taylor's
    % hypothesis of frozen turbulence
    ku{np} = fu{np}.*2.*pi./vel{np};
    Eu{np} = pu{np}.*vel{np};

    % Inertial sub range -5/3 line
    B = 10000; % B is a fudge to adjust the y intercept
    fit_wp = logspace(log10(min(ku{np}(2:end))),6,100);
    fit_Epw{np} = B.*fit_wp.^(-5/3);

    % Plot spectrum, non-dimming large scales
    figure()
    %loglog(ku{np},Eu{np}, 'k', "LineWidth",lw), hold on % dimensional version
    loglog(ku{np}.* Lturb, Eu{np}./(Lturb.*uturb(np).^2), 'k', "LineWidth",lw), hold on  % Non dim version
    % -5/3 inertial subrange
    %loglog(fit_wp.* Lurb, fit_Epw{np}./(Lturb.*uturb(np).^2), 'k--', 'LineWidth',lw-0.5) % Dimensional version
    loglog(fit_wp.* Lturb, fit_Epw{np}./(Lturb.*uturb(np).^2), 'k--', 'LineWidth',lw-0.5) % non dim version
    ylabel('$\frac{E_u}{L_{turb} u\prime^2}$', 'Interpreter', 'latex', 'FontSize',22)
    xlabel('$\kappa L_{turb}$', 'Interpreter','latex','FontSize',22)
    %xlim([1./Lturb, 10./Lk(np)]) % dimmed xlims
    xlim([10^-1, 10^2.5]) % non dim xlims
    title(['Probe ',num2str(np), ' - TKE Spectrum, Large Scales'], 'Interpreter','latex','FontSize',20)
    grid on
    set(gcf, 'Position', [100,100,650,1200])
%    saveas(gcf,([figures_dir, 'Spectrum_pwelch_LargeScales_Probe', num2str(np),'.pdf']))

    % Plot spectrum non-dimming small scales
    figure()
    %loglog(ku{np},Eu{np}, 'k', "LineWidth",lw), hold on % dimensional version
    loglog(ku{np}.* Lk(np), Eu{np}./(Lk(np).*uk(np).^2), 'k', "LineWidth",lw), hold on  % Non dim version
    % -5/3 inertial subrange
    %loglog(fit_wp.* Lurb, fit_Epw{np}./(Lturb.*uturb(np).^2), 'k--', 'LineWidth',lw-0.5) % Dimensional version
    loglog(fit_wp.* Lk(np), fit_Epw{np}./(Lk(np).*uk(np).^2), 'k--', 'LineWidth',lw-0.5) % non dim version
    ylabel('$\frac{E_u}{L_k u_k^2} $', 'Interpreter', 'latex', 'FontSize',22)
    xlabel('$\kappa L_k$', 'Interpreter','latex','FontSize',22)
    %xlim([1./Lturb, 10./Lk(np)]) % dimmed xlims
    xlim([10^-2.5, 10^1]) % non dim xlims
    title(['Probe ',num2str(np), ' - TKE Spectrum, Small Scales'], 'Interpreter','latex','FontSize',20)
    grid on
    set(gcf, 'Position', [100,100,650,1200])
%    saveas(gcf,([figures_dir, 'Spectrum_pwelch_SmallScales_Probe', num2str(np),'.pdf']))
end
%%
return



