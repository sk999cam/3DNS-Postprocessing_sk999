classdef dnsProbe < handle
    % DNSPROBE Class containing 3DNS probe data and metadata
    %   

    properties
        nb;
        i;
        j;
        k;
        x;
        y;
        t;
        ro;
        u;
        v;
        w;
        Et;
        gas;
        nSkip;
        nSamples;
    end

    properties (Dependent = true)
        tke;
        Tu;
    end

    methods
        function obj = dnsProbe(rundir, nProbe, nSkip, blkData, gas)
            obj.nSkip = nSkip;
            obj.x = blkData.x;
            obj.y = blkData.y;
            obj.nb = blkData.nb;
            obj.i = blkData.i;
            obj.j = blkData.j;
            obj.k = blkData.k;
            obj.gas = gas;
            
            fid = fopen(fullfile(rundir,['probe_' num2str(nProbe)]),'r');
            A = fscanf(fid, '%f %f %f %f %f %f', [6 inf]);
            fclose(fid);
            
            obj.t = A(1,:);
            obj.ro = A(2,:);
            obj.u = A(3,:)./obj.ro;
            obj.v = A(4,:)./obj.ro;
            obj.w = A(5,:)./obj.ro;
            obj.Et = A(6,:);

            obj.nSamples = length(obj.t);
        end

        function value = get.tke(obj)
            udash = obj.u-mean(obj.u);
            vdash = obj.v-mean(obj.v);
            wdash = obj.w-mean(obj.w);
            value = 0.5*(udash.^2 + vdash.^2 + wdash.^2);
        end

        function value = get.Tu(obj)
            %TURB_INTENSITY 
            %   Calculate turbulence intensuty
            
            U = mean(sqrt(obj.u.^2 + obj.v.^2 + obj.w.^2));
            value = sqrt(2*mean(obj.tke)/3)/U;

        end

        function plot_spectrum(obj,ax)
            if nargin<2 || isempty(ax)
                ax = gca;
            end
            N = length(obj.tke);
            Ts = (obj.t(end)-obj.t(1))/N;
            fs = 1/Ts;
            Y = fft(obj.tke);
            P2 = abs(Y/N);
            P1 = P2(1:N/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            psd = P1/fs;
            f = fs*(0:(N/2))/N;

            loglog(f, psd);
            grid on
            xlabel('f (Hz)')
            ylabel('PSD (dB/Hz)')
            set(gca,'fontSize',12)
        end
    end
end