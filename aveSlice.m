classdef aveSlice < flowSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        p0in;
        Uinf;
        roinf;
        muinf;
        dsdyThresh = -50;
        use_unsflo = false;
        blEdgeMode;
    end

    properties (Dependent = true)
        Msurf;          % Surface Mach No
        theta;          % Momentum thickness
        theta_incomp;
        H;              % Shape factor
        H_incomp;       % Incompresisble shape factor
        %delta99;        % BL thickness
        %delta99_unsflo;
        delStar;        % Displacemnt thickness
        delta99_incomp; % BL thickness, incompressible definition
        delStar_incomp; % Displacemnt thickness, incompressible definition
        dsdy;           % Wall normal entropy gradient
        %BLedgeInd;      % j index of detected BL edge
        U;              % Wall-parallel velocity
        Res;            % Surface distance Reynolds No
        blPr;           % Componant of cd due to production of tke
        tau_w;
        dUdy;
        dTdy;
        cf;
        Re_theta;
    end

    methods
        function obj = aveSlice(blk, gas)
            obj@flowSlice(blk, gas);
            disp('Constructing aveSlice')
        end

        function getBCs(obj, inlet_blocks)
            Mnow = obj.M;
            Unow = obj.vel;
            ronow = obj.ro;
            munow = obj.mu;
            %Mnow = Mnow{inlet_blocks};
            pnow = obj.p;
            %pnow = pnow{inlet_blocks};
            
            p0 = [];
            Uinf = [];
            muinf = [];
            roinf = [];
            for i=1:length(inlet_blocks)
                p0now = pnow{inlet_blocks(i)}.*(1+((obj.gas.gam - 1)/2)*Mnow{inlet_blocks(i)}.^2).^(obj.gas.gam/(obj.gas.gam-1));
                p0 = [p0 p0now(40:100,:)];
                Uinf = [Uinf Unow{inlet_blocks(i)}(40:100,:)];
                muinf = [muinf munow{inlet_blocks(i)}(40:100,:)];
                roinf = [roinf ronow{inlet_blocks(i)}(40:100,:)];
            end
            obj.p0in = mean(p0,'all');
            obj.Uinf = mean(Uinf,'all');
            obj.muinf = mean(muinf,'all');
            obj.roinf = mean(roinf,'all');
        end

        function value = get.dsdy(obj)
            %s = obj.oGridProp('s');

            %value = (s(:,2:end)-s(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('s');
        end

        function value = get.dUdy(obj)
            u = obj.oGridProp('U');
            value = (u(:,2:end)-u(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('U');
        end 

        function value = get.dTdy(obj)
            t = obj.oGridProp('T');
            value = (t(:,2:end)-t(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('T');
        end 

        function value = get.Msurf(obj)
            disp('Calculating surface M')
            pnow = obj.oGridProp('p');
            value = sqrt((2/(obj.gas.gam - 1)) * ( (pnow(:,1)/obj.p0in).^(-(obj.gas.gam-1)/obj.gas.gam) - 1));
        end

        function value = BLedgeInd(obj, mode)
            if nargin < 2
                mode = "sThresh";
            end
            u = obj.oGridProp('U');
            switch mode
                case "sGradThresh"
                    temp = obj.dsdy;
                    for i=1:size(temp,1)
                        j=size(temp,2);
                        while temp(i,j) > obj.dsdyThresh && j>1
                            j = j-1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        [~, value(i)] = max(u(i,1:value(i)));
                    end
                case "sThresh"
                    temp = obj.oGridProp('s');
                    for i=1:size(temp,1)
                        j=size(temp,2);
                        se = temp(i,end);
                        sw = temp(i,1);
                        se;
                        while temp(i,j) < se + 0.02*(sw-se) && j>1
                            j = j-1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        %[~, value(i)] = max(u(i,1:value(i)));
                    end

                case "unsflo"
                    del99 = obj.delta99("unsflo");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        [~, value(i)] = max(u(i,1:value(i)));
                        %value(i) = min(size(temp,2), indnow+1);
                    end
                case "sInteg"
                    del99 = obj.delta99("sInteg");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        [~, value(i)] = max(u(i,1:value(i)));
                        %value(i) = min(size(temp,2), indnow+1);
                    end
                case "p0Integ"
                    del99 = obj.delta99("p0Integ");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        %value(i) = min(size(temp,2), indnow+1);
                    end
            end

        end

        function [xedge, yedge] = getBLedgeCoords(obj,mode)
            if nargin < 2
                mode = "sGradThresh";
            end
            inds = obj.BLedgeInd(mode);
            for i=1:length(inds)
                xedge(i) = obj.xO(i,inds(i));
                yedge(i) = obj.yO(i,inds(i));
            end
        end

        

        function value = delta99_unsflo(obj)
            %inds = obj.BLedgeInd;
            
        end

        function value = delta99(obj, mode)
            if nargin < 2
                mode = "sThresh";
            end
            switch mode
                case "sGradThresh"
                    inds = obj.BLedgeInd("sGradThresh");
                    for i=1:size(obj.yBL,1)
                        value(i) = obj.yBL(i,inds(i));
                    end
                case "unsflo"
                    ugrad = abs(obj.dUdy);
                    tgrad = abs(obj.dTdy);
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(ugrad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (ugrad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "sInteg"
                    sgrad = abs(obj.blNormGrad('s'));
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(sgrad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (sgrad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "p0Integ"
                    p0grad = abs(obj.blNormGrad('p0'));
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(p0grad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (p0grad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "sThresh"
                    inds = obj.BLedgeInd("sThresh");
                    for i=1:size(obj.yBL,1)
                        value(i) = obj.yBL(i,inds(i));
                    end
            end
            
        end


        function value = get.U(obj)
            disp('Calculating U')
            unow = obj.oGridProp('u');
            vnow = obj.oGridProp('v');
            nnow = obj.n;
            value = zeros(size(obj.yBL));
            for i=1:size(obj.yBL,1)
                tang = [-nnow(2,i); nnow(1,i)];
                for j=1:size(obj.yBL,2)
                    velnow = [unow(i,j); vnow(i,j)];
                    value(i,j) = -dot(tang, velnow - nnow(:,i)*dot(nnow(:,i),velnow));
                end
            end
        end

        function value = get.delStar(obj)
            inds = obj.BLedgeInd(obj.blEdgeMode);
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value = zeros(1,length(inds));
            for i=1:size(obj.yBL,1)
                integrand = 1 - ronow(i,1:inds(i)).*Unow(i,1:inds(i))./(ronow(i,inds(i))*Unow(i,inds(i)));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
        end

        function value = get.delStar_incomp(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            value = zeros(1,length(inds));
            for i=1:size(obj.yBL,1)
                integrand = 1 - Unow(i,1:inds(i))./Unow(i,inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
        end

        function value = get.theta(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof/(ro0*U0)).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
        end

        function value = get.theta_incomp(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                Uprof = Unow(i,1:inds(i));
                U0 = Unow(i,inds(i));
                integrand = (Uprof/U0).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
        end

        function value = get.H(obj)
            value = obj.delStar./obj.theta;
        end

        function value = get.H_incomp(obj)
            value = obj.delStar_incomp./obj.theta_incomp;
        end

        function value = get.Res(obj)
            value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
        end

        function plot_BL_profile(obj,x,prop,ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            size(q);
            if string(prop) == "dsdy"
                plot(ax, q, obj.yBL(i,:))
                hold on
                scatter(ax, q(obj.BLedgeInd(i)), obj.yBL(i,obj.BLedgeInd(i)))
            else
                plot(ax, q, obj.yBL(i,:)/obj.yBL(i,end))
                hold on
                scatter(ax, q(obj.BLedgeInd(i)), obj.yBL(i,obj.BLedgeInd(i))/obj.yBL(i,end))
            end
        end

        function blDevPlot(obj, prop, ax, lims, xrange, fmt)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);

            if nargin > 5 && ~isempty(fmt)
                plot(ax,obj.xSurf,q,fmt);
            elseif nargin>4 && ~isempty(xrange)
                plot(ax,obj.xSurf(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)),q(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)));
            else
                plot(ax,obj.xSurf,q);
            end
            if nargin > 3 && ~isempty(lims)
                ylim(lims);
            end
        end

        function plot_H_Pr_locus(obj, ax, ploteq, xrange, fmt)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end

            H = obj.H_incomp;
            pr = obj.blPr;
            if nargin>3 && ~isempty(xrange)
                H = H(obj.xSurf>xrange(1)&obj.xSurf<xrange(2));
                pr = pr(obj.xSurf>xrange(1)&obj.xSurf<xrange(2));
            end
            
            if nargin > 4 && ~isempty(fmt)
                plot(ax,H,pr,fmt);
            else
                plot(ax,H,pr);
            end
            
            if ploteq
                xtmp = linspace(1,3,51);% linspace(min(H),max(H),51);
                ytmp = 0.02456*((xtmp-1)./xtmp).^3;
                hold on
                plot(xtmp,ytmp,'k:')
            end
        end

        function value = get.blPr(obj)
            inds = obj.BLedgeInd;
            Prnow = obj.oGridProp('Pr');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                Prprof = Prnow(i,1:inds(i));
                Ue = Unow(i,inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, Prprof)/Ue^3;
            end
        end

        function value = get.Re_theta(obj)
            inds = obj.BLedgeInd;
            munow = obj.oGridProp('mu');
            
            ronow = obj.oGridProp('ro');
            size(ronow)
            size(inds)
            Unow = obj.U;
            thnow = obj.theta;
            for i=1:length(inds)
                value(i) = ronow(i,inds(i)).*Unow(i,inds(i)).*(thnow(i))./munow(i,inds(i));
            end

        end

        function value = get.cf(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:length(inds)
                roe(i) = ronow(i,inds(i));
                Ue(i) = Unow(i,inds(i));
            end
            value = obj.tau_w'./(0.5*roe.*Ue.*Ue);
        end

        function value = get.tau_w(obj)

            Unow = obj.U(:,2);
            Y0 = obj.yBL(:,2);
            munow = obj.oGridProp('mu');
            value = munow(:,2).*Unow./Y0;

        end

        function [xplus,yplus,zplus] = wall_coords_offset(obj)
            dy = obj.yBL(:,2);
            size(dy)
            ds = obj.ssurf(2:end) - obj.ssurf(1:end-1);
            ds(end+1) = ds(end);
            ds = ds';
            size(ds)
            dz = ones(size(dy))*obj.span/(obj.nk-1);

            munow = obj.oGridProp('mu');
            munow = munow(:,2);
            ronow = obj.oGridProp('ro');
            ronow = ronow(:,2);

            xplus = ds.*sqrt(abs(obj.tau_w).*ronow)./munow;
            yplus = dy.*sqrt(abs(obj.tau_w).*ronow)./munow;
            zplus = dz.*sqrt(abs(obj.tau_w).*ronow)./munow;
        end
        
        function blContour(obj, prop, ax, lims, fmt)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.oGridProp(prop);

            if nargin > 4 && ~isempty(fmt)
                pcolor(ax,obj.xO,obj.yO,q,fmt);
            else
                pcolor(ax,obj.xO,obj.yO,q);
            end
            shading interp
            axis equal
            if nargin > 3 && ~isempty(lims)
                caxis(lims);
            end

            cb = colorbar(ax,"southoutside");
            cb.Label.String = string(prop);

        end



        
    end
end
