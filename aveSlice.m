classdef aveSlice < flowSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        p0in;
        Uinf;
        roinf;
        muinf;
        dsdyThresh = -50;
    end

    properties (Dependent = true)
        Msurf;          % Surface Mach No
        theta;          % Momentum thickness
        H;              % Shape factor
        delta99;        % BL thickness
        delStar;        % Displacemnt thickness
        dsdy;           % Wall normal entropy gradient
        BLedgeInd;      % j index of detected BL edge
        U;              % Wall-parallel velocity
        Res;            % Surface distance Reynolds No
        blPr;           % Componant of cd due to production of tke
        tau_w;
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
            s = obj.oGridProp('s');
            value = (s(:,2:end)-s(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
        end       

        function value = get.Msurf(obj)
            disp('Calculating surface M')
            pnow = obj.oGridProp('p');
            value = sqrt((2/(obj.gas.gam - 1)) * ( (pnow(:,1)/obj.p0in).^(-(obj.gas.gam-1)/obj.gas.gam) - 1));
        end

        function value = get.BLedgeInd(obj)
            temp = obj.dsdy;
            for i=1:size(temp,1)
                j=size(temp,2);
                while temp(i,j) > obj.dsdyThresh && j>1
                    j = j-1;
                end
                value(i) = j+1;
            end
        end

        function value = get.delta99(obj)
            inds = obj.BLedgeInd;
            for i=1:size(obj.yBL,1)
                value(i) = obj.yBL(i,inds(i));
            end
        end

        function value = get.U(obj)
            disp('Calculating U')
            unow = obj.oGridProp('u');
            vnow = obj.oGridProp('v');
            nnow = obj.n;
            value = zeros(size(obj.yBL));
            for i=1:size(obj.yBL,1)
                for j=1:size(obj.yBL,2)
                    velnow = [unow(i,j); vnow(i,j)];
                    value(i,j) = norm(velnow - nnow(:,i)*dot(nnow(:,i),velnow));
                end
            end
        end

        function value = get.delStar(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value = zeros(1,length(inds));
            for i=1:size(obj.yBL,1)
                integrand = 1 - ronow(i,1:inds(i)).*Unow(i,1:inds(i))./(ronow(i,inds(i))*Unow(i,inds(i)));
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

        function value = get.H(obj)
            value = obj.delStar./obj.theta;
        end

        function value = get.Res(obj)
            value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
        end

        function blDevPlot(obj, prop, ax, lims, fmt)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);

            if nargin > 4 && ~isempty(fmt)
                plot(ax,obj.xSurf,q,fmt);
            else
                plot(ax,obj.xSurf,q);
            end
            if nargin > 3 && ~isempty(lims)
                ylim(lims);
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


        
    end
end