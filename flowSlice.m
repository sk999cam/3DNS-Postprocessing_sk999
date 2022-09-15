classdef flowSlice < handle
    %FLOWSLICE Generic class containing flow on a k slice

    properties
        NB;
        gas;
        ro;
        u;
        v;
        w;
        Et;
        xSurf;
        yBL;
        xO;
        yO;
        oblocks;
        oblocks_flip;
        blk;
        iLE;
        iTE;
        n;
        ssurf;          % Surface distance fron LE
        vortZ;          % Z vorticity
    end

    properties (Dependent = true)
        %T;              % Temperature
        %p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        mu;             % Viscosity
        nu;             % Kinematic viscosity
        p0;
        schlieren;      % |grad(ro)|/ro
        cellSize;
    end

    methods
        function obj = flowSlice(blk, gas)
            %FLOWSLICE Construct a flowSlice object
            disp('Constructing flowSlice')
            if nargin > 0
                if size(blk.blockdims,1) > 1
                    obj.oblocks = blk.oblocks;
                    obj.oblocks_flip = blk.oblocks_flip;
                    obj.gas = gas;
                    obj.blk = blk;
                    %obj.gas.mu_ref = 0.0008748708280693193;
                    obj.gas.rgas = obj.gas.cp*(1-1/obj.gas.gam);
                    obj.NB = size(blk.blockdims,1);
                    xo = [];
                    yo = [];
                    for i=1:length(obj.oblocks)
                        xtmp = blk.x{obj.oblocks(i)}(:,:);
                        ytmp = blk.y{obj.oblocks(i)}(:,:);
                        xtmp = flip(xtmp,2);
                        ytmp = flip(ytmp,2);
                        if obj.oblocks_flip(i) == 1
                            xtmp = flip(xtmp);
                            ytmp = flip(ytmp);
                        end
                        xo = [xo; xtmp(1:end-1,:)];
                        yo = [yo; ytmp(1:end-1,:)];
                    end
                    xsurf = xo(:,1);
                    [~, obj.iLE] = min(xsurf);
                    [~, obj.iTE] = max(xsurf);
                    obj.xSurf = xsurf(obj.iLE:obj.iTE);
                    obj.xO = xo(obj.iLE:obj.iTE,:);
                    obj.yO = yo(obj.iLE:obj.iTE,:);

                    %obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
                    %size(obj.xSurf)
                    R = [0 -1; 1 0];

                    obj.yBL = zeros(obj.iTE-obj.iLE+1,size(xo,2));
                    obj.ssurf = zeros(1,size(obj.yBL,1));
                    for i = obj.iLE:obj.iTE

                        s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                        s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                        n1 = R*s1/norm(s1);
                        n2 = R*s2/norm(s2);
                        nnow = (0.5*(n1+n2)/norm(0.5*(n1+n2)));
                        obj.n(:,i+1-obj.iLE) = nnow;
                        for j=1:size(xo,2)
                            dx = xo(i,j) - xo(i,1);
                            dy = yo(i,j) - yo(i,1);
                            obj.yBL(i+1-obj.iLE,j) = dot(nnow, [dx;dy]);
                        end
                        if i>obj.iLE
                            dx = xo(i,1) - xo(i-1,1);
                            dy = yo(i,1) - yo(i-1,1);
                            ds = sqrt(dx^2 + dy^2);
                            obj.ssurf(i+1-obj.iLE) = obj.ssurf(i-obj.iLE) + ds;
                        end
                    end
                else
                    obj.gas = gas;
                    obj.gas.rgas = obj.gas.cp*(1-1/obj.gas.gam);
                    obj.NB = size(blk.blockdims,1);
                end
            end
        end         % End of constructor




        function value = get.vel(obj)
            disp('Calculating vel')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.p0(obj)
            pnow = obj.p;
            Mnow = obj.M;
            for nb = 1:obj.NB
                value{nb} = pnow{nb}.*(1+0.5*(obj.gas.gam-1)*Mnow{nb}.^2).^(obj.gas.gam/(obj.gas.gam-1));
            end
        end

        function value = get.M(obj)
            disp('Calculating M')
            value = cell(1,obj.NB);
            pnow = obj.p;
            for nb = 1:obj.NB
                %pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow{nb}./(obj.ro{nb}*obj.gas.rgas);
                velnow = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
                value{nb} = velnow./sqrt(obj.gas.gam*obj.gas.rgas*Tnow);
            end
        end

        function value = get.s(obj)
            disp('Calculating s')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.cp*log(Tnow/300) - obj.gas.rgas*log(pnow/1e5);
            end
        end

        function value = get.mu(obj)
            disp('Calcualting mu')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = obj.p{nb};%(obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.mu_ref*(Tnow/obj.gas.mu_tref).^(3/2) .* (obj.gas.mu_cref + obj.gas.mu_tref)./(obj.gas.mu_cref + Tnow);
            end
        end

        function value = get.nu(obj)
            munow = obj.mu;
            ronow = obj.ro;
            for nb = 1:obj.NB
                value{nb} = munow{nb}./ronow{nb};
            end
        end

        

        function value = get.cellSize(obj)
            fprintf('Calculating Cell Sizes')
            dz = obj.blk.span/(obj.blk.nk{1}-1);
            
            value = {};
            for ib = 1:obj.NB
                ni = size(obj.blk.x{ib},1);
                nj = size(obj.blk.x{ib},2);
                area = zeros(ni-1, nj-1);
                for i=1:ni-1
                    for j=1:nj-1
                        xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
                            obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
                        ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
                            obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];
                        area(i,j) = abs(polyarea(xnow,ynow));
                    end
                end
                area = dz*area;
                value{ib}(1,1) = area(1,1);
                value{ib}(1,nj) = area(1,nj-1);
                value{ib}(ni,1) = area(ni-1,1);
                value{ib}(ni,nj) = area(ni-1,nj-1);
                for i = 2:ni-1
                    value{ib}(i,1) = 0.5*(area(i-1,1)+area(i,1));
                    value{ib}(i,end) = 0.5*(area(i-1,end)+area(i,end));
                end
                for j = 2:nj-1
                    value{ib}(1,j) = 0.5*(area(1,j-1)+area(1,j));
                    value{ib}(end,j) = 0.5*(area(end,j-1)+area(end,j));
                end
                for i=2:ni-1
                    for j=2:nj-1
                        value{ib}(i,j) = 0.25*(area(i-1,j-1)+area(i-1,j)+area(i,j-1) +area(i,j));
                    end
                end
                value{ib} = value{ib}.^(1/3);
            end
        end

            


        function value = get.vortZ(obj)
            disp('Calculating z componant of vorticity')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                [~, dudy] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.u{nb});
                [dvdx, ~] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.v{nb});
                value{nb} = dvdx - dudy;
            end
        end

        function value = get.schlieren(obj)
            disp('calculating grad(ro)/ro')
            value = cell(1,obj.NB);
            for nb=1:obj.NB
                [drodx, drody] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.ro{nb});
                value{nb} = sqrt(drodx.^2 + drody.^2)./obj.ro{nb};
            end
        end

            
        

        function [profile, i] = BLprof(obj, x, prop)
            %BLPROF get a profile of prop across the BL at specified x
            [~, i] = min(abs(obj.xSurf-x));
            if any(strcmp(["dsdy" "U"],prop))
                profile = obj.(prop)(i,:);
                size(profile)
            else
                propfield = oGridProp(obj,prop);
                profile = propfield(i,:);
            end
        end

        function blfield = oGridProp(obj, prop)
            %OGRIDPROP Construct array of property in o grid for one surface
            %of blade
            
            if any(strcmp(["U","dsdy"], prop))
                blfield = obj.(prop);
            else
                blfield = [];
                propnow = obj.(prop);
                for i=1:length(obj.oblocks)
                    clear temp
                    temp = propnow{obj.oblocks(i)};
                    temp = flip(temp,2);
                    %size(temp)
                    if obj.oblocks_flip(i) == 1
                        temp = flip(temp);
                    end
                    blfield = [blfield; temp(1:end-1,:)];
                end
                blfield = blfield(obj.iLE:obj.iTE,:);
                %size(blfield,1)
                %blfield = blfield([obj.iLE:-1:1 end:-1:obj.iTE],:);
                %size(blfield,1)
            end
        end

        function value = blNormGrad(obj,prop)
            q = obj.oGridProp(prop);
            value = (q(:,2)-q(:,1))./(obj.yBL(:,2)-obj.yBL(:,1));
            value = [value (q(:,3:end)-q(:,1:end-2))./(obj.yBL(:,3:end)-obj.yBL(:,1:end-2))];
            value = [value (q(:,end)-q(:,end-1))./(obj.yBL(:,end)-obj.yBL(:,end-1))];

        end

        function plot_BL_profile(obj,x,prop,ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            if string(prop) == "dsdy"
                plot(ax, q, obj.yBL(i,2:end))
            else
                plot(ax, q, obj.yBL(i,:)/obj.yBL(i,end))
            end
        end

        

        function ind = x2ind(obj,x)
            [~, ind] = min(abs(obj.xSurf-x));
        end

        function getSize(obj)
            props = properties(obj);
            totSize = 0; 
            for ii=1:length(props) 
                currentProperty = obj.(props{ii});
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end

        function kPlot(obj,prop,ax,lims)
            
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis([-0.2 2 -0.5 0.5])
            axis equal
            axis off
            cb = colorbar('southoutside');
            if nargin == 5
                caxis(lims)
            end
        end

        function BLkPlot(obj,prop,ax,lims)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis([-0.1 1.1 0 0.15])
            axis equal
            cb = colorbar;
            if nargin == 5
                caxis(lims)
            end
        end

        function [xsurf, ysurf, n] = getSurfCoordsNorms(obj)
            xo = [];
            yo = [];
            for i=1:length(obj.oblocks)
                xtmp = obj.blk.x{obj.oblocks(i)}(:,:);
                ytmp = obj.blk.y{obj.oblocks(i)}(:,:);
                xtmp = flip(xtmp,2);
                ytmp = flip(ytmp,2);
                if obj.oblocks_flip(i) == 1
                    xtmp = flip(xtmp);
                    ytmp = flip(ytmp);
                end
                xo = [xo; xtmp(1:end-1,:)];
                yo = [yo; ytmp(1:end-1,:)];
            end
            xsurf = xo(:,1);
            ysurf = yo(:,1);

            %obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
            %size(obj.xSurf)
            R = [0 -1; 1 0];
            for i = 1:size(xo,1)
                
                if i==1
                    s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(end,1)); (yo(i,1)-yo(end,1))];
                elseif i==size(xo,1)
                    s1 = [(xo(1,1)-xo(i,1)); (yo(1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                else
                    s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                end
                n1 = R*s1/norm(s1);
                n2 = R*s2/norm(s2);
                nnow = (0.5*(n1+n2)/norm(0.5*(n1+n2)));
                n(:,i) = nnow;
            end
        end
    end
end
