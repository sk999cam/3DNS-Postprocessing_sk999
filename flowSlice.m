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
        oblocks;
        oblocks_flip;
        iLE;
        iTE;
        n;
        ssurf;          % Surface distance fron LE

    end

    properties (Dependent = true)
        T;              % Temperature
        p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        mu;             % Viscosity
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


        function value = get.p(obj)
            disp('Calculating p')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            disp('Calculating T')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                value{nb} = pnow./(obj.ro{nb}*obj.gas.rgas);
            end
        end

        function value = get.vel(obj)
            disp('Calculating vel')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.M(obj)
            disp('Calculating M')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
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
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.mu_ref*(Tnow/obj.gas.mu_tref).^(3/2) .* (obj.gas.mu_cref + obj.gas.mu_tref)./(obj.gas.mu_cref + Tnow);
            end
        end


        function [profile, i] = BLprof(obj, x, prop)
            %BLPROF get a profile of prop across the BL at specified x
            [~, i] = min(abs(obj.xSurf-x));
            if any(strcmp(["dsdy" "U"],prop))
                profile = obj.(prop)(i,:);
            else
                propfield = oGridProp(obj,prop);
                profile = propfield(i,:);
            end
        end

        function blfield = oGridProp(obj, prop)
            %OGRIDPROP Construct array of property in o grid for one surface
            %of blade
            
            if any(strcmp(["U"], prop))
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
            cb = colorbar;
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
    end
end