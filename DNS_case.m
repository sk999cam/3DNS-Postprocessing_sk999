classdef DNS_case < handle
    %  DNS_CASE Class containing 3DNS case information and methods to read
    %  and postprocess results

    properties
        casename;
        casepath;
        runpath;
        run;
        blk;
        solver;
        bcs;
        gas;
        x;
        y;
        nSlices;
        NB;
        inlet_width;
        nj_inlet;
        meanFlow = meanSlice.empty;
        instFlow = volFlow.empty;
        iSlices;
        jSlices = jSlice.empty;
        kSlices = kSlice.empty;
        inflowTurb = volTurbulence.empty;
        RANSSlices;
        trip;
        iTrip = false;
    end

    properties (Dependent = true)
        ftrip;     
    end

    methods
        function obj = DNS_case(casename,run)
            %DNS_CASE Construct an instance of this class
            %   Detailed explanation goes here
            obj.casename = casename;
            obj.casepath = fullfile(pwd,obj.casename);
            
            if nargin > 1
                obj.run = run;
                obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run)]);
            else
                obj.runpath = obj.casepath;
            end
            disp(obj.casepath)

            obj.blk.oblocks = [4 6 9 5];
            obj.blk.oblocks_flip = [0 0 1 1];
            if isfile(fullfile(obj.runpath,'slice_time.txt'))
                obj.nSlices = size(readmatrix(fullfile(obj.runpath,'slice_time.txt')),1);
            end
            cpu_file_path = fullfile(obj.runpath, 'input_cpu.txt');
            fid = fopen(cpu_file_path);
            obj.NB = str2num(fgetl(fid));
            fprintf('NB = %d\n', obj.NB);
            obj.blk.blockdims = zeros(obj.NB,3);
            for nb=1:obj.NB
                nijk = str2num(fgetl(fid));
                obj.blk.blockdims(nb,:) = nijk;
                fgetl(fid);
                nskip = sum(str2num(fgetl(fid)) == 0);
                for i=1:nskip
                    fgetl(fid);
                end
            end
            ncorner = str2num(fgetl(fid));
            for i=1:ncorner
                temp = str2num(fgetl(fid));
                skip = temp(1);
                for i=1:skip
                    fgetl(fid);
                end
            end
            temp = str2num(fgetl(fid));
            obj.solver.niter = temp(1);
            obj.solver.nwrite = temp(2);
            obj.solver.ncut = temp(3);

            temp = str2num(fgetl(fid));
            obj.solver.cfl = temp(1);
            obj.solver.sigma = temp(2);
            obj.solver.ifsplit = temp(3);
            obj.solver.ifsat = temp(4);
            obj.solver.ifLES = temp(5);

            temp = str2num(fgetl(fid));
            obj.bcs.Toin = temp(1);
            obj.bcs.Poin = temp(2);
            obj.bcs.pexit = temp(3);
            obj.bcs.vin = temp(4);
            obj.bcs.alpha = temp(5);
            obj.bcs.cax = temp(6);
            obj.bcs.aturb = temp(7);
            obj.bcs.lturb = temp(8);
            obj.bcs.ilength = temp(9);
            obj.bcs.radprof = temp(10);
            obj.bcs.gamma = 0.0;
            obj.bcs.g_z = 0.0;

            temp = str2num(fgetl(fid));
            obj.gas.gam = temp(1);
            obj.gas.cp = temp(2);
            obj.gas.mu_ref = temp(3);
            obj.gas.mu_tref = temp(4);
            obj.gas.mu_cref = temp(5);
            obj.gas.pr = temp(6);

            temp = str2num(fgetl(fid));
            obj.solver.span = temp(1);
            obj.solver.fexpan = temp(2);
            obj.blk.span = temp(1);

            temp = str2num(fgetl(fid));
            obj.solver.irestart = temp(1);
            obj.solver.istats = temp(2);
            
            n_inlets = str2num(fgetl(fid));
            for nin=1:n_inlets
                obj.blk.inlet_blocks{nin} = [];
                n_inlet_blocks = str2num(fgetl(fid));
                for i=1:n_inlet_blocks
                    obj.blk.inlet_blocks{nin}(end+1) = str2num(fgetl(fid));
                end
            end
            fgetl(fid);
            temp = str2num(fgetl(fid));
            obj.bcs.twall = temp(3);

            fclose(fid);

            [blk_tmp, ~] = read_grid(casename);
            
            for i=1:obj.NB
                obj.blk.x{i} = blk_tmp{i}.x;
                obj.blk.y{i} = blk_tmp{i}.y;
                obj.blk.nk{i} = blk_tmp{i}.nk;
            end

            if length(obj.blk.inlet_blocks) == 1
                obj.nj_inlet = 0;
                for i=1:length(obj.blk.inlet_blocks{1})
                    max_ys(i) = max(obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1));
                    min_ys(i) = min(obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1));
                    obj.nj_inlet = obj.nj_inlet + obj.blk.blockdims(obj.blk.inlet_blocks{1}(i),2);
                end
                obj.nj_inlet = obj.nj_inlet - (length(obj.blk.inlet_blocks{1})-1);
                obj.inlet_width = max(max_ys) - min(min_ys);
            end
            obj.trip = [];
            obj.iTrip = false;
            slices = dir(fullfile(obj.runpath,'k_cuts','kcu2_1_*'));
            obj.nSlices = length(slices);

            obj.blk.io_surfaces.blks = [1 2 2 7 11 11 10 12 12 8 3 3];
            obj.blk.io_surfaces.types = [1 1 3 3 3 2 2 2 4 4 4 1];

        end

        function slice = readSingleKSlice(obj,numslice)
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =kSlice(obj.runpath,slicenum,obj.blk.blockdims,obj.gas);
        end

        function slice = readSingleJSlice(obj,numslice)
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
        end

        function readKSlices(obj, runs, numslices)
            %READKSLICES Read in instantaneous k slices
            % Optional: numslices - only read last n slices if there are many
            
            exist("runs",'var')
            exist("numslices",'var')
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slices = dir(fullfile(obj.runpath,'k_cuts','kcu2_1_*'));
            obj.nSlices = length(slices);
            if nargin == 1
                for i=1:obj.nSlices
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.runpath,slicenum,obj.blk.blockdims,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("numslices",'var') == 1
                for i=numslices
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.runpath,slicenum,obj.blk.blockdims,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("runs",'var') == 1
                runs
                obj.nSlices = 0
                slicenums=[];
                for nrun=runs
                    runpath = fullfile(obj.casepath,['run' num2str(nrun)]);
                    slices = dir(fullfile(runpath,'k_cuts','kcu2_1_*'));
                    obj.nSlices = obj.nSlices + length(slices);
                    read = 0;
                    for i=1:length(slices)
                        read = read+1;
                        slicenum = str2num(slices(i).name(8:end));
                        slicenums(end+1) = slicenum;
                        obj.kSlices(end+1) = kSlice(runpath,slicenum,obj.blk.blockdims,obj.gas);
                    end
                end
                [~,inds] = sort([obj.kSlices.nSlice]);
                obj.kSlices = obj.kSlices(inds);

            end
        end

        function readJSlices(obj, runs, numslices)
            %READJSLICES Read in instantaneous j slices
            % Optional: numslices - only read last n slices if there are many
            
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slices = dir(fullfile(obj.runpath,'j_cuts','jcu2_4_*'));
            %obj.nSlices = length(slices);
            obj.nSlices = length(slices);
            if nargin == 1
                for i=1:obj.nSlices
                    slicenum = str2num(slices(i).name(8:end));
                    obj.jSlices(i) = jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("numslices",'var') == 1
                for i=1:numslices
                    slicenum = str2num(slices(i+obj.nSlices-numslices).name(8:end));
                    obj.jSlices(i) = jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("runs",'var') == 1
                runs
                obj.nSlices = 0;
                slicenums=[];
                for nrun=runs
                    runpath = fullfile(obj.casepath,['run' num2str(nrun)]);
                    slices = dir(fullfile(runpath,'j_cuts','jcu2_4_*'));
                    obj.nSlices = obj.nSlices + length(slices);
                    read = 0;
                    for i=1:length(slices)
                        read = read+1;
                        slicenum = str2num(slices(i).name(8:end));
                        slicenums(end+1) = slicenum;
                        obj.jSlices(end+1) = jSlice(runpath,slicenum,obj.blk,obj.gas);
                    end
                end
                [~,inds] = sort([obj.jSlices.nSlice]);
                obj.jSlices = obj.jSlices(inds);

            end
        end

        function readInstFlow(obj)
            %READINSTFLOW Read in instantaneous 3D flow
            obj.instFlow = volFlow(obj.runpath,obj.blk,obj.gas);
            
        end

        function readMeanFlow(obj)
            %READMEANFLOW Read in 2D mean flow
            obj.meanFlow = meanSlice(obj.runpath,obj.blk,obj.gas);
        end

        function readInflowTurb(obj)
            %READINFLOWTURB Read inflow turbulence file
            obj.inflowTurb = volTurbulence(obj.nj_inlet,obj.blk.nk{1},obj.inlet_width,obj.casepath);
        end

        function kPlot(obj,slice,prop,ax,lims,label)
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            slice
            if isempty(slice)
                disp('here')
                q = obj.(prop);
            else
                q = slice.(prop);
            end
            hold on
            for i=1:slice.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis([-0.6 2 -0.5 0.5])
            axis equal
            cb = colorbar;
            if nargin > 4
                caxis(lims)
            end
            nargin
            if nargin > 5
                label
                cb.Label.String = label;
            end
        end

        function BLkPlot(obj,slice,prop,ax,lims,label)
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            if isempty(slice)
                disp('here')
                q = obj.(prop);
            else
                q = slice.(prop);
            end
            hold on
            for i=1:slice.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            
            axis equal
            pbaspect([0.2 0.15 1])

            axis([0.6 0.8 0 0.15])
            cb = colorbar;%('southoutside');
            if nargin > 4 && ~isempty(lims)
                caxis(lims)
            end
            nargin
            if nargin > 5 && ~isempty(label)
                label
                cb.Label.String = label;
            end
        end

        function jPlot(obj,slice,prop,ax,lims,label)
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            if nargin < 5 || isempty(lims)
                lims = [];
            end
            if nargin < 6 || isempty(label)
                label = [];
            end
            slice.jPlot(prop,ax,lims,label)
            
%             if isempty(slice)
%                 disp('here')
%                 q = obj.(prop);
%             else
%                 q = slice.(prop);
%             end
% 
%             pcolor(ax, slice.X, slice.Z, q);
%             
%             shading('interp')
%             %axis([-0.2 2 -0.5 0.5])
%             axis equal
%             cb = colorbar('southoutside');
%             if nargin > 4
%                 caxis(lims)
%             end
%             nargin
%             if nargin > 5
%                 label
%                 cb.Label.String = label;
%             end
        end

        function jkPlot(obj, slice, jProp, kProp, jLims, kLims)
            subplot(3,5,1:10)
            if nargin > 4
                obj.kPlot(obj.kSlices(slice),kProp,kLims)
                subplot(3,5,11:14)
                dx = abs(obj.jSlices(slice).X(end,1)-obj.jSlices(slice).X(1,1));
                dy = abs(obj.jSlices(slice).Z(1,end)-obj.jSlices(slice).Z(1,1));
                obj.jPlot(obj.jSlices(slice),jProp,jLims);
                pbaspect(gca, [dx dy dx]);
            else
                obj.kPlot(obj.kSlices(slice),kProp)
                subplot(3,5,11:14)
                dx = abs(obj.jSlices(slice).X(end,1)-obj.jSlices(slice).X(1,1));
                dy = abs(obj.jSlices(slice).Z(1,end)-obj.jSlices(slice).Z(1,1));
                obj.jPlot(obj.jSlices(slice),jProp);
                pbaspect(gca, [dx dy dx]);
            end
        end

        function plot_surf_prop(obj,slice,prop,ax,lims)

            if nargin < 4 || isempty(ax)
                ax = gca;
            end

            q = slice.(prop);
            xsurf = [];
            for i=1:length(obj.blk.oblocks)
                temp = obj.blk.x{obj.blk.oblocks(i)}(:,end);
                if obj.blk.oblocks_flip(i) == 1
                    temp = flip(temp);
                end
                xsurf = [xsurf temp'];
            end

            [~,iLE] = min(xsurf);
            [~,iTE] = max(xsurf);
            xsurf = xsurf(iLE:iTE);
            q = q(iLE:iTE);

            plot(ax,xsurf,q)

            if nargin == 5
                ylim(lims)
            end
        end

        function plot_BL_profile(obj,slice,x,prop,ax)
            
        end

        function setTrip(obj, x, amp, scale)
            obj.trip.amp = amp;
            obj.trip.scale = scale;
            [obj.trip.nb, obj.trip.i] = obj.x2point(x);
            obj.trip.x = obj.blk.x{obj.trip.nb}(obj.trip.i,end);
            obj.trip.y = obj.blk.y{obj.trip.nb}(obj.trip.i,end);
            obj.iTrip = true;
        end

        function [nb, ni] = x2point(obj, x)
            % Find closest point on top surface
            for i=1:length(obj.blk.oblocks)
                [mindist(i), ind(i)] = min(abs(obj.blk.x{obj.blk.oblocks(i)}(:,end).*double(obj.blk.y{obj.blk.oblocks(i)}(:,end)>0) - x));
            end
            [~,n] = min(mindist);
            nb = obj.blk.oblocks(n);
            ni = ind(n);
        end

        function value = get.ftrip(obj)
            if ~obj.iTrip
                disp('Setting trip with default params')
                obj.setTrip(0.1, 1, 100)
            end
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                dist2 = ((obj.blk.x{nb} - obj.trip.x).^2 + (obj.blk.y{nb} - obj.trip.y).^2);
                if obj.trip.mode == 1
                    f = obj.trip.amp*exp(-dist2*obj.trip.scale);
                    value{nb} = min(f,1.0);
                else
                    value{nb} = double(sqrt(dist2)<obj.trip.scale);
                end

            end
        end

        function clearTrip(obj)

            obj.trip.amp
            obj.trip.scale
            [obj.trip.nb, obj.trip.i]
            obj.trip.x
            obj.trip.y
            obj.iTrip = false;
        end

        function plotTripFactor(obj,ax)
            if ~obj.iTrip
                disp('Must setup trip first')
            else
                ynow = obj.meanFlow.yBL(obj.meanFlow.x2ind(obj.trip.x),:);
                f = obj.ftrip{obj.trip.nb}(obj.trip.i,end:-1:1);
                U = obj.meanFlow.U(obj.meanFlow.x2ind(obj.trip.x),:);
                if nargin > 1
                    plot(ax, f, ynow)
                    hold on
                    plot(ax, U/max(U), ynow)
                else
                    plot(f, ynow)
                    hold on
                    plot(U/max(U), ynow)
                end
            end
        end

        function readTrip(obj)
            filepath = fullfile(obj.casepath,'trip.txt');
            if ~exist(filepath)
                disp("Can't find trip.txt")
            else
                fid = fopen(filepath);
                obj.trip.mode = str2num(fgetl(fid));
                temp = fgetl(fid)
                temp = str2num(temp);
                if obj.trip.mode == 1
                    obj.trip.amp = temp(1);
                    obj.trip.scale = temp(2);
                    xtmp = temp(3);
                else
                    obj.trip.scale = temp(1);
                    xtmp = temp(2);
                end

                fclose(fid);
                [obj.trip.nb, obj.trip.i] = obj.x2point(xtmp);
                obj.trip.x = obj.blk.x{obj.trip.nb}(obj.trip.i,end);
                obj.trip.y = obj.blk.y{obj.trip.nb}(obj.trip.i,end);
                obj.iTrip = true;
            end
        end

        function writeTripInput(obj, dir)
            if ~obj.iTrip
                disp('Must setup trip first')
            else
                if nargin > 1
                    fid = fopen(fullfile(dir, 'trip.txt'), 'w');
                else
                    fid = fopen(fullfile(obj.casepath, 'trip.txt'), 'w');
                end
                fprintf(fid,'%20.16e %20.16e %20.16e %20.16e', obj.trip.amp, obj.trip.scale, obj.trip.x, obj.trip.y);
                fclose(fid);
            end
        end

        function readRANSSlice(obj,turb,mod)
            ransdir = fullfile('RANS','cwl90');
            obj.RANSSlices{turb.mod} = RANSSlice(ransdir,data,obj.blk,obj.gas);
        end

        function [e_s, e_phi, e_irrev, e_N] = entropy_balance(obj)

            e_s = 0.0;
            for ii=1:length(obj.blk.io_surfaces.blks)
                
                prop(:,:,1) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
                    obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
                    obj.meanFlow.u{obj.blk.io_surfaces.blks(ii)};

                prop(:,:,2) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
                    obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
                    obj.meanFlow.v{obj.blk.io_surfaces.blks(ii)};
                
                temp(ii) = obj.surface_integral(prop, obj.blk.io_surfaces.blks(ii), ...
                    obj.blk.io_surfaces.types(ii));

                e_s = e_s+temp(ii);

                if obj.blk.io_surfaces.blks(ii) == 3
                    temp(ii)
                end

                clear prop

            end

            e_phi = 0.0;
            e_irrev = 0.0;
            
            parfor ib = 1:obj.NB
                ib
                prop_phi = obj.meanFlow.diss_T{ib};
                prop_phi = obj.meanFlow.diss{ib}./obj.meanFlow.T{ib};
                prop_irrev = obj.meanFlow.irrev_gen{ib};
                e_phi = e_phi + obj.area_integral(prop_phi, ib);
                e_irrev = e_irrev + obj.area_integral(prop_irrev, ib);
            end

            e_N = e_s - e_phi - e_irrev;

            
        end

        function value = surface_integral(obj, prop, ib, type)
            if ib == 3
                ib
                type
            end
            switch type
                case 1
                    xnow = obj.blk.x{ib}(1,end:-1:1);
                    ynow = obj.blk.y{ib}(1,end:-1:1);
                    qx = prop(1,end:-1:1,1);
                    qy = prop(1,end:-1:1,2);
                case 2 
                    xnow = obj.blk.x{ib}(end,:);
                    ynow = obj.blk.y{ib}(end,:);
                    qx = prop(end,:,1);
                    qy = prop(end,:,2);
                case 3
                    xnow = obj.blk.x{ib}(:,1);
                    ynow = obj.blk.y{ib}(:,1);
                    qx = prop(:,1,1);
                    qy = prop(:,1,2);
                case 4
                    xnow = obj.blk.x{ib}(end:-1:1,end);
                    ynow = obj.blk.y{ib}(end:-1:1,end);
                    qx = prop(end:-1:1,end,1);
                    qy = prop(end:-1:1,end,2);
            end

            xnow = reshape(xnow,1,[]);
            ynow = reshape(ynow,1,[]);
            qx = reshape(qx,1,[]);
            qy = reshape(qy,1,[]);

            dx = xnow(2:end)-xnow(1:end-1);
            dy = ynow(2:end)-ynow(1:end-1);
            ds = sqrt(dx.^2+dy.^2);
            s = [0 cumsum(ds)];
            
            n(1,:) = [dy(1)/ds(1) -dx(1)/ds(1)];
            n(length(xnow),:) = [dy(end)/ds(end) -dx(end)/ds(end)];

            for i=2:length(xnow)-1
                n1 = [dy(i-1)/ds(i-1) -dx(i-1)/ds(i-1)];
                n2 = [dy(i)/ds(i) -dx(i)/ds(i)];
                n(i,:) = (n1+n2)/norm(n1+n2);
            end
            

            integrand = n(:,1).*qx' + n(:,2).*qy';
            value = trapz(s',integrand);


        end

        function value = area_integral(obj, prop, nb)
            value = 0.0;
            for i=1:obj.blk.blockdims(nb,1)-1
                for j=1:obj.blk.blockdims(nb,2)-1
                    xnow = [obj.blk.x{nb}(i,j) obj.blk.x{nb}(i+1,j) ...
                        obj.blk.x{nb}(i+1,j+1) obj.blk.x{nb}(i,j+1)];
                    ynow = [obj.blk.y{nb}(i,j) obj.blk.y{nb}(i+1,j) ...
                        obj.blk.y{nb}(i+1,j+1) obj.blk.y{nb}(i,j+1)];

                    q = 0.25*(prop(i,j)+prop(i+1,j)+prop(i+1,j+1)+prop(i,j+1));

                    area = polyarea(xnow,ynow);
                    value = value+area*q;
                end
            end
        end


        

    end
end