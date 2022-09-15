
classdef DNS_case < handle
    %  DNS_CASE Class containing 3DNS case information and methods to read
    %  and postprocess results

    properties
        casename;
        casepath;
        runpath;
        runpaths;
        run;
        blk;
        next_block;
        next_patch;
        corner;
        solver;
        bcs;
        gas;
        x;
        y;
        nSlices;
        nProbes;
        nSkip;
        probes;
        NB;
        y_inlet;
        inlet_width;
        nj_inlet;
        meanFlow = meanSlice.empty;
        meanFlows = {};
        instFlow = volFlow.empty;
        iSlices;
        jSlices = jSlice.empty;
        kSlices = kSlice.empty;
        inflowTurb = volTurbulence.empty;
        RANSSlices;
        trip;
        iTrip = false;
        e_unst = [];
    end

    properties (Dependent = true)
        ftrip;
        Re_k        % Trip Reynolds number
    end

    methods
        function obj = DNS_case(casename,run)
            %DNS_CASE Construct an instance of this class
            %   Detailed explanation goes here
            obj.casename = casename;
            obj.casepath = fullfile(pwd,obj.casename);
            
            if nargin < 2 || isempty(run)
                obj.runpath = obj.casepath;
            elseif length(run) == 1
                obj.run = run;
                obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run)]);
            else
                obj.run = run;
                obj.runpaths = {};
                obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run(end))]); 
                for ir = 1:length(run)
                    obj.runpaths{ir} = fullfile(obj.casepath,['run' num2str(obj.run(ir))]);
                end
            end
            disp(obj.casepath)

            rcase = read_case(casename);
            obj.NB = rcase.NB;
            obj.blk = rcase.blk;
            obj.next_block = rcase.next_block;
            obj.next_patch = rcase.next_patch;
            obj.corner = rcase.corner;
            obj.bcs = rcase.bcs;
            obj.gas = rcase.gas;
            obj.solver = rcase.solver;
            obj.blk.inlet_blocks{1} = rcase.inlet_blocks;
            obj.blk.z = linspace(0, obj.solver.span, obj.blk.nk{1});
            
            if obj.NB == 12
                obj.blk.oblocks = [4 6 9 5];
                obj.blk.oblocks_flip = [0 0 1 1];
            end

            if isfile(fullfile(obj.runpath,'slice_time.txt'))
                obj.nSlices = size(readmatrix(fullfile(obj.runpath,'slice_time.txt')),1);
            end
%             cpu_file_path = fullfile(obj.runpath, 'input_cpu.txt');
%             fid = fopen(cpu_file_path);
%             obj.NB = str2num(fgetl(fid));
%             fprintf('NB = %d\n', obj.NB);
%             obj.blk.blockdims = zeros(obj.NB,3);
%             for nb=1:obj.NB
%                 nijk = str2num(fgetl(fid))
%                 obj.blk.blockdims(nb,:) = nijk;
%                 procdims = str2num(fgetl(fid))
%                 obj.solver.npp = floor(nijk(1)/procdims(1));
%                 if nb==1, obj.blk.npp = nijk(1)/procdims(1); end
%                 nskip = sum(str2num(fgetl(fid)) == 0);
%                 for i=1:nskip
%                     fgetl(fid);
%                 end
%             end
%             ncorner = str2num(fgetl(fid));
%             for i=1:ncorner
%                 temp = str2num(fgetl(fid));
%                 skip = temp(1);
%                 for i=1:skip
%                     fgetl(fid);
%                 end
%             end
%             temp = str2num(fgetl(fid));
%             obj.solver.niter = temp(1);
%             obj.solver.nwrite = temp(2);
%             obj.solver.ncut = temp(3);
% 
%             temp = str2num(fgetl(fid));
%             obj.solver.cfl = temp(1);
%             obj.solver.sigma = temp(2);
%             obj.solver.ifsplit = temp(3);
%             obj.solver.ifsat = temp(4);
%             obj.solver.ifLES = temp(5);
% 
%             temp = str2num(fgetl(fid));
%             obj.bcs.Toin = temp(1);
%             obj.bcs.Poin = temp(2);
%             obj.bcs.pexit = temp(3);
%             obj.bcs.vin = temp(4);
%             obj.bcs.alpha = temp(5);
%             obj.bcs.cax = temp(6);
%             obj.bcs.aturb = temp(7);
%             obj.bcs.lturb = temp(8);
%             obj.bcs.ilength = temp(9);
%             obj.bcs.radprof = temp(10);
%             obj.bcs.gamma = 0.0;
%             obj.bcs.g_z = 0.0;
% 
%             temp = str2num(fgetl(fid));
%             obj.gas.gam = temp(1);
%             obj.gas.cp = temp(2);
%             obj.gas.mu_ref = temp(3);
%             obj.gas.mu_tref = temp(4);
%             obj.gas.mu_cref = temp(5);
%             obj.gas.pr = temp(6);
% 
%             temp = str2num(fgetl(fid));
%             obj.solver.span = temp(1);
%             obj.solver.fexpan = temp(2);
%             obj.blk.span = temp(1);
% 
%             temp = str2num(fgetl(fid));
%             obj.solver.irestart = temp(1);
%             obj.solver.istats = temp(2);
%             
%             n_inlets = str2num(fgetl(fid));
%             for nin=1:n_inlets
%                 obj.blk.inlet_blocks{nin} = [];
%                 n_inlet_blocks = str2num(fgetl(fid));
%                 for i=1:n_inlet_blocks
%                     obj.blk.inlet_blocks{nin}(end+1) = str2num(fgetl(fid));
%                 end
%             end
%             temp = fgetl(fid);
%             obj.solver.istability = str2num(temp(1));
% 
%             if obj.solver.istability == 2
%                 obj.iTrip = true;
%                 %obj.setTrip()
%             end
% 
%             temp = str2num(fgetl(fid));
%             obj.bcs.twall = temp(3);
% 
%             fclose(fid);

%             [blk_tmp, ~] = read_grid(casename);
%             
%             for i=1:obj.NB
%                 obj.blk.x{i} = blk_tmp.x{i};
%                 obj.blk.y{i} = blk_tmp.y{i};
%                 obj.blk.nk{i} = blk_tmp.nk{i};
%             end

            obj.solver.nk = obj.blk.nk{1};

            if length(obj.blk.inlet_blocks) == 1
                obj.nj_inlet = 0;
                y_inlet = [];
                for i=1:length(obj.blk.inlet_blocks{1})
                    y_inlet = [y_inlet obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1)];
                    ymidnow(i) = obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,ceil(end/2));
                end
                obj.y_inlet = sort(unique(y_inlet));
                obj.nj_inlet = length(obj.y_inlet);
                obj.inlet_width = max(obj.y_inlet) - min(obj.y_inlet);
                
                [~,inds] = sort(ymidnow);
                obj.blk.inlet_blocks{1} = obj.blk.inlet_blocks{1}(inds);
            end
            

            obj.trip = [];
            if obj.solver.istability == 2
                obj.iTrip = true;
                obj.readTrip;
            else
                obj.iTrip = false;
            end
            slices = dir(fullfile(obj.runpath,'k_cuts','kcu2_1_*'));
            obj.nSlices = length(slices);
            if obj.NB == 12
                obj.blk.io_surfaces.blks = [1 2 2 7 11 11 10 12 12 8 3 3];
                obj.blk.io_surfaces.types = [1 1 3 3 3 2 2 2 4 4 4 1];
            end

        end

        function slice = readSingleKSlice(obj,numslice)
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =kSlice(obj.runpath,slicenum,obj.blk,obj.gas);
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
                    fprintf('Reading slice %d/%d\n',[i obj.nSlices])
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.runpath,slicenum,obj.blk,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("numslices",'var') == 1
                for i=numslices
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.runpath,slicenum,obj.blk,obj.gas);
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

        function readJSlices(obP    , runs, numslices)
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

        function writeInstFlow(obj, path)
            %WRITEINSTFLOW Write out intantaneous 3D flow to path

            if nargin < 2
                path = obj.casepath;
            end

            obj.instFlow.writeFlow(path);
        end

        function readMeanFlow(obj)
            %READMEANFLOW Read in 2D mean flow
            obj.meanFlow = meanSlice(obj.runpath,obj.blk,obj.gas);
        end

        function readMeanFlows(obj)
            %READMEANFLOWS Read in and average 2D mean flows for multiple
            %runs
            regions = obj.getIntRegions;
            for ir = 1:length(obj.run)
                fprintf('Reading meanSlice %d/%d (run %d)\n', [ir length(obj.run) obj.run(ir)])
                mF = meanSlice(obj.runpaths{ir},obj.blk,obj.gas);
                if ir == 1
                    obj.meanFlow = mF;
                    [e_unst_s, ~, ~, ~, ~] = entropy_balance(mF, regions);
                else
                    obj.meanFlow.addSlice(mF)
                end
                if ir == length(obj.run)
                    [e_unst_e, ~, ~, ~, ~] = entropy_balance(mF, regions);
                end
                clear mF
            end
            obj.e_unst = e_unst_e - e_unst_s;
        end

        function readInflowTurb(obj)
            %READINFLOWTURB Read inflow turbulence file
            obj.inflowTurb = volTurbulence(obj.casepath,obj.bcs.ilength,obj.blk.nk{1},obj.bcs.lturb,obj.y_inlet,obj.solver.span);
        end

        function readProbes(obj)
            f = fopen(fullfile(obj.runpath,'probe.txt'), 'r');
            temp = str2num(char(split(fgetl(f))));
            obj.nProbes = temp(1);
            obj.nSkip = temp(2);
            for nProbe = 1:obj.nProbes
                temp = str2num(char(split(fgetl(f))));
                blkData.nb = temp(1);
                blkData.i = temp(2);
                blkData.j = temp(3);
                blkData.k = temp(4);
                blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
                blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
                obj.probes{nProbe} = dnsProbe(obj.runpath, nProbe, obj.nSkip, blkData, obj.gas);
            end
            fclose(f);
        end

        function plot_probes(obj, ax)
            if nargin<2 || isempty(ax)
                ax = gca;
            end
            hold on
            C = colororder;
            sz = 25;
            for ip = 1:obj.nProbes
                scatter(ax, obj.probes{ip}.x, obj.probes{ip}.y, sz, C(mod(ip-1,7)+1,:), 'filled')
            end
        end

        function kPlot(obj,slice,prop,ax,lims,label)
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            %slice
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
            if string(prop) == "schlieren"
                colormap(gray)
                map = colormap;
                map = flip(map,1);
                colormap(map);
                if nargin < 6
                    label = '$|\nabla \rho|/\rho$';
                end
            end
            cb = colorbar;
            if nargin > 4 && ~isempty(lims)
                caxis(lims)
            end
            if nargin > 5 || exist("label",'var')
                label;
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end
            
            set(ax, 'FontSize', 12)
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
            pbaspect([8 1.5 1])

            axis([0.15 0.95 0 0.15])
            cb = colorbar('southoutside');
            if nargin > 4 && ~isempty(lims)
                caxis(lims)
            end
            nargin
            if nargin > 5 && ~isempty(label)
                label;
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end

            set(ax, 'FontSize', 12)
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
            if string(prop) == "cf" && isempty(slice.pdyn)
                    obj.setPdyn(slice);
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
                qtemp = q{obj.blk.oblocks(i)}(:,end);
                if obj.blk.oblocks_flip(i) == 1
                    temp = flip(temp);
                    qtemp = flip(qtemp);
                end
                xsurf = [xsurf temp'];
                qsurf = [qsurf]
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
                temp = fgetl(fid);
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

%         function [e_s, e_phi, e_irrev, e_N] = entropy_balance(obj)
% 
%             e_s = 0.0;
%             for ii=1:length(obj.blk.io_surfaces.blks)
%                 
%                 prop(:,:,1) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.u{obj.blk.io_surfaces.blks(ii)};
% 
%                 prop(:,:,2) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.v{obj.blk.io_surfaces.blks(ii)};
%                 
%                 temp(ii) = obj.surface_integral(prop, obj.blk.io_surfaces.blks(ii), ...
%                     obj.blk.io_surfaces.types(ii));
% 
%                 e_s = e_s+temp(ii);
% 
%                 if obj.blk.io_surfaces.blks(ii) == 3
%                     temp(ii)
%                 end
% 
%                 clear prop
% 
%             end
% 
%             e_phi = 0.0;
%             e_irrev = 0.0;
%             
%             parfor ib = 1:obj.NB
%                 ib
%                 prop_phi = obj.meanFlow.diss_T{ib};
%                 prop_phi = obj.meanFlow.diss{ib}./obj.meanFlow.T{ib};
%                 prop_irrev = obj.meanFlow.irrev_gen{ib};
%                 e_phi = e_phi + obj.area_integral(prop_phi, ib);
%                 e_irrev = e_irrev + obj.area_integral(prop_irrev, ib);
%             end
% 
%             e_N = e_s - e_phi - e_irrev;
% 
%             
%         end

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

%         function [e_s, e_phi, e_irrev, e_rev] = entropy_balance(obj)
%             
%             pool = gcp('nocreate');
%             if isempty(pool)
%                 parpool;
%             end
%             
%             %
%             
%             celldims = obj.blk.blockdims(:,1:2) - 1;
%             cells.x = zeros([sum(prod(celldims,2)) 1]);
%             cells.y = zeros([sum(prod(celldims,2)) 1]);
%             cells.area = zeros([sum(prod(celldims,2)) 1]);
%             cells.unst = zeros([sum(prod(celldims,2)) 1]);
%             cells.conv = zeros([sum(prod(celldims,2)) 1]);
%             cells.diss = zeros([sum(prod(celldims,2)) 1]);
%             cells.irrev = zeros([sum(prod(celldims,2)) 1]);
%             cells.rev = zeros([sum(prod(celldims,2)) 1]);
%             
%             
%             %
%             
%             parfor ib=1:obj.NB
%                 nentries = celldims(ib,1)*celldims(ib,2);
%             
%             
%                 [drousdx,~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rous{ib});
%                 [~,drovsdy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rovs{ib});
%                 [dqx_Tdx,~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_x{ib});
%                 [~,dqy_Tdy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_y{ib});
%             
%                 unst_prop = obj.meanFlow.ro{ib}.*dsdt;
%             
%                 conv_prop = obj.meanFlow.ro{ib}.*(obj.meanFlow.u{ib}.*dsdx + ...
%                    obj.meanFlow.v{ib}.*dsdy);
% 
%                 conv_prop = drousdx + drovsdy;
%             
%                 diss_prop = obj.meanFlow.diss_T{ib};
%                 irrev_prop = obj.meanFlow.irrev_gen{ib};
%                 rev_prop = dqx_Tdx + dqy_Tdy;
%             
%                 xtmp{ib} = zeros(nentries, 1);
%                 ytmp{ib} = zeros(nentries, 1);
%                 areatmp{ib} = zeros(nentries, 1);
%                 unsttmp{ib} = zeros(nentries,1);
%                 convtmp{ib} = zeros(nentries, 1);
%                 disstmp{ib} = zeros(nentries, 1);
%                 irrevtmp{ib} = zeros(nentries, 1);
%                 revtmp{ib} = zeros(nentries, 1);
%                 
%                 for i=1:celldims(ib,1)
%                     if mod(i, 50) == 0
%                          sprintf('Block %d, i=%d/%d', ib, i, celldims(ib,1))
%                     end
%                     for j=1:celldims(ib,2)
%                         pos = celldims(ib,2)*(i-1)+j;
%                         xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
%                             obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
%                         ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
%                             obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];
%             
%                         xtmp{ib}(pos) = mean(xnow);
%                         ytmp{ib}(pos) = mean(ynow);
%             
%                         area = polyarea(xnow,ynow);
%                         areatmp{ib}(pos) = area;
%             
%                         unsttmp{ib}(pos) = 0.25*(unst_prop(i,j)+unst_prop(i+1,j)+unst_prop(i+1,j+1)+unst_prop(i,j+1))*area;
%                         convtmp{ib}(pos) = 0.25*(conv_prop(i,j)+conv_prop(i+1,j)+conv_prop(i+1,j+1)+conv_prop(i,j+1))*area;
%                         disstmp{ib}(pos) = 0.25*(diss_prop(i,j)+diss_prop(i+1,j)+diss_prop(i+1,j+1)+diss_prop(i,j+1))*area;
%                         irrevtmp{ib}(pos) = 0.25*(irrev_prop(i,j)+irrev_prop(i+1,j)+irrev_prop(i+1,j+1)+irrev_prop(i,j+1))*area;
%                         revtmp{ib}(pos) = 0.25*(rev_prop(i,j)+rev_prop(i+1,j)+rev_prop(i+1,j+1)+rev_prop(i,j+1))*area;
%             
%                     end
%                 end
%             end
%             
%             for ib = 1:obj.NB
%             
%                 offset = sum(prod(celldims(1:ib,:),2));
%                 nentries = celldims(ib,1)*celldims(ib,2);
%             
%                 cells.x(offset-nentries+1:offset) = xtmp{ib};
%                 cells.y(offset-nentries+1:offset) = ytmp{ib};
%                 cells.area(offset-nentries+1:offset) = areatmp{ib};
%                 %cells.unst(offset-nentries+1:offset) = unsttmp{ib};
%                 cells.conv(offset-nentries+1:offset) = convtmp{ib};
%                 cells.diss(offset-nentries+1:offset) = disstmp{ib};
%                 cells.irrev(offset-nentries+1:offset) = irrevtmp{ib};
%                 cells.rev(offset-nentries+1:offset) = revtmp{ib};
%             end
%             
%             %%
%             
%             indices = cells.x > 0.75 & cells.x < 1.2 & cells.y > -0.1 & cells.y < 0.1;
%             
%             %e_unst = sum(cells.unst(indices));
%             e_s = sum(cells.conv(indices));
%             e_phi = sum(cells.diss(indices));
%             e_irrev = sum(cells.irrev(indices));
%             e_rev = sum(cells.rev(indices));
% 
%         end

        function plot_mesh(obj, skip, ax)
            if nargin<2 || isempty(skip)
                skip=8;
            end

            if nargin<3 || isempty(ax)
                ax = gca;
            end
            
            hold on
            for ib=1:obj.NB
                [ni, nj] = size(obj.blk.x{ib});
                for j=[1:skip:ni ni]
                    if (j==1) || (j==ni)
                        plot(ax, obj.blk.x{ib}(j,:),obj.blk.y{ib}(j,:),'r','LineWidth',1)
                    else
                        plot(ax, obj.blk.x{ib}(j,:),obj.blk.y{ib}(j,:),'k')
                    end
                end
                for j=[1:skip:nj nj]
                    if (j==1) || (j==nj)
                        plot(ax, obj.blk.x{ib}(:,j),obj.blk.y{ib}(:,j),'r','LineWidth',1)
                    else
                        plot(ax, obj.blk.x{ib}(:,j),obj.blk.y{ib}(:,j),'k')
                    end
                end
            end
            
            axis equal
        end

        function setPdyn(obj,slice)
            if isempty(obj.meanFlow)
                obj.readMeanFlow;
            end
            slice.pdyn = obj.meanFlow.pdyn;
        end

        function mesh_analysis(obj, skip, ax)

            if nargin<2 || isempty(skip), skip=8; end

            mesh_analysis(obj.blk, skip)

%             if nargin < 3 || isempty(ax)
%                 ax = gca;
%             end
%             C = colororder;
%             
%             npp = obj.blk.npp;
%             
%             nij_pts = 0;
%             hold on
%             for i=1:obj.NB
%                 x = obj.blk.x{i};
%                 y = obj.blk.y{i};
%                 x = x(:,:,1);
%                 y = y(:,:,1);
%                 ni = size(x,1);
%                 nj = size(x,2);
%                 nij_pts = nij_pts + ni*nj;
%                 %nk = blk{i}.nk;
%                 for j=[1:skip:ni ni]
%                     if (j==1) || (j==ni)
%                         plot(x(j,:),y(j,:),'r','LineWidth',1)
%                     else
%                         plot(x(j,:),y(j,:),'k')
%                     end
%                 end
%                 for j=[1:skip:nj nj]
%                     if (j==1) || (j==nj)
%                         plot(x(:,j),y(:,j),'r','LineWidth',1)
%                     else
%                         plot(x(:,j),y(:,j),'k')
%                     end
%                 end
%             end

        end

        function update_input_files(obj)
            if obj.NB == 12
                topology = 2;
            else
                topology = 1;
            end
            write_input_files(obj.casename,obj.blk,obj.next_block,obj.next_patch,obj.corner,obj.bcs,obj.gas,obj.solver,topology);
        end

        function [e_unst, e_s, e_phi, e_irrev, e_rev] = getEntropyTerms(obj, mF)

            if nargin < 2
                mF = obj.meanFlow;
            end

            regions = obj.getIntRegions;

%             if length(obj.run) == 1
            if isempty(mF.e_s)
                [~, e_s, e_phi, e_irrev, e_rev] = entropy_balance(mF, regions);
                mF.e_s = e_s;
                mF.e_phi = e_phi;
                mF.e_irrev = e_irrev;
                mF.e_rev = e_rev;
            else
                e_s = mF.e_s;
                e_phi = mF.e_phi;
                e_irrev = mF.e_irrev;
                e_rev = mF.e_rev;
            end
            if isempty(obj.e_unst)
                e_unst = zeros(size(e_s));
            else
                e_unst = obj.e_unst;
            end
%                e_unst = zeros(size(e_unst));
%             else
%                 es_unst = {};
%                 es_s = {};
%                 es_phi = {};
%                 es_irrev = {};
%                 es_rev = {};
%                 for ir = 1:length(obj.run)
%                     [es_unst{ir}, es_s{ir}, es_phi{ir}, es_irrev{ir}, es_rev{ir}] = entropy_balance(mF, regions);
%                 end
%                 for i = 1:length(regions)
%                     e_unst(i) = es_unst{end}(i)-es_unst{1}(i);
%                     e_s(i) = 0; e_phi(i) = 0; e_irrev(i) = 0; e_rev(i) = 0;
%                     for j = 1:length(obj.run)
%                         e_s(i) = e_s(i) + es_s{j}(i);
%                         e_phi(i) = e_phi(i) + es_phi{j}(i);
%                         e_irrev(i) = e_irrev(i) + es_irrev{j}(i);
%                         e_rev(i) = e_rev(i) + es_rev{j}(i);
%                     end
%                     e_s(i) = e_s(i)/length(obj.run);
%                     e_phi(i) = e_phi(i)/length(obj.run);
%                     e_irrev(i) = e_irrev(i)/length(obj.run);
%                     e_rev(i) = e_rev(i)/length(obj.run);
%                 end
%             end
        end

        function [e_Pr, e_diss] = getTkeTerms(obj, mF)

            if nargin < 2
                mF = obj.meanFlow;
            end

            regions = obj.getIntRegions;

%             if length(obj.run) == 1
            if isempty(mF.e_Pr)
                [e_Pr, e_diss] = tke_balance(mF, regions);
                mF.e_Pr = e_Pr;
                mF.e_diss = e_diss;
            else
                e_Pr = mF.e_Pr;
                e_diss = mF.e_diss;
            end

        end

        function plot_diss_regions(obj)
            [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            enow = e_s/e_s(1);
            X = categorical({'Pre-shock', 'Post-shock', 'Trailing edge', 'Wake'});
            X = reordercats(X,{'Pre-shock', 'Post-shock', 'Trailing edge', 'Wake'});
            bar(X,enow)
            ylabel('\epsilon_S/\epsilon_{S, Pre-shock}')
            set(gca, 'FontSize', 12)
            set(gca, 'XTickLabelRotation', 20)
        end

        function plot_entropy_stackup(obj, normalise)
            if nargin < 2
                normalise = false;
            end

            [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            if ~isempty(obj.e_unst)
                e_unst = obj.e_unst
            end
            normalise

            if normalise
                for i=1:length(e_s)
                    e_unst(i) = e_unst(i)/abs(e_s(i));
                    e_phi(i) = e_phi(i)/abs(e_s(i));
                    e_irrev(i) = e_irrev(i)/abs(e_s(i));
                    e_rev(i) = e_rev(i)/abs(e_s(i));
                    e_s(i) = e_s(i)/abs(e_s(i));
                end
            else
                factor = abs(e_phi(1));
                for i=1:length(e_s)
                    e_unst(i) = e_unst(i)/factor;
                    e_phi(i) = e_phi(i)/factor;
                    e_irrev(i) = e_irrev(i)/factor;
                    e_rev(i) = e_rev(i)/factor;
                    e_s(i) = e_s(i)/factor;
                end
            end

            for ng = 1:length(e_s)
                stackData(ng, 1, 1) = e_s(ng);
                stackData(ng, 1, 2) = e_unst(ng);
                stackData(ng, 1, 3:5) = 0;
            
                stackData(ng, 2, 1:2) = 0;
                stackData(ng, 2, 3) = e_phi(ng);
                stackData(ng, 2, 4) = e_irrev(ng);
                stackData(ng, 2, 5) = e_rev(ng);
            end
            
            groupLabels = {"Pre-shock", "Post-shock", "Trailing edge", "Wake"};
            h = plotBarStackGroups(stackData, groupLabels);
            c = colororder;
            c = c(1:5,:);
            c = repelem(c,size(h,1),1); 
            c = mat2cell(c,ones(size(c,1),1),3);
            set(h,{'FaceColor'},c);
            legend(h(1,:),'\epsilon_S','\epsilon_{unst}','\epsilon_\phi','\epsilon_{irrev}','\epsilon_{rev}','Location','northeast')
            set(gca,'FontSize',12)
            if normalise
                ylabel('\epsilon/|\epsilon_\Phi|')
            else
                ylabel('\epsilon/|\epsilon_{\Phi, Pre-Shock}|')
            end
            set(gca, 'XTickLabelRotation',20)
        end
        

        function value = get.Re_k(obj)
            if isempty(obj.trip)
                obj.setTrip;
            end
            tripInd = obj.meanFlow.x2ind(obj.trip.x);
            inds = obj.meanFlow.BLedgeInd;
            nunow = obj.meanFlow.nu_e;
            Unow = obj.meanFlow.U;
            Ue = Unow(tripInd,inds(tripInd));
            nu_k = nunow(tripInd);
            value = Ue*obj.trip.scale/nu_k;
        end

        function regions = getIntRegions(obj)
            xrangePreShock = [0.16 0.59];
            xrangePostShock = [0.70 0.98];
            xMaxWake = 1.4;

            regions = {};
            regions{1}.nb = 6;
            [~, regions{1}.is] = min(abs(obj.blk.x{6}(:,end) - xrangePreShock(1)));
            [~, regions{1}.ie] = min(abs(obj.blk.x{6}(:,end) - xrangePreShock(2)));
            regions{1}.js = 1;
            regions{1}.je = obj.blk.blockdims(6,2);
            
            regions{2}.nb = 6;
            [~, regions{2}.is] = min(abs(obj.blk.x{6}(:,end) - xrangePostShock(1)));
            [~, regions{2}.ie] = min(abs(obj.blk.x{6}(:,end) - xrangePostShock(2)));
            regions{2}.js = 1;
            regions{2}.je =  obj.blk.blockdims(6,2);
            
            regions{3}.nb = 9;
            regions{3}.is = 1; regions{3}.ie = size(obj.blk.x{9},1);
            regions{3}.js = 1; regions{3}.je = size(obj.blk.x{9},2);

            regions{4}.nb = 10;
            regions{4}.is = 1;
            [~, regions{4}.ie] = min(abs(obj.blk.x{10}(:,obj.blk.blockdims(10,2)/2) - xMaxWake));
            regions{4}.js = 1; regions{4}.je = obj.blk.blockdims(10,2);
        end

    end
end