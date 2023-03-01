classdef DNS_case < handle
    %  DNS_CASE Class containing 3DNS case information and methods to read
    %  and postprocess results

    properties
        casename;
        casepath;
        casetype;
        topology;
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
        cell_area = [];
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
                obj.run = [];
            elseif length(run) == 1
                obj.run = run;
                obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run)]);
                obj.runpaths = obj.runpath;
            else
                obj.run = run;
                obj.runpaths = {};
                obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run(end))]); 
                for ir = 1:length(run)
                    obj.runpaths{ir} = fullfile(obj.casepath,['run' num2str(obj.run(ir))]);
                end
            end
            disp(obj.casepath)

            if exist(fullfile(obj.casepath,'body.txt'),'file')
                obj.casetype = 'gpu';
            else
                obj.casetype = 'cpu';
            end
            fprintf('Case type: %s\n', obj.casetype);

            rcase = read_case(casename, obj.casetype, obj.run);
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
            obj.blk.viewarea = [];
            if obj.NB == 9
                obj.topology = 1;
                obj.blk.oblocks = [3 5 7 4];
                obj.blk.oblocks_flip = [0 0 1 1];
            elseif obj.NB == 12
                obj.topology = 2;
                obj.blk.oblocks = [4 6 9 5];
                obj.blk.oblocks_flip = [0 0 1 1];
                obj.blk.viewarea = [-0.6 2 -0.5 0.5];
            else
                obj.topology = 3;
                obj.blk.oblocks = [1];
                obj.blk.oblocks_flip = [0];
                obj.blk.viewarea = [0 1 0 0.25];
            end
            
            if ~isempty(obj.blk.viewarea)
                obj.blk.aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                    (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];
            end

            if obj.NB == 12
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

        function writeCase(obj, nkproc)
            
            write_case(obj.casename, obj.blk, obj.next_block, obj.next_patch, ...
                obj.corner, obj.bcs, obj.gas, obj.solver, obj.topology, nkproc);
        end




        function slice = readSingleKSlice(obj,numslice)
            switch obj.casetype
                case 'cpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
                case 'gpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'kslice_time.txt'));
            end
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =kSlice(obj.blk,obj.gas,obj.runpath,slicenum);
        end

        function slice = readSingleJSlice(obj,numslice)
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
        end

        function readKSlices(obj, slicenums, runs)
            %READKSLICES Read in instantaneous k slices
            % Optional: numslices - only read last n slices if there are many
            
            exist("runs",'var')
            exist("numslices",'var')
            obj.kSlices = kSlice.empty;
            obj.nSlices = [];
            switch obj.casetype
                case 'cpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
                    if strcmp(obj.casepath, obj.runpath)
                        slices = dir(fullfile(obj.runpath,'kcu2_1_*'));
                        ishere = true;
                    else
                        slices = dir(fullfile(obj.runpath,'k_cuts','kcu2_1_*'));
                        ishere = false;
                    end
                case 'gpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'kslice_time.txt'));
                    if strcmp(obj.casepath, obj.runpath)
                        slices = dir(fullfile(obj.runpath,'kcut_1_*'));
                        ishere = true;
                    else
                        slices = dir(fullfile(obj.runpath,'k_cuts','kcut_1_*'));
                        ishere = false;
                    end
            end
            obj.nSlices = length(slices);
            for i=1:length(slices)
                inds(i) = str2num(slices(i).name(8:end));
            end
            [~,inds] = sort(inds);
            slices = slices(inds);
            if nargin == 1
                for i=1:obj.nSlices
                    fprintf('Reading slice %d/%d\n',[i obj.nSlices])
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.blk,obj.gas,obj.runpath,slicenum,obj.casetype,ishere);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("slicenums",'var')
                n=0;
                for i=slicenums
                    n=n+1;
                    fprintf('Reading slice %d/%d\n',[n length(slicenums)])
                    slices(i)
                    slicenum = str2num(slices(i).name(8:end));
                    obj.kSlices(i) = kSlice(obj.blk,obj.gas,obj.runpath,slicenum,obj.casetype,ishere);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("runs",'var')
                runs
                obj.nSlices = 0;
                slicenums=[];
                for nrun=runs
                    runpath = fullfile(obj.casepath,['run' num2str(nrun)]);
                    switch obj.casetype
                        case 'cpu'
                            slicetime = readmatrix(fullfile(runpath,'slice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcu2_1_*'));
                        case 'gpu'
                            slicetime = readmatrix(fullfile(runpath,'kslice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcut_1_*'));
                    end
                    obj.nSlices = obj.nSlices + length(slices);
                    read = 0;
                    for i=1:length(slices)
                        fprintf('Run %d: Reading slice %d/%d\n',[nrun i obj.nSlices])
                        read = read+1;
                        slicenum = str2num(slices(i).name(8:end))
                        slicenums(end+1) = slicenum;
                        obj.kSlices(end+1) = kSlice(obj.blk,obj.gas,runpath,slicenum,obj.casetype,ishere);
                    end
                end
                [~,inds] = sort([obj.kSlices.nSlice]);
                obj.kSlices = obj.kSlices(inds);

            end
        end

        function readJSlices(obj, runs, slicenums)
            %READJSLICES Read in instantaneous j slices
            % Optional: numslices - only read last n slices if there are many
            
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slices = dir(fullfile(obj.runpath,'j_cuts','jcu2_4_*'));
            %obj.nSlices = length(slices);
            obj.nSlices = length(slices);
            for i=1:length(slices)
                inds(i) = str2num(slices(i).name(8:end));
            end
            [~,inds] = sort(inds);
            slices = slices(inds);
            if nargin == 1
                for i=1:obj.nSlices
                    fprintf('Reading slice %d/%d\n',[i obj.nSlices])
                    slicenum = str2num(slices(i).name(8:end));
                    obj.jSlices(i) = jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
                    %slices(i).time = slicetime(i,2);
                end
            elseif exist("slicenums",'var') == 1
                n=0;
                for i=slicenums
                    n=n+1;
                    fprintf('Reading slice %d/%d\n',[n length(slicenums)])
                    slicenum = str2num(slices(i+obj.nSlices-slicenums(i)).name(8:end));
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
            obj.instFlow = volFlow(obj.runpath,obj.blk,obj.gas,obj.casetype);
            
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

        function readMeanFlows(obj, calc_s_unst)
            %READMEANFLOWS Read in and average 2D mean flows for multiple
            %runs
            if nargin < 2
                calc_s_unst = false;
            end
            regions = obj.getIntRegions;
            for ir = 1:length(obj.run)
                fprintf('Reading meanSlice %d/%d (run %d)\n', [ir length(obj.run) obj.run(ir)])
                mF = meanSlice(obj.runpaths{ir},obj.blk,obj.gas);
                if ir == 1
                    obj.meanFlow = mF;
                    if calc_s_unst
                        [e_unst_s, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                else
                    obj.meanFlow.addSlice(mF)
                end
                if ir == length(obj.run)
                    if calc_s_unst
                        [e_unst_e, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                end
                clear mF
            end
            if calc_s_unst
                obj.e_unst = e_unst_e - e_unst_s;
            end
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
            fclose(f);

%             for nProbe = 1:obj.nProbes
%                 temp = str2num(char(split(fgetl(f))));
%                 blkData.nb = temp(1);
%                 blkData.i = temp(2);
%                 blkData.j = temp(3);
%                 blkData.k = temp(4);
%                 blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
%                 blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
%                 obj.probes{nProbe} = dnsProbe(obj.runpath, nProbe, obj.nSkip, blkData, obj.gas);
%             end

            
            [obj.probes, obj.nProbes, obj.nSkip] = obj.readRunProbes(obj.runpaths{1});
            if length(obj.run) > 1
                for ir = 2:length(obj.run)

                    [probesNow, nProbesNow, nSkipNow] = obj.readRunProbes(obj.runpaths{ir});
                    if nProbesNow ~= obj.nProbes
                        fprintf('nProbes mismatch\n');
                        return
                    elseif nSkipNow ~= obj.nSkip
                        fprintf('nSkip mismatch\n');
                        return
                    else
                        for ip = 1:obj.nProbes
                            obj.probes{ip}.concatenate(probesNow{ip});
                        end
                    end
                end
            end
        end
        
        function [probes, nProbes, nSkip] = readRunProbes(obj, path)
            path
            f = fopen(fullfile(path,'probe.txt'), 'r');
            temp = str2num(char(split(fgetl(f))));
            nProbes = temp(1);
            nSkip = temp(2);
            for nProbe = 1:obj.nProbes
                temp = str2num(char(split(fgetl(f))));
                blkData.nb = temp(1);
                blkData.i = temp(2);
                blkData.j = temp(3);
                blkData.k = temp(4);
                blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
                blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
                probes{nProbe} = dnsProbe(path, nProbe, nSkip, blkData, obj.gas);
            end
            fclose(f);
        end

        function addProbe(obj, nb, i, j, k)
            blkData.nb = nb;
            blkData.i = i;
            blkData.j = j;
            if nargin < 5
                blkData.k = obj.solver.nk/2;
            else
                blkData.k = k;
            end
            blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
            blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
            obj.probes{end+1} = dnsProbe([], obj.nProbes+1, obj.nSkip, blkData, obj.gas);
            obj.nProbes = obj.nProbes + 1;
        end

        function addProbes(obj, newProbes)
            for i = 1:length(newProbes)
                blkData.nb = newProbes{i}.nb;
                blkData.i = newProbes{i}.i;
                blkData.j = newProbes{i}.j;
                blkData.k = newProbes{i}.k;
                blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
                blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
                obj.probes{end+1} = dnsProbe([], obj.nProbes+1, obj.nSkip, blkData, obj.gas);
            end
            obj.nProbes = obj.nProbes + length(newProbes);
        end

        function removeProbes(obj, rmProbes)
            for i = rmProbes
                obj.probes(i) = [];
                obj.nProbes = obj.nProbes-1;
            end
        end

        function plot_probes(obj, nProbe, ax)
            if nargin<3 || isempty(ax)
                ax = gca;
            end
            
            if nargin < 2 || isempty(nProbe)
                probes_to_plot = 1:obj.nProbes;
            else
                probes_to_plot = nProbe;
            end

            hold on
            C = colororder;
            sz = 25;
            for ip = probes_to_plot
                scatter(ax, obj.probes{ip}.x, obj.probes{ip}.y, sz, C(mod(ip-1,7)+1,:), 'filled')
            end
        end

        function writeProbeInput(obj, path)
            if nargin < 2
                path = obj.casepath;
            end
            f = fopen(fullfile(path, 'probe.txt'), 'w');

            fprintf(f,'%d %d\n',[obj.nProbes obj.nSkip]);
            for ip = 1:obj.nProbes
                fprintf(f,'%d %d %d %d\n', [obj.probes{ip}.nb obj.probes{ip}.i obj.probes{ip}.j obj.probes{ip}.k]);
            end
            
            fclose(f);
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
                s = pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis equal
            if ~isempty(obj.blk.viewarea)
                pbaspect(obj.blk.aspect);
                axis(obj.blk.viewarea);
            end
            if string(prop) == "schlieren"
                colormap(gray)
                map = colormap;
                map = flip(map,1);
                colormap(map);
                if nargin < 6
                    label = '$|\nabla \rho|/\rho$';
                end
            elseif string(prop) == "vortZ"
                colormap(redblue)
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

        function flipbook(obj,slices, prop, ax, lims, label)
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            if nargin < 5 || isempty(lims)
                lims = [];
                switch prop
                    case 'M'
                        lims = [0 1.6];
                    case 'tau_w'
                        lims = [0 800];
                    case 'vortZ'
                        lims = 1e5*[-0.7 0.7];
                end
            end
            if nargin < 6 || isempty(label)
                label = [];
                switch prop
                    case 'M'
                        label = '$M$';
                    case 'tau_w'
                        label = '$\tau_w$';
                    case 'vortZ'
                        label = '$\omega_z$';
                end
            end
            for i=1:length(slices)
                switch class(slices(i))
                    case 'kSlice'
                        obj.kPlot(slices(i),prop,ax,lims,label)
                    case 'jSlice'
                        obj.jPlot(slices(i),prop,ax,lims,label)
                end
                title(sprintf('Slice %d/%d', i, length(slices)))
                pause
            end
        end
                

        function kContour(obj,slice,prop,ax,levels,label,fmt,linew,falpha)
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            if nargin < 7 || isempty(fmt)
                fmt = 'k';
            end
            if nargin < 8 || isempty(linew)
                linew = 0.2;
            end
            if nargin < 9 || isempty(falpha)
                falpha = 1.0;
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
                s = contour(ax, obj.blk.x{i}, obj.blk.y{i}, q{i}, levels, 'k', 'LineWidth', linew, 'EdgeAlpha', falpha);
            end
            shading('interp')
            axis([-0.6 2 -0.5 0.5])
            axis equal
            if string(prop) == "schlieren"
                if nargin < 6
                    %label = '$|\nabla \rho|/\rho$';
                end
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
            pbaspect([5 2 1])

            axis([0.4 0.9 0 0.2])
            cb = colorbar('southoutside');
            if nargin > 4 && ~isempty(lims)
                caxis(lims)
            end
            nargin
            if nargin > 5 && ~isempty(label)
                label;
                %cb.Label.Interpreter = 'latex';
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

        function plot_blade(obj,fmt)
            if nargin < 2
                fmt = 'k';
            end
            xsurf = [];
            ysurf = [];
            for i=1:length(obj.blk.oblocks)
                xtemp = obj.blk.x{obj.blk.oblocks(i)}(:,end);
                ytemp = obj.blk.y{obj.blk.oblocks(i)}(:,end);
                if obj.blk.oblocks_flip(i) == 1
                    xtemp = flip(xtemp);
                    ytemp = flip(ytemp);
                end
                xsurf = [xsurf xtemp'];
                ysurf = [ysurf ytemp'];
            end
            plot(xsurf,ysurf,fmt)
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

        function setTrip(obj)
            mode = input('Trip mode:');
            obj.trip.mode = mode;
            x = input('x location:');
            switch mode
                case 1
                    obj.trip.amp = input('trip_amp:');
                    obj.trip.scale = input('trip_scale:');
                case {2, 3}
                    obj.trip.scale = input('trip_scale:');
                case 4
                    obj.trip.scale = input('ttrip_scale:');
                    obj.trip.nbumps = input('nbumps:');
            end

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
            if nargin < 2 || isempty(ax)
                ax = gca;
            end
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
                switch obj.trip.mode
                    case 1
                        obj.trip.amp = temp(1);
                        obj.trip.scale = temp(2);
                        xtmp = temp(3);
                    case {2,3}
                        obj.trip.scale = temp(1);
                        xtmp = temp(2);
                    case {4,5}
                        obj.trip.scale = temp(1);
                        xtmp = temp(2);
                        if length(temp) < 4
                            obj.trip.nbumps = str2num(fgetl(fid));
                        else
                            obj.trip.nbumps = temp(4);
                        end
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
                switch obj.trip.mode
                    case 1
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %20.16e\n', obj.trip.amp, obj.trip.scale, obj.trip.x, obj.trip.y);
                    case {2, 3}
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %20.16e\n', obj.trip.scale, obj.trip.x, obj.trip.y);
                    case 4
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %d\n', obj.trip.scale, obj.trip.x, obj.trip.y, obj.trip.nbumps);
                end
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

        function writeInputFiles(obj)
            write_input_files(obj.casename,obj.blk,obj.next_block,obj.next_patch,obj.corner,obj.bcs,obj.gas,obj.solver,obj.topology);
        end

        function writeGridFiles(obj)
            write_grid(obj.casename, obj.blk)
        end

        function get_cell_areas(obj)
            obj.cell_area = {};
            fprintf('Calculating cell areas\n');
            for ib = 1:obj.NB
                fprintf('Block %d\n',ib);
                for i=1:obj.blk.blockdims(ib,1)-1
                    for j=1:obj.blk.blockdims(ib,2)-1
                        xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
                            obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
                        ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
                            obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];

                        obj.cell_area{ib}(i,j) = abs(polyarea(xnow,ynow));
                    end
                end
            end
        end

        function [e_unst, e_conv] = mass_balance(obj)
            if isempty(obj.cell_area)
                obj.get_cell_areas;
            end
            
            prop = {};
            for ib=1:obj.NB
                [drou_dx, ~] = gradHO(obj.blk.x{ib}, obj.blk.y{ib}, obj.meanFlow.ro{ib}.*obj.meanFlow.u{ib});
                [~, drov_dy] = gradHO(obj.blk.x{ib}, obj.blk.y{ib}, obj.meanFlow.ro{ib}.*obj.meanFlow.v{ib});

                conv_prop{ib} = drou_dx + drov_dy;
                
                prop{ib}(:,:,1) = obj.meanFlow.ro{ib}.*obj.meanFlow.u{ib};
                prop{ib}(:,:,2) = obj.meanFlow.ro{ib}.*obj.meanFlow.v{ib};
                
            end

            regions = obj.getIntRegions;
            [conv_vol, conv_surf] = prop_balance(prop, obj.blk, obj.cell_area, regions);
            [conv_vol2, conv_surf2] = prop_balance(conv_prop, obj.blk, obj.cell_area, regions);

        end

        function [e_unst, e_s, e_phi, e_irrev, e_rev] = entropy_balance(obj, surface_int)
            if nargin < 2
                surface_int = false;
            end

            if isempty(obj.cell_area)
                obj.get_cell_areas;
            end

            for ib = 1:obj.NB
                s_flux{ib}(:,:,1) = obj.meanFlow.rous{ib};
                s_flux{ib}(:,:,2) = obj.meanFlow.rovs{ib};
                q_T_flux{ib}(:,:,1) = obj.meanFlow.rev_gen_x{ib};
                q_T_flux{ib}(:,:,2) = obj.meanFlow.rev_gen_y{ib};
            end

            regions = obj.getIntRegions;
            [e_s_vol, e_s_surf] = prop_balance(s_flux, obj.blk, obj.cell_area, regions);
            [e_phi, ~] = prop_balance(obj.meanFlow.diss_T, obj.blk, obj.cell_area, regions);
            [e_irrev, ~] = prop_balance(obj.meanFlow.irrev_gen, obj.blk, obj.cell_area, regions);
            [e_rev_vol, e_rev_surf] = prop_balance(q_T_flux, obj.blk, obj.cell_area, regions);
            e_unst = {};
            
            if surface_int
                e_s = e_s_surf;
                e_rev = e_rev_surf;
            else
                e_s = e_s_vol;
                e_rev = e_rev_vol;
            end 
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
%             [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            [~, e_s, e_phi, e_irrev, e_rev] = obj.entropy_balance;
            e_unst = obj.e_unst;
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

%             [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            [~, e_s, e_phi, e_irrev, e_rev] = obj.entropy_balance;
            if ~isempty(obj.e_unst)
                e_unst = obj.e_unst
            else
                e_unst = zeros(size(e_s));
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

        function interpInstFlow(obj, newcase)
            for ib = 1:obj.NB
                oldFlow = volFlowBlock(obj.runpath, ib, obj.blk, obj.gas, obj.casetype);
                newFlow = oldFlow.interpOntoNewGrid(newcase, ib);
                newFlow.writeFlow(newcase.casepath)
                clear oldFlow newFlow
            end
        end

        function updateBlk(obj, newblk)
            obj.blk.x = newblk.x;
            obj.blk.y = newblk.y;
            obj.blk.nk = newblk.nk;
            obj.blk.npp = newblk.npp;
            obj.blk.z = linspace(0, obj.solver.span, obj.blk.nk{1});
            obj.blk.blockdims = newblk.blockdims;

            y_inlet = [];
            for i=1:length(obj.blk.inlet_blocks{1})
                y_inlet = [y_inlet obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1)];
                ymidnow(i) = obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,ceil(end/2));
            end
            obj.y_inlet = sort(unique(y_inlet));
            obj.nj_inlet = length(obj.y_inlet);
            obj.inlet_width = max(obj.y_inlet) - min(obj.y_inlet);
        
        end

        function [i,j] = find_ij(obj, nb, x, y)
            xnow = obj.blk.x{nb};
            ynow = obj.blk.y{nb};
            ds = sqrt((xnow-x).^2 + (ynow-y).^2);
            [~,ind] = min(ds,[],'all');
            [i,j] = ind2sub(size(xnow),ind);
        end

        function newProbes = interp_probes(obj, oldCase)
            oldProbes = oldCase.probes;
            newProbes = {};
            for ip = 1:oldCase.nProbes
                [i,j] = obj.find_ij(oldProbes{ip}.nb,oldProbes{ip}.x,oldProbes{ip}.y);
                blkData.i = i; blkData.j = j;
                blkData.k = obj.solver.nk/2;
                blkData.x = obj.blk.x{oldProbes{ip}.nb}(i,j);
                blkData.y = obj.blk.y{oldProbes{ip}.nb}(i,j);
                blkData.nb = oldProbes{ip}.nb;
                newProbes{ip} = dnsProbe([],ip,obj.nSkip,blkData,obj.gas);
            end
        end


        function writeMovie(obj, slices, prop, lims, label, area)
            if nargin < 6 || isempty(area)
                area = obj.blk.viewarea;
                aspect = obj.blk.aspect;
            else
                aspect = [(area(2)-area(1))...
                    (area(4)-area(3)) 1];
            end
            if nargin < 4 || isempty(lims)
                switch prop
                    case 'M'
                        lims = [0 1.6];
                    case 'tau_w'
                        lims = [0 800];
                    case 'vortZ'
                        lims = 1e5*[-0.7 0.7];
                end
            end
            if nargin < 5 || isempty(label)
                switch prop
                    case 'M'
                        label = 'M';
                    case 'tau_w'
                        label = '\tau_w';
                    case 'vortZ'
                        label = '\omega_z';
                end
            end

            p = gcp('nocreate');
%             if isempty(p)
%                 parpool(16);
%             end
            switch class (slices(1))
                case 'kSlice'
                    var = ['k_' prop];
                case 'jSlice'
                    var = ['j_' prop];
            end
            imgfolder = fullfile(obj.runpath,'animation_images',var);
            %%
            if ~exist(imgfolder, 'dir')
                   mkdir(imgfolder);
            end
            
            %%
            outfolder = obj.runpath;
            for i=1:length(slices)
                fprintf('Plotting slice %d/%d\n',i,length(slices))
                if ~exist(fullfile(imgfolder,sprintf('img_%03d.png',slices(i).nSlice)),'file')
                    switch class(slices(i))
                        case 'kSlice'
                            slice2kPlot(slices(i), obj.blk, prop, fullfile(imgfolder,sprintf('img_%03d.png',slices(i).nSlice)), lims, label, area, aspect);
                        case 'jSlice'
		                    slice2jPlot(slices(i), prop, fullfile(imgfolder, sprintf('img_%03d.png',slices(i).nSlice)), lims, label);
                    end
                end
            end

            system(['ffmpeg -framerate 5 -pattern_type glob -i "' imgfolder '/*.png" '...
                '-c:v libx264 -profile:v high -pix_fmt yuv420p -vf ' ...
                '"pad=ceil(iw/2)*2:ceil(ih/2)*2" "' ...
                 outfolder '/run' num2str(obj.run(end)) '_' var '.mp4"']);



        end

        function dim = generateInflowTurb(obj, write)
            
            if nargin < 2
                write=true;
            end

            L = 1.0;            % Lengthscale of largest eddies
            M = 10000;          % Number of Fourier modes
            Tu = 0.005;         % Turbulence intensity

            [xb,yb,zb,dx] = writeTurbGrid(obj.blk.inlet_blocks{1}, obj.solver.span, obj.casepath);
            [~, lMax, minL, ~] = gridSpac(xb,yb,zb);
            [ni,nj,nk] = size(xb);
            dim = [ni,nj,nk];

            rgas = obj.gas.cp*(obj.gas.gam-1)/obj.gas.gam;
            Tin = obj.bcs.Toin - obj.bcs.vin^2/(2*obj.gas.cp);
            ronow = (obj.bcs.Poin/(rgas*Tin)) * (Tin/obj.bcs.Toin)^(obj.gas.gam/(obj.gas.gam-1));
            nunow = 1/ronow*obj.gas.mu_ref*((obj.gas.mu_tref+obj.gas.mu_cref)/(Tin+obj.gas.mu_cref))*(Tin/obj.gas.mu_tref)^1.5;

            if write
                dir = pwd;
                cd(fullfile(obj.casepath,'turb-gen'))
                f = fopen('turbin_params_grid.txt','w');
                fprintf(f,'%s\n','turb-grid.dat');
                fprintf(f,'%d, %d, %d\n',[ni nj nk]);
                fprintf(f,'%0.6e\n', lMax);
                fprintf(f,'%e\n', minL);
                fprintf(f,'%d\n',M);
                fprintf(f,'%e\n', nunow);
                fprintf(f,'%f, %f\n', [Tu obj.bcs.vin]);
                fprintf(f,'%f, %e\n', [-1 L]);
                fprintf(f,'%s\n', 'inflow_turb_gridTurb.dat');
                fclose(f)
                system('turbInGrid > out &')
                cd(dir)
            end
            
            fprintf('Set ilength to %d\n',ni);
            fprintf('Set lturb to: %9.7f\n', dx);

        end

        function [u,v,w] = conditionInflowTurb(obj,dim)
            gridPath = py.str(fullfile(obj.casepath,'turb-gen','turb-grid.dat'));
            namePath = py.str(fullfile(obj.casepath,'turb-gen','inflow_turb_gridTurb.dat'));
            DIM = py.tuple(int32(dim));
            turb_postproc = py.importlib.import_module('turb_postproc');
            U = turb_postproc.conditionAxSlc(namePath,gridPath,DIM,true,false,false);
            u = U(1); v = U(2); w = U(3);
            namePath = py.str(fullfile(obj.casepath,'turb-gen','inflow_turb.dat'));
            turb_postproc.write3DNSorder(namePath,u,v,w)
        end

        function write_HYDRA_mesh(obj)
            path = fullfile(obj.casepath,[obj.casename '_HYDRA.xyz']);
            nk = 6;
            span = (nk-1)*obj.solver.span/(obj.solver.nk-1);
            write_plot3d_extruded(obj.blk, path, nk, span);
        end

        function write_2d_plot3d_mesh(obj)
            path = fullfile(obj.casepath,[obj.casename '_2d.xyz']);
            nk = 1;
            span = (nk-1)*obj.solver.span/(obj.solver.nk-1);
            write_plot3d_2d(obj.blk, path);
        end
    end
end
