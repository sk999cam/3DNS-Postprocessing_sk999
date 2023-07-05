classdef DNS_channel < DNS_case
    %DNS_CHANNEL Subclass of DNS_case contaning methods and properties
    ... specific to channel flow cases

    properties
        nbi;
        nbj;
    end

    properties (Dependent = true)
        Re_theta_in
    end

    methods
        function obj = DNS_channel(casename,run)
            %DNS_CHANNEL Construct an instance of this class

            args.topology = 3;
            if nargin > 0
                if nargin < 2
                    run = [];
                end
            else
                casename = [];
                run = [];
            end
            obj@DNS_case(casename,run,args);

            if nargin > 0 && ~isempty(casename)
                obj.compute_blk_metadata;
            end
        end

        function compute_blk_metadata(obj)
            obj.blk.viewarea = [inf -inf inf -inf];
            for ib = 1:obj.NB
                obj.blk.viewarea(1) = min(obj.blk.viewarea(1), min(obj.blk.x{ib},[],'all'));
                obj.blk.viewarea(2) = max(obj.blk.viewarea(2), max(obj.blk.x{ib},[],'all'));
                obj.blk.viewarea(3) = min(obj.blk.viewarea(3), min(obj.blk.y{ib},[],'all'));
                obj.blk.viewarea(4) = max(obj.blk.viewarea(4), max(obj.blk.y{ib},[],'all'));
            end
            obj.blk.aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];

            obj.nbi = 1;
            while obj.blk.next_block{obj.nbi}.ip ~= 0
                obj.nbi = obj.nbi+1;
            end
            obj.nbj = obj.NB/obj.nbi;
        end

        function split_domain(obj, newCase)
            newFlow = volFlow();
            newFlow.flowpath = newCase.casepath;
            if isempty(obj.instFlow)
                obj.readInstFlow;
            end
            ims = ones(newCase.nbi, newCase.nbj);
            ips = ones(newCase.nbi, newCase.nbj);
            jms = ones(newCase.nbi, newCase.nbj);
            jps = ones(newCase.nbi, newCase.nbj);

            X = obj.concat_prop(obj.blk.x);
            Y = obj.concat_prop(obj.blk.y);

            for j = 1:newCase.nbj
                ib = 1+(j-1)*newCase.nbi;
                xc = newCase.blk.x{ib}(end,end);
                yc = newCase.blk.y{ib}(end,end);
                dist = sqrt((X-xc).^2+(Y-yc).^2);
                [~, ind] = min(dist,[],'all');
                [ic, jc] = ind2sub(size(X), ind);
                jps(:,j) = jc;
            end

            for i = 1:newCase.nbi
                ib = i;
                xc = newCase.blk.x{ib}(end,end);
                yc = newCase.blk.y{ib}(end,end);
                dist = sqrt((X-xc).^2+(Y-yc).^2);
                [~, ind] = min(dist,[],'all');
                [ic, jc] = ind2sub(size(X), ind);
                ips(i,:) = ic;
            end
            jms(:,2:end) = jps(:,1:end-1);
            ims(2:end,:) = ips(1:end-1,:);

            fkc = linspace(0,1,obj.blk.nk);
            fk = linspace(0,1,newCase.blk.nk);
            props = {'ro','u','v','w','Et'};
            for ip = 1:length(props)
                prop = props{ip};
                fprintf('Interpolating %s\n', prop);
                propnow = obj.concat_prop(obj.instFlow.(prop));
                
                for j = 1:newCase.nbj
                    for i = 1:newCase.nbi
                        ib = i+(j-1)*newCase.nbi;
                        fprintf('Block %d\n', ib);
                        Xn = X(ims(i,j):ips(i,j),jms(i,j):jps(i,j));
                        Yn = Y(ims(i,j):ips(i,j),jms(i,j):jps(i,j));
                        Vn = propnow(ims(i,j):ips(i,j),jms(i,j):jps(i,j),:);
                        [fic, fjc] = newFlow.get_spacing(Xn,Yn);
                        [fi,fj] = newFlow.get_spacing(newCase.blk.x{ib}, newCase.blk.y{ib});
                        [Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
                        [J,I,K] = meshgrid(fj,fi,fk);
                        newFlow.(prop){ib} = interp3(Jc,Ic,Kc,Vn,J,I,K);
                    end
                end
                clear propnow
            end
            newFlow.blk = newCase.blk;
            newFlow.gas = newCase.gas;
            newFlow.NB = newCase.NB;
            newFlow.casetype = obj.casetype;
            newCase.instFlow = newFlow;
        end

        function arr = concat_prop(obj, prop)
            arr = [];
            for j = 1:obj.nbj
                row = [];
                for i=1:obj.nbi
                    ib = i+(j-1)*obj.nbi;
                    row = [row prop{ib}];
                end
                arr = [arr; row];
            end
        end

        function newCase = instantiate(obj)
            newCase = DNS_channel;
        end

        function write_hydra_inlet_bc(obj)
            kprof = obj.meanFlow.BLprof(0.25,'k');
            omprof = obj.meanFlow.BLprof(0.25,'omega_opt');
            [vel_prof, po_prof, To_prof, T_prof] = blasius_bl(obj.bcs.Toin, obj.bcs.vin, obj.bcs.theta, obj.blk.y{1}(1,:), obj.gas);
            Toin = To_prof*obj.bcs.Toin;
            Tin = obj.bcs.Toin - obj.bcs.vin^2/(2*obj.gas.cp);
            Poin = po_prof*obj.bcs.Poin;
            rgas = obj.gas.cp*(obj.gas.gam-1)/obj.gas.gam;
            Mprof = obj.bcs.vin*vel_prof./sqrt(obj.gas.gam*rgas*Tin*T_prof);

            f = fopen(fullfile(obj.casepath, 'hydra_inlet_bc_data.txt'),'w');
            for j=1:length(kprof)
                data = fprintf(f,'%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                    [obj.blk.y{1}(1,j) Toin(j) Poin(j) 0.0 Mprof(j) kprof(j) omprof(j)]);
            end
            fclose(f);
        end

        function boundaries = getBoundaries(obj)


            % BCs:
            % 3: Wall
            % 4: Pressure inlet
            % 5: Pressure outlet
            % 12: Periodic
            % 8: Periodic shadow

            boundaries = {};
            b.label = "Inlet";
            blks = [];
            for j = 1:obj.nbj
                ib = 1+(j-1)*obj.nbi;
                blks = [blks ib];
            end
            b.blocks = blks;
            b.patches(1:length(blks)) = 1;
            b.type = 4;
            boundaries{end+1} = b;
            
            b.label = "Wall";
            b.blocks = 1:obj.nbi;
            b.patches(1:obj.nbi) = 3;
            b.type = 3;
            boundaries{end+1} = b;
            
            b.label = "Pre-shock";
            xmid = obj.blk.x{obj.nbi}(end,1)/2;
            blks = [];
            psblks = [];
            for ib=(obj.NB-obj.nbi)+1:obj.NB
                if obj.blk.x{ib}(floor(obj.blk.blockdims(ib,1)/2),end) < xmid
                    blks = [blks ib];
                else
                    psblks = [psblks ib];
                end
            end
            b.blocks = blks;
            b.patches(1:length(blks)) = 4;
            b.type = 5;
            boundaries{end+1} = b;

            b.label = "Post-shock";
            b.blocks = psblks;
            b.patches(1:length(psblks)) = 4;
            b.type = 5;
            boundaries{end+1} = b;
            
            b.label = "Outlet";
            blks = [];
            for j = 1:obj.nbj
                ib = j*obj.nbi;
                blks = [blks ib];
            end
            b.blocks = blks;
            b.patches(1:length(blks)) = 2;
            b.type = 5;
            boundaries{end+1} = b;
        end


    end
end