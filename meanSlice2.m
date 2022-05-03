classdef meanSlice2
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        NB;
        gam;
        cp;
        rgas;
        ro;
        u;
        v;
        w;
        Et;
        time;
        nMean;
        meanTime;
        oblocks;
        oblocks_flip;
        p0in;
    end

    properties (Dependent = true)
        T;              % Temperature
        p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        Msurf;          % Surface Mach No
        Pr;             % Turbulence production
        
    end

    methods
        function obj = meanSlice2(casedir, blk, gas)

            nstats = 17;
            obj.oblocks = blk.oblocks;
            obj.oblocks_flip = blk.oblocks_flip;
            
            if nargin > 0
                obj.gam = gas.gam;
                obj.cp = gas.cp;
                obj.rgas = obj.cp*(1-1/obj.gam);
    
                obj.NB = size(blk.blockdims,1);
                fid = fopen(fullfile(casedir, 'mean_flo', 'mean_time.txt'));
                while ~feof(fid) % Use lastest mean files
                    temp=fgetl(fid);
                end
                fclose(fid);
                temp = str2num(temp);
                obj.nMean = temp(1);
                obj.meanTime = temp(3);
                for nb = 1:obj.NB

                    flopath = fullfile(casedir, 'mean_flo',  ['mean2_' num2str(nb) '_' num2str(obj.nMean)]);
                    flofile = fopen(flopath,'r');
                    nodfile = fopen(fullfile(casedir, 'mean_flo', ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]),'r');
                    A = fread(flofile,inf,'float64');
                    A = reshape(A,nstats,length(A)/nstats);
                    
                    B = fread(nodfile,inf,'uint32');
                    B = reshape(B,3,length(B)/3);
            
                    fclose(flofile);
                    fclose(nodfile);
    
                    rodt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rudt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rvdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rwdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    Etdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    ro2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rou2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rov2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    row2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    rouv = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rouw = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rovw = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    p2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    p = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    ros = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    diss = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    
                    for n=1:size(A,2)
                        i = B(1,n);
                        j = B(2,n);
                        k = B(3,n);
                        rodt(i,j) = A(1,n);
                        rudt(i,j) = A(2,n);
                        rvdt(i,j) = A(3,n);
                        rwdt(i,j) = A(4,n);
                        Etdt(i,j) = A(5,n);

                        ro2(i,j) = A(6,n);
                        rou2(i,j) = A(7,n);
                        rov2(i,j) = A(8,n);
                        row2(i,j) = A(9,n);

                        rouv(i,j) = A(10,n);
                        rouw(i,j) = A(11,n);
                        rovw(i,j) = A(12,n);

                        p2(i,j) = A(13,n);
                        p(i,j) = A(14,n);
                        T(i,j) = A(15,n);
                        ros(i,j) = A(16,n);
                        diss(i,j) = A(17,n);

                    end
    
                    obj.ro{nb} = rodt/obj.meanTime;
                    obj.u{nb} = rudt./(rodt);
                    obj.v{nb} = rvdt./(rodt);
                    obj.w{nb} = rwdt./(rodt);
                    obj.Et{nb} = Etdt/obj.meanTime;
                end

            end
            Mnow = obj.M;
            %Mnow = Mnow{blk.inlet_blocks{1}};
            pnow = obj.p;
            %pnow = pnow{blk.inlet_blocks{1}};
            
            p0 = [];
            for i=1:length(blk.inlet_blocks{1})
               p0now = pnow{blk.inlet_blocks{1}(i)}.*(1+((obj.gam - 1)/2)*Mnow{blk.inlet_blocks{1}(i)}.^2).^(obj.gam/(obj.gam-1));
               p0 = [p0 p0now(40:100,:)];
            end
            obj.p0in = mean(p0,'all');
        end

        function value = get.p(obj)
            disp('Calculating p')
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            disp('Calculating T')
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                value{nb} = pnow./(obj.ro{nb}*obj.rgas);
            end
        end

        function value = get.vel(obj)
            disp('Calculating vel')
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.M(obj)
            disp('Calculating M')
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                pnow = (obj.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.rgas);
                velnow = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
                value{nb} = velnow./sqrt(obj.gam*obj.rgas*Tnow);
            end
        end

        function value = get.s(obj)
            disp('Calculating s')
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                pnow = (obj.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.rgas);
                value{nb} = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
            end
        end

        function value = get.Msurf(obj)
            disp('Calculating surface M')
            psurf = [];
            pnow = obj.p;
            size(pnow{4})
            for i=1:length(obj.oblocks)
                clear temp
                temp = pnow{obj.oblocks(i)}(:,end);
                size(temp)
                if obj.oblocks_flip(i) == 1
                    i
                    temp = flip(temp);
                end
                psurf = [psurf temp'];
            end
            value = sqrt((2/(obj.gam - 1)) * ( (psurf/obj.p0in).^(-(obj.gam-1)/obj.gam) - 1));
        end


        function getSize(obj)
            props = properties(obj); 
            totSize = 0; 
            
           
            for ii=1:length(props) 
                currentProperty = getfield(obj, char(props(ii))); 
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end
    end
end