classdef kSlice
    % KSLICE Contains a 2D slice of the flow at k-boundary
    %   Detailed explanation goes here

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
        nSlice;
        blk;
        velDash;
    end

    properties (Dependent = true)
        T;
        p;
        M;
        s;
        vel;
        vortZ;
    end

    methods
        function obj = kSlice(casedir, nSlice, blk, gas)
        
            if nargin > 0
                obj.gam = gas.gam;
                obj.cp = gas.cp;
                obj.rgas = obj.cp*(1-1/obj.gam);
                blockdims = blk.blockdims;
                obj.NB = size(blockdims,1);
                obj.blk = blk;
    
                if nargin > 2
                    obj.nSlice = nSlice;
        
                    for nb = 1:obj.NB
                        flopath = fullfile(casedir, 'k_cuts',  ['kcu2_' num2str(nb) '_' num2str(nSlice)]);
                        flofile = fopen(flopath,'r');
                        nodfile = fopen(fullfile(casedir, 'k_cuts', ['knd2_' num2str(nb) '_' num2str(nSlice)]),'r');
                        A = fread(flofile,inf,'float64');
                        A = reshape(A,5,length(A)/5);
                        
                        B = fread(nodfile,inf,'uint32');
                        B = reshape(B,3,length(B)/3);
                
                        fclose(flofile);
                        fclose(nodfile);
        
                        ro = zeros(blockdims(nb,1),blockdims(nb,2));
                        ru = zeros(blockdims(nb,1),blockdims(nb,2));
                        rv = zeros(blockdims(nb,1),blockdims(nb,2));
                        rw = zeros(blockdims(nb,1),blockdims(nb,2));
                        Et = zeros(blockdims(nb,1),blockdims(nb,2));
        
                        for n=1:size(A,2)
                            i = B(1,n);
                            j = B(2,n);
                            k = B(3,n);
                            ro(i,j) = A(1,n);
                            ru(i,j) = A(2,n);
                            rv(i,j) = A(3,n);
                            rw(i,j) = A(4,n);
                            Et(i,j) = A(5,n);
                        end
        
                        obj.ro{nb} = ro;
                        obj.u{nb} = ru./ro;
                        obj.v{nb} = rv./ro;
                        obj.w{nb} = rw./ro;
                        obj.Et{nb} = Et;
                    end
                end
            end
        end

        function value = get.p(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.p{nb}./(obj.ro{nb}*obj.rgas);
            end
        end

        function value = get.vel(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.M(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.vel{nb}./sqrt(obj.gam*obj.rgas*obj.T{nb});
            end
        end

        function value = get.s(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.cp*log(obj.T{nb}/300) - obj.rgas*log(obj.p{nb}/1e5);
            end
        end

        function value = get.vortZ(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                [~,DUDY] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.u{nb});
                [DVDX,~] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.v{nb});
                value{nb} = DVDX-DUDY;
            end
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