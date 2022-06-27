classdef volFlow
    % VOLFLOW Contains a 3D instntaneous flow from DNS case

        
    properties
        NB;
        gas;
        blk;
        gam;
        cp;
        rgas;
        ro;
        u;
        v;
        w;
        Et;
        time;
        tau;
        div;
        vort_x;
        vort_y;
        vort_z;
        rmut;
        nk;
    end

    properties (Dependent = true)
        T;              % Temperature
        p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        Msurf;          % Surface Mach No
        MSlice;
        
        
    end

    methods
        function obj = volFlow(casedir, blk, gas)
            
            if nargin > 0
                blockdims = blk.blockdims;
                obj.blk = blk;
                obj.gas = gas;
                obj.nk = blockdims(1,3);
                obj.gam = gas.gam;
                obj.cp = gas.cp;
                obj.rgas = obj.cp*(1-1/obj.gam);
    
                obj.NB = size(blockdims,1);
                for nb = 1:obj.NB
                    

                    flopath = fullfile(casedir,  ['flo2_' num2str(nb)])
                    flofile = fopen(flopath,'r');
                    nodfile = fopen(fullfile(casedir, ['nod2_' num2str(nb)]),'r');
                    %viscpath = fullfile(casedir,  ['visc_' num2str(nb)]);
                    %viscfile = fopen(viscpath,'r');

                    A = fread(flofile,inf,'float64');
                    A = reshape(A,5,length(A)/5);
                    
                    B = fread(nodfile,inf,'uint32');
                    B = reshape(B,3,length(B)/3);

%                     C = fread(viscfile,inf,'float64');
%                     C = reshape(C,13,length(C)/13 );

            
                    fclose(flofile);
                    fclose(nodfile);
    
                    ro = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    ru = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    rv = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    rw = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    Et = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xx = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_yy = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_zz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xy = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_yz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
    
                    for n=1:size(A,2)
                        i = B(1,n);
                        j = B(2,n);
                        k = B(3,n);
                        ro(i,j,k) = A(1,n);
                        ru(i,j,k) = A(2,n);
                        rv(i,j,k) = A(3,n);
                        rw(i,j,k) = A(4,n);
                        Et(i,j,k) = A(5,n);
                    end
    
                    obj.ro{nb} = ro;
                    obj.u{nb} = ru./ro;
                    obj.v{nb} = rv./ro;
                    obj.w{nb} = rw./ro;
                    obj.Et{nb} = Et;
                end
            else
                obj.ro = {};
                obj.u = {};
                obj.v = {};
                obj.w = {};
                obj.Et = {};
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
            pool = gcp('nocreate');
            value = cell(1,obj.NB);
            if isempty(pool)
                for nb =1:obj.NB
                    pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                    Tnow = pnow./(obj.ro{nb}*obj.rgas);
                    value{nb} = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                end
            else
                parfor nb =1:obj.NB
                    pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb})
                    Tnow = pnow./(obj.ro{nb}*obj.rgas);
                    value{nb} = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                end
            end
        end

        function value = get.Msurf(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.cp*log(obj.T{nb}/300) - obj.rgas*log(obj.p{nb}/1e5);
            end
        end

        function value = k_ave(obj,prop)
            pool = gcp('nocreate');
            if isempty(pool)
                for nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                            Tnow = pnow./(obj.ro{nb}*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value{nb} = mean(propnow,3);
                end
            else
                parfor nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                            Tnow = pnow./(obj.ro{nb}*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value{nb} = mean(propnow,3);
                end
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

        function value = slice(obj,k)
            if nargin == 1
                k = floor(obj.nk/2);
            end

            value = flowSlice(obj.blk, obj.gas);
            for nb=1:obj.NB
                nb
                value.ro{nb} = obj.ro{nb}(:,:,k);
                value.u{nb} = obj.u{nb}(:,:,k);
                value.v{nb} = obj.v{nb}(:,:,k);
                value.w{nb} = obj.w{nb}(:,:,k);
                value.Et{nb} = obj.Et{nb}(:,:,k);
            end
        end

        function write_tecplot_files(obj)
            for nb = 1:obj.NB
                if(nb==4 | nb==6)
                    tdata=[];
                    tdata.Nvar=9;
                    tdata.varnames={'x','y','z','p','T','ro','u','v','w'};
                    tdata.cubes(1).zonename=['block ',num2str(nb)];
                    zz(1,1,1:obj.blk.nk{nb}) = linspace(0,obj.blk.span,obj.blk.nk{nb}); 
                    z = repmat(zz,[obj.blk.blockdims(nb,1) obj.blk.blockdims(nb,2) 1]);
                    x = repmat(obj.blk.x{nb},[1 1 obj.blk.nk{nb}]);
                    y = repmat(obj.blk.y{nb},[1 1 obj.blk.nk{nb}]);
                    tdata.cubes(1).x=x(end:-1:1,:,:);
                    tdata.cubes(1).y=y(end:-1:1,:,:);
                    tdata.cubes(1).z=z(end:-1:1,:,:);
                    tdata.cubes(1).v(1,:,:,:)=obj.p{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(2,:,:,:)=obj.T{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(3,:,:,:)=obj.ro{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(4,:,:,:)=obj.u{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(5,:,:,:)=obj.v{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(6,:,:,:)=obj.w{nb}(end:-1:1,:,:);
                    tdata.vformat(1:9) = 1; 
                    tdata.cubes.solutiontime=0;
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])

                else
                    tdata=[];
                    tdata.Nvar=9;
                    tdata.varnames={'x','y','z','p','T','ro','u','v','w'};
                    tdata.cubes(1).zonename=['block ',num2str(nb)];
                    zz(1,1,1:obj.blk.nk{nb}) = linspace(0,obj.blk.span,obj.blk.nk{nb}); 
                    z = repmat(zz,[obj.blk.blockdims(nb,1) obj.blk.blockdims(nb,2) 1]);
                    x = repmat(obj.blk.x{nb},[1 1 obj.blk.nk{nb}]);
                    y = repmat(obj.blk.y{nb},[1 1 obj.blk.nk{nb}]);
                    tdata.cubes(1).x=x;
                    tdata.cubes(1).y=y;
                    tdata.cubes(1).z=z;
                    tdata.cubes(1).v(1,:,:,:)=obj.p{nb};
                    tdata.cubes(1).v(2,:,:,:)=obj.T{nb};
                    tdata.cubes(1).v(3,:,:,:)=obj.ro{nb};
                    tdata.cubes(1).v(4,:,:,:)=obj.u{nb};
                    tdata.cubes(1).v(5,:,:,:)=obj.v{nb};
                    tdata.cubes(1).v(6,:,:,:)=obj.w{nb};
                    tdata.vformat(1:9) = 1; 
                    tdata.cubes.solutiontime=0;
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])
        
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])
                end
            end
        end

        function newFlow = interpOntoNewGrid(obj, newblk, nk)

            newFlow = volFlow();
            ronow = {};


            for ib = 1:obj.NB
                sprintf('interpolting block: %d', ib)
                if ~ismember(ib, obj.blk.oblocks)
                    ib
                    if obj.blk.x{ib}(1,1) < 0.2; iq = 1; else iq = size(obj.blk.x{ib},1); end
                    if obj.blk.y{ib}(1,1) < 0; jq = 1; else jq = size(obj.blk.x{ib},2); end
                    fic = (obj.blk.x{ib}(:,jq) - obj.blk.x{ib}(1,jq))/(obj.blk.x{ib}(end,jq) - obj.blk.x{ib}(1,jq));
                    fjc = (obj.blk.y{ib}(iq,:) - obj.blk.y{ib}(iq,:))/(obj.blk.y{ib}(iq,end) - obj.blk.y{ib}(iq,1));
                    fic = fic';

                    

                    if newblk{ib}.x(1,1) < 0.2; iq = 1; else iq = size(newblk{ib}.x,1); end
                    if newblk{ib}.y(1,1) < 0; jq = 1; else jq = size(newblk{ib}.y,2); end
                    fi = (newblk{ib}.x(:,jq) - newblk{ib}.x(1,jq))/(newblk{ib}.x(end,jq) - newblk{ib}.x(1,jq));
                    fj = (newblk{ib}.y(iq,:) - newblk{ib}.y(iq,1))/(newblk{ib}.y(iq,end) - newblk{ib}.y(iq,1));
                    fi = fi';
                else
                    fic = linspace(0,1,size(obj.blk.x{ib},1));
                    fjc = linspace(0,1,size(obj.blk.x{ib},2));
                    fi = linspace(0,1,size(newblk{ib}.x,1));
                    fj = linspace(0,1,size(newblk{ib}.x,2));
                end

                fkc = linspace(0,1,obj.blk.nk{1});
                fk = linspace(0,1,nk);


                [Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
                [J,I,K] = meshgrid(fj,fi,fk);

                ronow{ib} = interp3(Jc,Ic,Kc,obj.ro{ib},J,I,K);
                unow{ib} = interp3(Jc,Ic,Kc,obj.u{ib},J,I,K);
                vnow{ib} = interp3(Jc,Ic,Kc,obj.v{ib},J,I,K);
                wnow{ib} = interp3(Jc,Ic,Kc,obj.w{ib},J,I,K);
                Etnow{ib} = interp3(Jc,Ic,Kc,obj.Et{ib},J,I,K);
            end
            newFlow.ro = ronow;
            newFlow.u = unow;
            newFlow.v = vnow;
            newFlow.w = wnow;
            newFlow.Et = Etnow;
        end

        function writeFlow(obj, path)


            for nb = 1:obj.NB
                ronow = obj.ro{nb};
                runow = ronow*obj.u{nb};
                rvnow = ronow*obj.v{nb};
                rwnow = ronow*obj.w{nb};
                Etnow = obj.Et{nb};

                fpos = 0;
                for i=1:size(obj.ro{nb},1)
                    for j=1:size(obj.ro{nb},1)
                        for k=1:size(obj.ro{nb},1)
                            fpos = fpos+1;
                            B(1, fpos) = i;
                            B(2, fpos) = j;
                            B(3,fpos) = k;

                            A(1,n) = ronow(i,j,k);
                            A(2,n) = runow(i,j,k);
                            A(3,n) = rvnow(i,j,k);
                            A(4,n) = rwnow(i,j,k);
                            A(5,n) = Etnow(i,j,k);
                        end
                    end
                end

                flopath = fullfile(path,  ['flo2_' num2str(nb)]);
                flofile = fopen(flopath,'w');
                nodfile = fopen(fullfile(path, ['nod2_' num2str(nb)]),'w');
                %viscpath = fullfile(casedir,  ['visc_' num2str(nb)]);
                %viscfile = fopen(viscpath,'r');
    
                fwrite(flofile,A,'float64');
                fwrite(nodfile,B,'uint32');
    
        
                fclose(flofile);
                fclose(nodfile);
                    
            end
        end 
    end
end