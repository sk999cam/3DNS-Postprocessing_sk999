classdef volFlowBlock
    % VOLFLOW Contains a 3D instntaneous flow from DNS case

        
    properties
        ib;
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
        span;
        z;
        Qstore;
    end

    properties (Dependent = true)
        T;              % Temperature
        p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        Msurf;          % Surface Mach No
        MSlice;
        Q;              % Q criterion
        
    end

    methods
        function obj = volFlowBlock(casedir, ib, blk, gas)
            
            if nargin > 0
                obj.ib = ib;
                obj.blk.blockdims = blk.blockdims(ib,:);
                obj.blk.x = blk.x{ib};
                obj.blk.y = blk.y{ib};
                obj.gas = gas;
                obj.nk = obj.blk.blockdims(1,3);
                obj.gam = gas.gam;
                obj.cp = gas.cp;
                obj.rgas = obj.cp*(1-1/obj.gam);
                    
                fprintf('Reading block %d\n',ib)
                flopath = fullfile(casedir,  ['flo2_' num2str(ib)]);
                flofile = fopen(flopath,'r');
                nodfile = fopen(fullfile(casedir, ['nod2_' num2str(ib)]),'r');
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

                ro = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                ru = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                rv = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                rw = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                Et = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_xx = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_yy = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_zz = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_xy = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_xz = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));
                tau_yz = zeros(obj.blk.blockdims(1),obj.blk.blockdims(2),obj.blk.blockdims(3));

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

                obj.ro = ro;
                obj.u = ru./ro;
                obj.v = rv./ro;
                obj.w = rw./ro;
                obj.Et = Et;
                
            else
                obj.ro = {};
                obj.u = {};
                obj.v = {};
                obj.w = {};
                obj.Et = {};
            end
        end

        function value = get.p(obj)
            value = (obj.gam -1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
        end

        function value = get.T(obj)
            value = obj.p./(obj.ro*obj.rgas);
        end

        function value = get.vel(obj)
            value = sqrt(obj.u.^2 + obj.v.^2 + obj.w.^2);
        end

        function value = get.M(obj)
            value = obj.vel./sqrt(obj.gam*obj.rgas*obj.T);
        end

        function value = get.s(obj)
            pnow = (obj.gam -1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
            Tnow = pnow./(obj.ro*obj.rgas);
            value = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
        end

        function value = get.Q(obj)
            if ~isempty(obj.Qstore)
                value = obj.Qstore;
            else
                value = cell(1,obj.NB);
                for nb = obj.blk.oblocks
                    fprintf('Calculating Q criterion in block %d\n',obj.ib)
                    value = Q_criterion(obj.blk.x, obj.blk.y, obj.blk.z, obj.u, obj.v, obj.w);
                end
                obj.Qstore = value;
            end
        end

        function s = plot_Q_criterion(obj, thresh)
%             X = repmat(obj.blk.x{1},[1 1 obj.nk]);
%             Y = repmat(obj.blk.y{1},[1 1 obj.nk]);
%             Z = repmat(obj.blk.z,[size(obj.blk.x{1}) 1]);
            Qnow = obj.Q;
            for nb = 5%obj.blk.oblocks
                %X = repmat(obj.blk.x,1,1,length(obj.blk.z));
                X = linspace(0,1,size(obj.blk.x,1));
                %Y = repmat(obj.blk.y,1,1,length(obj.blk.z));
                Y = linspace(0,1,size(obj.blk.y,2));
                %znow = reshape(obj.blk.z,1,1,[]);
                %Z = repmat(znow, size(X,1), size(X,2));
                Z = linspace(0,1,length(obj.blk.z));
                %[X, Y, Z] = meshgrid(X,Y,Z);
                Qblock = permute(Qnow,[2 1 3]);
                s = isosurface(X, Y, Z, Qblock, thresh);
                %p = patch(s);
            end
            %axis equal
        end

        function value = k_ave(obj,prop)
            pool = gcp('nocreate');
            if isempty(pool)
                for nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
                            Tnow = pnow./(obj.ro*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value = mean(propnow,3);
                end
            else
                parfor nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
                            Tnow = pnow./(obj.ro*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value = mean(propnow,3);
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

        function write_tecplot_files(obj)
            for nb = 1:obj.NB
                if(nb==4 | nb==6)
                    tdata=[];
                    tdata.Nvar=9;
                    tdata.varnames={'x','y','z','p','T','ro','u','v','w'};
                    tdata.cubes(1).zonename=['block ',num2str(nb)];
                    zz(1,1,1:obj.nk) = linspace(0,obj.blk.span,obj.nk); 
                    z = repmat(zz,[obj.blk.blockdims(1) obj.blk.blockdims(2) 1]);
                    x = repmat(obj.blk.x,[1 1 obj.nk]);
                    y = repmat(obj.blk.y,[1 1 obj.nk]);
                    tdata.cubes(1).x=x(end:-1:1,:,:);
                    tdata.cubes(1).y=y(end:-1:1,:,:);
                    tdata.cubes(1).z=z(end:-1:1,:,:);
                    tdata.cubes(1).v(1,:,:,:)=obj.p(end:-1:1,:,:);
                    tdata.cubes(1).v(2,:,:,:)=obj.T(end:-1:1,:,:);
                    tdata.cubes(1).v(3,:,:,:)=obj.ro(end:-1:1,:,:);
                    tdata.cubes(1).v(4,:,:,:)=obj.u(end:-1:1,:,:);
                    tdata.cubes(1).v(5,:,:,:)=obj.v(end:-1:1,:,:);
                    tdata.cubes(1).v(6,:,:,:)=obj.w(end:-1:1,:,:);
                    tdata.vformat(1:9) = 1; 
                    tdata.cubes.solutiontime=0;
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])

                else
                    tdata=[];
                    tdata.Nvar=9;
                    tdata.varnames={'x','y','z','p','T','ro','u','v','w'};
                    tdata.cubes(1).zonename=['block ',num2str(nb)];
                    zz(1,1,1:obj.nk) = linspace(0,obj.blk.span,obj.nk); 
                    z = repmat(zz,[obj.blk.blockdims(1) obj.blk.blockdims(2) 1]);
                    x = repmat(obj.blk.x,[1 1 obj.nk]);
                    y = repmat(obj.blk.y,[1 1 obj.nk]);
                    tdata.cubes(1).x=x;
                    tdata.cubes(1).y=y;
                    tdata.cubes(1).z=z;
                    tdata.cubes(1).v(1,:,:,:)=obj.p;
                    tdata.cubes(1).v(2,:,:,:)=obj.T;
                    tdata.cubes(1).v(3,:,:,:)=obj.ro;
                    tdata.cubes(1).v(4,:,:,:)=obj.u;
                    tdata.cubes(1).v(5,:,:,:)=obj.v;
                    tdata.cubes(1).v(6,:,:,:)=obj.w;
                    tdata.vformat(1:9) = 1; 
                    tdata.cubes.solutiontime=0;
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])
        
                    mat2tecplot(tdata,['tec_flow_',num2str(nb),'.plt'])
                end
            end
        end

        function newFlow = interpOntoNewGrid(obj, newcase, ib)


            newFlow = volFlowBlock();
            ronow = {};

            fprintf('interpolting block: %d\n', obj.ib)

            [fic, fjc] = obj.get_spacing(obj.blk.x,obj.blk.y);
            [fi, fj] = obj.get_spacing(newcase.blk.x{obj.ib},newcase.blk.y{obj.ib});
            fkc = linspace(0,1,obj.blk.nk);
            fk = linspace(0,1,newcase.blk.nk);
            
%                 if ~ismember(ib, obj.blk.oblocks)
%                     ib
%                     if obj.blk.x{obj.ib}(1,1) < 0.2; iq = 1; else iq = size(obj.blk.x{obj.ib},1); end
%                     if obj.blk.y{obj.ib}(1,1) < 0; jq = 1; else jq = size(obj.blk.x{obj.ib},2); end
%                     fic = (obj.blk.x{obj.ib}(:,jq) - obj.blk.x{obj.ib}(1,jq))/(obj.blk.x{obj.ib}(end,jq) - obj.blk.x{obj.ib}(1,jq));
%                     fjc = (obj.blk.y{obj.ib}(iq,:) - obj.blk.y{obj.ib}(iq,:))/(obj.blk.y{obj.ib}(iq,end) - obj.blk.y{obj.ib}(iq,1));
%                     fic = fic';
% 
%                     
% 
%                     if newblk{obj.ib}.x(1,1) < 0.2; iq = 1; else iq = size(newblk{obj.ib}.x,1); end
%                     if newblk{obj.ib}.y(1,1) < 0; jq = 1; else jq = size(newblk{obj.ib}.y,2); end
%                     fi = (newblk{obj.ib}.x(:,jq) - newblk{obj.ib}.x(1,jq))/(newblk{obj.ib}.x(end,jq) - newblk{obj.ib}.x(1,jq));
%                     fj = (newblk{obj.ib}.y(iq,:) - newblk{obj.ib}.y(iq,1))/(newblk{obj.ib}.y(iq,end) - newblk{obj.ib}.y(iq,1));
%                     fi = fi';
%                 else
%                     fic = linspace(0,1,size(obj.blk.x{obj.ib},1));
%                     fjc = linspace(0,1,size(obj.blk.x{obj.ib},2));
%                     fi = linspace(0,1,size(newblk{obj.ib}.x,1));
%                     fj = linspace(0,1,size(newblk{obj.ib}.x,2));
%                 end

            [Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
            [J,I,K] = meshgrid(fj,fi,fk);

            if ndims(I) == 3
                newFlow.ro = interp3(Jc,Ic,Kc,obj.ro,J,I,K);
                newFlow.u = interp3(Jc,Ic,Kc,obj.u,J,I,K);
                newFlow.v = interp3(Jc,Ic,Kc,obj.v,J,I,K);
                newFlow.w = interp3(Jc,Ic,Kc,obj.w,J,I,K);
                newFlow.Et = interp3(Jc,Ic,Kc,obj.Et,J,I,K);
            else
                newFlow.ro = interp2(Jc,Ic,obj.ro,J,I);
                newFlow.u = interp2(Jc,Ic,obj.u,J,I);
                newFlow.v = interp2(Jc,Ic,obj.v,J,I);
                newFlow.w = interp2(Jc,Ic,obj.w,J,I);
                newFlow.Et = interp2(Jc,Ic,obj.Et,J,I);
            end
            newFlow.ib = obj.ib;
            

        end

        function writeFlow(obj, path, casetype)
            if nargin < 3
                casetype = 'cpu';
            end
            fprintf('Writing flow in block %d\n',obj.ib)

            ronow = obj.ro;
            runow = ronow.*obj.u;
            rvnow = ronow.*obj.v;
            rwnow = ronow.*obj.w;
            Etnow = obj.Et;

            fpos = 0;
            switch casetype
                case 'cpu'
                    for i=1:size(obj.ro,1)
                        for j=1:size(obj.ro,2)
                            for k=1:size(obj.ro,3)
                                fpos = fpos+1;
                                B(1, fpos) = i;
                                B(2, fpos) = j;
                                B(3,fpos) = k;
        
                                A(1,fpos) = ronow(i,j,k);
                                A(2,fpos) = runow(i,j,k);
                                A(3,fpos) = rvnow(i,j,k);
                                A(4,fpos) = rwnow(i,j,k);
                                A(5,fpos) = Etnow(i,j,k);
                            end
                        end
                    end
        
                    flopath = fullfile(path,  ['flo2_' num2str(obj.ib)]);
                    flofile = fopen(flopath,'w');
                    nodfile = fopen(fullfile(path, ['nod2_' num2str(obj.ib)]),'w');
                    %viscpath = fullfile(casedir,  ['visc_' num2str(nb)]);
                    %viscfile = fopen(viscpath,'r');
        
                    fwrite(flofile,A,'float64');
                    fwrite(nodfile,B,'uint32');
        
            
                    fclose(flofile);
                    fclose(nodfile);

                case 'gpu'
                    for k=1:size(obj.ro,3)
                        for j=1:size(obj.ro,2)
                            for i=1:size(obj.ro,1)
                                fpos = fpos+1;
                                A(1,fpos) = ronow(i,j,k);
                                A(2,fpos) = runow(i,j,k);
                                A(3,fpos) = rvnow(i,j,k);
                                A(4,fpos) = rwnow(i,j,k);
                                A(5,fpos) = Etnow(i,j,k);
                            end
                        end
                    end
        
                    flopath = fullfile(path,  ['flow_' num2str(obj.ib)]);
                    flofile = fopen(flopath,'w');
                    fwrite(flofile,A,'float64');
                    fclose(flofile);
            end
        end

        
    end
    methods (Static)
        function [fi, fj] = get_spacing(x,y)
            iq = ceil(size(x,1)/2);
            xj = x(iq,:);
            yj = y(iq,:);

            jq = ceil(size(x,2)/2);
            xi = x(:,jq);
            yi = y(:,jq);

            si = path_length(xi,yi);
            sj = path_length(xj,yj);

            fi = si/si(end);
            fj = sj/sj(end);
        end
    end

end