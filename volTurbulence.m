classdef volTurbulence < handle
    % VOLTURBULENCE Contains inflow turbulence fluctions for DNS case

    properties
        ni;
        nj;
        nk;
        inlet_width;
        blk;
        u;
        v;
        w;
        x;
        y;
        casepath;
    end

    properties (Dependent = true)
        Q;
        k;
        lturb;
    end

    methods
        function obj = volTurbulence(path, y, ilength, nk, lturb, span)

            if nargin == 6
                
                path = fullfile(path, 'inflow_turb.dat');
                turbfile = fopen(path,'r');
                q = fread(turbfile,inf,'float64');
                q = reshape(q,3,length(q)/3);
                size(q)
                fclose(turbfile);
                obj.ni = ilength;
                obj.nj = size(q,2)/nk/ilength;
                obj.nk = nk;
    
                if length(y) == 2
                    y = linspace(y(1), y(end), obj.nj);
                end
                
                obj.nk = nk;
                obj.inlet_width = abs(y(end) - y(1));
                x = linspace(0,(ilength-1)*lturb,ilength);
    
                obj.u = zeros(obj.ni,obj.nj,nk);
                obj.v = zeros(obj.ni,obj.nj,nk);
                obj.w = zeros(obj.ni,obj.nj,nk);
                obj.blk.x = zeros(obj.ni,obj.nj);
                obj.blk.y = zeros(obj.ni,obj.nj);
                obj.blk.z = linspace(0,span,nk);
                obj.blk.blockdims = [obj.ni obj.nj, obj.nk];
%                 n=1;
%                 for i = 1:obj.ni
%                     fprintf('i=%d/%d\n',i,obj.ni)
%                     for j = 1:obj.nj
%                         for k = 1:obj.nk
%                             obj.u(i,j,k) = q(1,n);
%                             obj.v(i,j,k) = q(2,n);
%                             obj.w(i,j,k) = q(3,n);
%                             n = n+1;
%                             obj.blk.x(i,j) = x(i);
%                             obj.blk.y(i,j) = y(j);
%                         end
%                     end
%                 end

                obj.blk.x = repmat(x', [1 obj.nj]);
                obj.blk.y = repmat(y, [obj.ni 1]);

                fid_f = fopen(path,'r'); %
                A = fread(fid_f,inf,'float64');
                A = reshape(A,3,obj.nj,nk,obj.ni);
                
                
                u = squeeze(A(1,:,:,:));
                v = squeeze(A(2,:,:,:));
                w = squeeze(A(3,:,:,:));
                u = permute(u,[3,1,2]);
                v = permute(v,[3,1,2]);
                w = permute(w,[3,1,2]);

                obj.u = u;
                obj.v = v;
                obj.w = w;

            elseif nargin == 1
                if exist(fullfile(path, 'inflow_turb.dat'),'file')
                    turbfile = fopen(fullfile(path, 'inflow_turb.dat'),'r');
                elseif exist(fullfile(path,'turbid_files','inflow_turb.dat'),'file')
                    obj.casepath = path;
                    turbfile = fopen(fullfile(path,'turbid_files','inflow_turb.dat'),'r');
                else
                    fprintf('inflow_turb.dat not found\n');
                    return
                end

                if exist(fullfile(path, 'turbid.txt'),'file')
                    inputfile = fopen(fullfile(path, 'turbid.txt'),'r');
                elseif exist(fullfile(path,'turbid_files','turbid.txt'),'file')
                    inputfile = fopen(fullfile(path,'turbid_files','turbid.txt'),'r');
                else
                    fprintf('turbid.txt not found\n');
                    return
                end

                temp = str2num(char(split(fgetl(inputfile))));
                nx = temp(1); ny = temp(2); nz = temp(3);
                temp = str2num(char(split(fgetl(inputfile))));
                lsx = temp(1); lsy = temp(2); lsx = temp(3);
                temp = str2num(char(split(fgetl(inputfile))));
                nfx = temp(1); nfy = temp(2); nfx = temp(3);
                fclose(inputfile);

                

                obj.ni = nx; obj.nj = ny; obj.nk = nz;
                obj.inlet_width = 1;
                obj.blk.y = repmat(linspace(-0.5, 0.5, ny), [nx 1]);
                obj.blk.x = repmat(linspace(0, (nx-1)/(ny-1), nx)', [1 ny]);
                obj.blk.z = linspace(0, (nz-1)/(ny-1), nz);
                obj.blk.blockdims = [nx ny nz];

                A = fread(turbfile,inf,'float64');
                A = reshape(A,3,ny,nz,nx);
                u = squeeze(A(1,:,:,:));
                v = squeeze(A(2,:,:,:));
                w = squeeze(A(3,:,:,:));
                obj.u = permute(u,[3,1,2]);
                obj.v = permute(v,[3,1,2]);
                obj.w = permute(w,[3,1,2]);

                fclose(turbfile);

%                 obj.u = zeros(nx,ny,nz);
%                 obj.v = zeros(nx,ny,nz);
%                 obj.w = zeros(nx,ny,nz);
% 
%                 q = fread(turbfile,inf,'float64');
%                 q = reshape(q,3,length(q)/3);
% 
%                 
%                 n = 1;
% 
%                 for i = 1:obj.ni
%                     if mod(i, 20) == 0
%                         fprintf('i=%d/%d\n',i,obj.ni)
%                     end
%                     for j = 1:obj.nj
%                         for k = 1:obj.nk
%                             obj.u(i,j,k) = q(1,n);
%                             obj.v(i,j,k) = q(2,n);
%                             obj.w(i,j,k) = q(3,n);
%                             n = n+1;
%                         end
%                     end
%                 end
                
            elseif nargin == 3 && ischar(y)
                ni = ilength(1); nj = ilength(2); nk = ilength(3);
                obj.ni = ni; obj.nj = nj; obj.nk = nk;
                gridName = y;

                fid_f = fopen(path,'r'); %
                A = fread(fid_f,inf,'float64');
                A = reshape(A,3,nj,nk,ni);
                
                
                u = squeeze(A(1,:,:,:));
                v = squeeze(A(2,:,:,:));
                w = squeeze(A(3,:,:,:));
                obj.u = permute(u,[3,1,2]);
                obj.v = permute(v,[3,1,2]);
                obj.w = permute(w,[3,1,2]);
                
                
                fid_g = fopen(gridName,'r'); %
                G = fread(fid_g,inf,'float64');       
                
                
                
                G = reshape(G,3,length(G)/3);  
                xb = reshape(G(1,:),[ni,nj,nk]);
                yb = reshape(G(2,:),[ni,nj,nk]);
                zb = reshape(G(3,:),[ni,nj,nk]);
                fclose(fid_g);
                obj.blk.x = squeeze(xb(:,:,1));
                obj.blk.y = squeeze(yb(:,:,1));
                obj.blk.z = squeeze(zb(1,1,:));
                obj.blk.blockdims = [ni nj nk];



            end
        end

        function value = slice(obj,k, gas, bcs)
            if nargin == 1 || isempty(k)
                k = floor(obj.nk/2);
            end

            value = flowSlice(obj.blk, gas, bcs);
            value.u{1} = obj.u(:,:,k);
            value.v{1} = obj.v(:,:,k);
            value.w{1} = obj.w(:,:,k);
            value.blk.x = {};
            value.blk.y = {};
            value.blk.x{1} = obj.blk.x;
            value.blk.y{1} = obj.blk.y;
        end

        function apply_window(obj, vector)
            if length(vector) ~= obj.nj
                disp('Vector incorrect length')
            else
                B = repmat(vector, obj.ni, 1, obj.nk);
                disp('processing u')
                obj.u = obj.u.*B;
                disp('processing v')
                obj.v = obj.v.*B;
                disp('processing w')
                obj.w = obj.w.*B;
            end
        end

        function write_turb(obj, path)
            disp('Writing inflow turbulence')
            path = fullfile(path, 'inflow_turb_new.dat');
            
            buf = zeros(3*obj.ni*obj.nj*obj.nk,1);
            n = 1;
            for i=1:obj.ni
                fprintf('i=%d/%d\n',i,obj.ni)
                for k=1:obj.nk
                    for j=1:obj.nj
                        buf(n)=obj.u(i,j,k);
                        n=n+1;
                        buf(n)=obj.v(i,j,k);
                        n=n+1;
                        buf(n)=obj.w(i,j,k);
                        n=n+1;
                    end
                end
            end
            turbfile = fopen(path,'w');
            fwrite(turbfile,buf,'float64');
            fclose(turbfile);
        end

        function value = get.Q(obj)
            value = Q_criterion(obj.blk.x,obj.blk.y,obj.blk.z,obj.u,obj.v,obj.w);
        end

        function value = get.k(obj)
            value = 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2);
        end

        function value = Tu(obj, vref)
            value = sqrt(2*mean(obj.k,'all')/3)/vref;
        end

        function plot_Q_criterion(obj, thresh)
%             X = repmat(obj.blk.x,[1 1 obj.nk]);
%             Y = repmat(obj.blk.y,[1 1 obj.nk]);
%             Z = repmat(obj.blk.z,[size(obj.blk.x) 1]);
            X = obj.blk.x(:,1)';
            Y = obj.blk.y(1,:);
            Z = obj.blk.z;
            [X, Y, Z] = meshgrid(X,Y,Z);
            Q = permute(obj.Q,[2 1 3]);
            isosurface(X, Y, Z, Q, thresh)
            axis equal
        end

        function [ilength_new, lturb_new] = interp_turb(obj, new_case, write)
            
            if nargin < 3
                write = false;
            end
            yc = obj.blk.y{1}(1,:);

            fic = linspace(0,1,obj.ni);
            fjc = (yc - yc(1))/(yc(end)-yc(1));
            fkc = linspace(0,1,obj.nk);

            lturb_new = new_case.inlet_width/(new_case.nj_inlet-1);
            ilength_new = floor(obj.blk.x(end,1)/(lturb_new))+1;

            y_new = new_case.y_inlet;
            fi = linspace(0,1,ilength_new);
            fj = (y_new - y_new(1))/(y_new(end)-y_new(1));
            fk = linspace(0,1,new_case.blk.nk);
            
                
            %[Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
            [J,I,K] = meshgrid(fj,fi,fk);
            fprintf('Interpolating u\n')
            u_new = interp3(fjc,fic,fkc,obj.u,J,I,K);
            fprintf('Inerpolating v\n')
            v_new = interp3(fjc,fic,fkc,obj.v,J,I,K);
            fprintf('Interpolating w\n')
            w_new = interp3(fjc,fic,fkc,obj.w,J,I,K);

            if write
                disp('Writing inflow turbulence')
                path = fullfile(new_case.casepath, 'inflow_turb_new.dat');
                
                buf = zeros(3*obj.ni*obj.nj*obj.nk,1);
                n = 1;
                clear Ic Jc Kc I J K
    
                for i=1:length(fi)
                    fprintf('i=%d/%d\n',i,length(fi))
                    for k=1:length(fk)
                        for j=1:length(fj)
                            buf(n)=u_new(i,j,k);
                            n=n+1;
                            buf(n)=v_new(i,j,k);
                            n=n+1;
                            buf(n)=w_new(i,j,k);
                            n=n+1;
                        end
                    end
                end

            

                clear u_new v_new w_new
                turbfile = fopen(path,'w');
                fwrite(turbfile,buf,'float64');
                fclose(turbfile);
            end
        end

        function plot_mesh(obj, skip, ax)
            if nargin<2 || isempty(skip)
                skip=8;
            end

            if nargin<3 || isempty(ax)
                ax = gca;
            end
            
            hold on
            [ni, nj] = size(obj.blk.x);
            for j=[1:skip:ni ni]
                if (j==1) || (j==ni)
                    plot(ax, obj.blk.x(j,:),obj.blk.y(j,:),'r','LineWidth',1)
                else
                    plot(ax, obj.blk.x(j,:),obj.blk.y(j,:),'k')
                end
            end
            for j=[1:skip:nj nj]
                if (j==1) || (j==nj)
                    plot(ax, obj.blk.x(:,j),obj.blk.y(:,j),'r','LineWidth',1)
                else
                    plot(ax, obj.blk.x(:,j),obj.blk.y(:,j),'k')
                end
            end
        
            
            axis equal
        end

        function value = get.lturb(obj)
            value = (obj.blk.x(end,1) - obj.blk.x(1,1))/(size(obj.blk.x,1)-1);
        end

        function value = get_aturb(obj,Tu,vref)
            value = Tu/obj.Tu(vref);
        end
    end
end