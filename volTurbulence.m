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
    end

    properties (Dependent = true)
        Q;
    end

    methods
        function obj = volTurbulence(path, ilength, nk, lturb, y, span)

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
                nk = length(z);
                x = linspace(0,(ilength-1)*lturb,ilength);
    
                obj.u = zeros(obj.ni,obj.nj,nk);
                obj.v = zeros(obj.ni,obj.nj,nk);
                obj.w = zeros(obj.ni,obj.nj,nk);
                obj.blk.x{1} = zeros(obj.ni,obj.nj);
                obj.blk.y{1} = zeros(obj.ni,obj.nj);
                obj.blk.z = linspace(0,span,nk);
                n=1;
                for i = 1:obj.ni
                    fprintf('i=%d/%d\n',i,obj.ni)
                    for j = 1:obj.nj
                        for k = 1:obj.nk
                            obj.u(i,j,k) = q(1,n);
                            obj.v(i,j,k) = q(2,n);
                            obj.w(i,j,k) = q(3,n);
                            n = n+1;
                            obj.blk.x{1}(i,j) = x(i);
                            obj.blk.y{1}(i,j) = y(j);
                        end
                    end
                end
            else

                

            end
        end

        function value = slice(obj,k)
            if nargin == 1
                k = floor(obj.nk/2);
            end

            value = flowSlice(obj.blk, obj.gas);
            value.u{1} = obj.u{nb}(:,:,k);
            value.v{1} = obj.v{nb}(:,:,k);
            value.w{1} = obj.w{nb}(:,:,k);
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
            value = Q_criterion(obj.blk.x{1},obj.blk.y{1},obj.blk.z,obj.u,obj.v,obj.w);
        end

        function plot_Q_criterion(obj, thresh)
%             X = repmat(obj.blk.x{1},[1 1 obj.nk]);
%             Y = repmat(obj.blk.y{1},[1 1 obj.nk]);
%             Z = repmat(obj.blk.z,[size(obj.blk.x{1}) 1]);
            X = obj.blk.x{1}(:,1)';
            Y = obj.blk.y{1}(1,:);
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
            ilength_new = floor(obj.blk.x{1}(end,1)/(lturb_new))+1;

            y_new = new_case.y_inlet;
            fi = linspace(0,1,ilength_new);
            fj = (y_new - y_new(1))/(y_new(end)-y_new(1));
            fk = linspace(0,1,new_case.blk.nk{1});
            
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
    end
end