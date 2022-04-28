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

    methods
        function obj = volTurbulence(nj_inlet, nk, inlet_width, path)
            obj.nj = nj_inlet
            obj.nk = nk
            obj.inlet_width = inlet_width;
            path = fullfile(path, 'inflow_turb.dat');
            turbfile = fopen(path,'r');
            q = fread(turbfile,inf,'float64');
            q = reshape(q,3,length(q)/3);
            size(q)
            fclose(turbfile);
            obj.ni = size(q,2)/nk/nj_inlet
            ilength = inlet_width*(obj.ni-1)/(nj_inlet-1);
            x = linspace(ilength,0,obj.ni);
            y = linspace(-inlet_width/2,inlet_width/2,nj_inlet);
            obj.u = zeros(obj.ni,nj_inlet,nk);
            obj.v = zeros(obj.ni,nj_inlet,nk);
            obj.w = zeros(obj.ni,nj_inlet,nk);
            obj.blk.x{1} = zeros(obj.ni,nj_inlet);
            obj.blk.y{1} = zeros(obj.ni,nj_inlet);
            n=1;
            for i = 1:obj.ni
                fprintf('i=%d/%d\n',i,obj.ni)
                for k = 1:obj.nk
                    for j = 1:obj.nj
                        obj.u(i,j,k) = q(1,n);
                        obj.v(i,j,k) = q(2,n);
                        obj.w(i,j,k) = q(3,n);
                        n = n+1;
                        obj.blk.x{1}(i,j) = x(i);
                        obj.blk.y{1}(i,j) = y(j);
                    end
                end
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
    end
end