classdef flowSlice < handle
    %FLOWSLICE Generic class containing a slice of the flow

    properties
        NB;
        gas;
        ro;
        u;
        v;
        w;
        Et;
%         xSurf;
%         yBL;
%         xO;
%         yO;
%         iO;
%         jO;
%         blkO;
%         oblocks;
%         oblocks_flip;
        blk;
%         iLE;
%         iTE;
%         n;
%         ssurf;          % Surface distance fron LE
%         vortZ;          % Z vorticity
    end

    properties (Abstract)
    end

    properties (Dependent = true)
        T;             % Temperature
        p;             % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        mu;             % Viscosity
        nu;             % Kinematic viscosity
        p0;
        schlieren;      % |grad(ro)|/ro
        cellSize;
        S_an_mag;   % Magnitude of Anisotripic part of strain rate tensor
    end

    methods (Abstract)
        plot
    end

    methods
        function obj = flowSlice(blk, gas)
            %FLOWSLICE Construct a flowSlice object
            disp('Constructing flowSlice')

            obj.gas = gas;
            obj.NB = size(blk.blockdims,1);
            obj.blk = blk;
        end         % End of constructor

        function value = get.p(obj)
            value = obj.get_p;
        end

        function value = get_p(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            value = obj.get_T;
        end

        function value = get_T(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.p{nb}./(obj.ro{nb}*obj.gas.rgas);
            end
        end

        function value = get.vel(obj)
            disp('Calculating vel')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.p0(obj)
            pnow = obj.p;
            Mnow = obj.M;
            for nb = 1:obj.NB
                value{nb} = pnow{nb}.*(1+0.5*(obj.gas.gam-1)*Mnow{nb}.^2).^(obj.gas.gam/(obj.gas.gam-1));
            end
        end

        function value = get.M(obj)
            disp('Calculating M')
            value = cell(1,obj.NB);
            pnow = obj.p;
            for nb = 1:obj.NB
                %pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow{nb}./(obj.ro{nb}*obj.gas.rgas);
                velnow = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
                value{nb} = velnow./sqrt(obj.gas.gam*obj.gas.rgas*Tnow);
            end
        end

        function value = get.s(obj)
            disp('Calculating s')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.cp*log(Tnow/300) - obj.gas.rgas*log(pnow/1e5);
            end
        end

        function value = get.mu(obj)
            disp('Calcualting mu')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = obj.p{nb};%(obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.mu_ref*(Tnow/obj.gas.mu_tref).^(3/2) .* (obj.gas.mu_cref + obj.gas.mu_tref)./(obj.gas.mu_cref + Tnow);
            end
        end

        function value = get.nu(obj)
            munow = obj.mu;
            ronow = obj.ro;
            for nb = 1:obj.NB
                value{nb} = munow{nb}./ronow{nb};
            end
        end

        

        function value = get.cellSize(obj)
            fprintf('Calculating Cell Sizes\n')
            dz = obj.blk.span/(obj.blk.nk{1}-1);
            
            value = {};
            for ib = 1:obj.NB
                ni = size(obj.blk.x{ib},1);
                nj = size(obj.blk.x{ib},2);
                area = zeros(ni-1, nj-1);
                for i=1:ni-1
                    for j=1:nj-1
                        xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
                            obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
                        ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
                            obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];
                        area(i,j) = abs(polyarea(xnow,ynow));
                    end
                end
                area = dz*area;
                value{ib}(1,1) = area(1,1);
                value{ib}(1,nj) = area(1,nj-1);
                value{ib}(ni,1) = area(ni-1,1);
                value{ib}(ni,nj) = area(ni-1,nj-1);
                for i = 2:ni-1
                    value{ib}(i,1) = 0.5*(area(i-1,1)+area(i,1));
                    value{ib}(i,end) = 0.5*(area(i-1,end)+area(i,end));
                end
                for j = 2:nj-1
                    value{ib}(1,j) = 0.5*(area(1,j-1)+area(1,j));
                    value{ib}(end,j) = 0.5*(area(end,j-1)+area(end,j));
                end
                for i=2:ni-1
                    for j=2:nj-1
                        value{ib}(i,j) = 0.25*(area(i-1,j-1)+area(i-1,j)+area(i,j-1) +area(i,j));
                    end
                end
                value{ib} = value{ib}.^(1/3);
            end
        end


        function value = get.schlieren(obj)
            disp('calculating grad(ro)/ro')
            value = cell(1,obj.NB);
            for nb=1:obj.NB
                [drodx, drody] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.ro{nb});
                value{nb} = sqrt(drodx.^2 + drody.^2)./obj.ro{nb};
            end
        end

        function getSize(obj)
            props = properties(obj);
            totSize = 0; 
            for ii=1:length(props) 
                currentProperty = obj.(props{ii});
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end

        

        function value = get.S_an_mag(obj)
            for ib = 1:obj.NB
                [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
                [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});
        
                S = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
                
                
                S(:,:,1,1) = 2*DUDX/3 - DVDY/3;
                S(:,:,2,2) = 2*DVDY/3 - DUDX/3;
                S(:,:,3,3) = -(DUDX+DVDY)/3;
        
                S(:,:,1,2) = 0.5*(DUDY+DVDX);
                S(:,:,2,1) = S(:,:,1,2);

                value{ib} = sqrt(abs(sum(sum(S.*S,4),3)));
            end

        end
            
    end
end
