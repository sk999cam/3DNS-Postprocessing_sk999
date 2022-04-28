classdef meanSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        time;
        nMean;
        meanTime;
        Pr;             % Turbulence production
        diss;
    end

    properties (Dependent = true)
        
    end

    methods
        function obj = meanSlice(casedir, blk, gas)
            blk.inlet_blocks{1}
            obj@aveSlice(blk, gas);
            disp('Constructing meanSlice')
            nstats = 17;

            if nargin > 0

                fullfile(casedir, 'mean_flo', 'mean_time.txt')
                fid = fopen(fullfile(casedir, 'mean_flo', 'mean_time.txt'));
                while ~feof(fid) % Use lastest mean files
                    temp=fgetl(fid);
                end
%                 temp = fgetl(fid);
                fclose(fid);
                temp = str2num(temp);
                obj.nMean = temp(1);
                obj.meanTime = temp(3);
                for nb = 1:obj.NB

                    flopath = fullfile(casedir, 'mean_flo',  ['mean2_' num2str(nb) '_' num2str(obj.nMean)]);
                    flofile = fopen(flopath,'r');
                    fullfile(casedir, 'mean_flo', ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]);
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

                    sz = size(ro2);
                    icount(1:prod(sz)) = 0;
                    for n=1:size(A,2)

                        i = B(1,n);
                        j = B(2,n);

                        if icount(sub2ind(sz,i,j)) == 0

                            icount(sub2ind(sz,i,j)) = 1;

                            rodt(i,j) = A(1,n);
                            rudt(i,j) = A(2,n);
                            rvdt(i,j) = A(3,n);
                            rwdt(i,j) = A(4,n);
                            Etdt(i,j) = A(5,n);
    
                            ro2(i,j) = A(6,n)/obj.meanTime;
                            rou2(i,j) = A(7,n)/obj.meanTime;
                            rov2(i,j) = A(8,n)/obj.meanTime;
                            row2(i,j) = A(9,n)/obj.meanTime;
    
                            rouv(i,j) = A(10,n)/obj.meanTime;
                            rouw(i,j) = A(11,n)/obj.meanTime;
                            rovw(i,j) = A(12,n)/obj.meanTime;
    
                            p2(i,j) = A(13,n)/obj.meanTime;
                            p(i,j) = A(14,n)/obj.meanTime;
                            T(i,j) = A(15,n)/obj.meanTime;
                            ros(i,j) = A(16,n)/obj.meanTime;
                            diss(i,j) = A(17,n)/obj.meanTime;
                        end
                    end
    
                    obj.ro{nb} = rodt/obj.meanTime;
                    obj.u{nb} = rudt./(rodt);
                    obj.v{nb} = rvdt./(rodt);
                    obj.w{nb} = rwdt./(rodt);
                    obj.Et{nb} = Etdt/obj.meanTime;
                    obj.diss{nb} = diss;

                    [DUDX,DUDY] = gradHO(blk.x{nb},blk.y{nb},obj.u{nb});
                    [DVDX,DVDY] = gradHO(blk.x{nb},blk.y{nb},obj.v{nb});

                    UdUd = rou2./obj.ro{nb} - obj.u{nb}.*obj.u{nb};
                    UdVd = rouv./obj.ro{nb} - obj.u{nb}.*obj.v{nb};
                    VdVd = rov2./obj.ro{nb} - obj.v{nb}.*obj.v{nb};

                    obj.Pr{nb} = -(UdUd.*DUDX + UdVd.*(DUDY+DVDX) + VdVd.*DVDY);
                    obj.diss{nb} = diss;
                end
                
            end
            obj.getBCs(blk.inlet_blocks{1});
%             Mnow = obj.M;
%             Unow = obj.vel;
%             ronow = obj.ro;
%             munow = obj.mu;
%             %Mnow = Mnow{blk.inlet_blocks{1}};
%             pnow = obj.p;
%             %pnow = pnow{blk.inlet_blocks{1}};
%             
%             p0 = [];
%             Uinf = [];
%             muinf = [];
%             roinf = [];
%             for i=1:length(blk.inlet_blocks{1})
%                 p0now = pnow{blk.inlet_blocks{1}(i)}.*(1+((obj.gas.gam - 1)/2)*Mnow{blk.inlet_blocks{1}(i)}.^2).^(obj.gas.gam/(obj.gas.gam-1));
%                 p0 = [p0 p0now(40:100,:)];
%                 Uinf = [Uinf Unow{blk.inlet_blocks{1}(i)}(40:100,:)];
%                 muinf = [muinf munow{blk.inlet_blocks{1}(i)}(40:100,:)];
%                 roinf = [roinf ronow{blk.inlet_blocks{1}(i)}(40:100,:)];
%             end
%             obj.p0in = mean(p0,'all');
%             obj.Uinf = mean(Uinf,'all');
%             obj.muinf = mean(muinf,'all');
%             obj.roinf = mean(roinf,'all');
        end

%         function value = get.dsdy(obj)
%             s = obj.oGridProp('s');
%             value = (s(:,2:end)-s(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
%         end
%         
% 
%         function value = get.Msurf(obj)
%             disp('Calculating surface M')
%             psurf = [];
%             pnow = obj.p;
%             size(pnow{4})
%             for i=1:length(obj.oblocks)
%                 clear temp
%                 temp = pnow{obj.oblocks(i)}(:,end);
%                 size(temp)
%                 if obj.oblocks_flip(i) == 1
%                     temp = flip(temp);
%                 end
%                 psurf = [psurf temp'];
%             end
%             value = sqrt((2/(obj.gas.gam - 1)) * ( (psurf/obj.p0in).^(-(obj.gas.gam-1)/obj.gas.gam) - 1));
%         end
% 
%         function value = get.BLedgeInd(obj)
%             temp = obj.dsdy;
%             for i=1:size(temp,1)
%                 j=3;
%                 temp(i,j);
%                 while temp(i,j) < obj.dsdyThresh
%                     j = j+1;
%                 end
%                 value(i) = j+1;
%             end
%         end
% 
%         function value = get.delta99(obj)
%             inds = obj.BLedgeInd;
%             for i=1:size(obj.yBL,1)
%                 value(i) = obj.yBL(i,inds(i));
%             end
%         end
% 
%         function value = get.U(obj)
%             disp('Calculating U')
%             unow = obj.oGridProp('u');
%             vnow = obj.oGridProp('v');
%             nnow = obj.n;
%             value = zeros(size(obj.yBL));
%             for i=1:size(obj.yBL,1)
%                 for j=1:size(obj.yBL,2)
%                     velnow = [unow(i,j); vnow(i,j)];
%                     value(i,j) = norm(velnow - nnow(:,i)*dot(nnow(:,i),velnow));
%                 end
%             end
%         end
% 
%         function value = get.delStar(obj)
%             inds = obj.BLedgeInd;
%             ronow = obj.oGridProp('ro');
%             Unow = obj.U;
%             value = zeros(1,length(inds));
%             for i=1:size(obj.yBL,1)
%                 integrand = 1 - ronow(i,1:inds(i)).*Unow(i,1:inds(i))/(ronow(i,inds(i))*Unow(i,inds(i)));
%                 ys = obj.yBL(1:inds(i));
%                 value(i) = trapz(ys, integrand);
%             end
%         end
% 
%         function value = get.Res(obj)
%             value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
%         end

        

        
    end
end