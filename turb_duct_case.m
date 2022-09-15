classdef turb_duct_case < DNS_case
    %TURB_DUCT_CASE Inflow turbulence duct case.
    % Derived from DNS_case

    properties
        
    end

    methods
        function obj = turb_duct_case(casename,run)
            %TURB_DUCT_CASE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                run = [];
            end
            obj@DNS_case(casename,run);
        end

        function turb = case2volTurb(obj,xLength)
            %CASE2VOLTURB Create volTurbulence object,
            % taking downstrem xLength of domain.

            if nargin < 2
                xLength = obj.inlet_width;
            end

            if isempty(obj.instFlow)
                obj.readInstFlow;
            end

            turb = volTurbulence();
            xCut = obj.blk.x{1}(end,1) - xLength;
            [~, iCut] = min(abs(obj.blk.x{1}(:,1) - xCut));
            xCut = obj.blk.x{1}(iCut,1);

            
            turb.blk.z = linspace(0,obj.solver.span,obj.blk.nk{1});
            
            turb.ni = size(obj.blk.x{1},1)+1-iCut;
            turb.nj = length(obj.y_inlet);
            turb.nk = obj.blk.nk{1};
            turb.inlet_width = obj.inlet_width;
                
            utmp = zeros(turb.ni,1,turb.nk);
            vtmp = zeros(turb.ni,1,turb.nk);
            wtmp = zeros(turb.ni,1,turb.nk);
            xtmp = zeros(turb.ni,1);
            ytmp = zeros(turb.ni,1);

            for ii = 1:length(obj.blk.inlet_blocks{1})
                ib = obj.blk.inlet_blocks{1}(ii);
                utmp = [utmp(:,1:end-1,:) obj.instFlow.u{ib}(iCut:end,:,:)];
                vtmp = [vtmp(:,1:end-1,:) obj.instFlow.v{ib}(iCut:end,:,:)];
                wtmp = [wtmp(:,1:end-1,:) obj.instFlow.w{ib}(iCut:end,:,:)];
                xtmp = [xtmp(:,1:end-1) obj.blk.x{ib}(iCut:end,:)];
                ytmp = [ytmp(:,1:end-1) obj.blk.y{ib}(iCut:end,:)];
            end

            turb.blk.x{1} = xtmp-xCut;
            turb.blk.y{1} = ytmp;
            turb.u = obj.normalize(utmp);
            turb.v = obj.normalize(vtmp);
            turb.w = obj.normalize(wtmp);

            

        end

    end
    methods (Static)
        function qhat = normalize(q)
            qbar = mean(q,1);
            q2bar = mean(q.^2,1);
            qrms = sqrt(q2bar - qbar.^2);

            qhat = (q-qbar)./qrms;
        end
    end
end