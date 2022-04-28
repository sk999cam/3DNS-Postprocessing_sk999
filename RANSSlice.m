classdef RANSSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        blk;
        StR;            % Strain rate
        mut;            % Eddy viscosity
        k;
        omega;
        Reth;
        gamma;
        trans;          % Transition model on/off
        mod;            % beta1 modification on/off
    end

    properties (Dependent = true)
        Pr;             % Turbulence production
        Pr_dist;
    end

    methods
        function obj = RANSSlice(basedir, trans, mod, blk, gas, nodes)
            obj@aveSlice(blk, gas);
            disp('Constructing RANSSlice')
            obj.blk = blk;
            obj.trans = trans;
            obj.mod = mod;
            transstr = ["turb","trans"];
            modstr = ["bsl","mod"];
            data = readtable(fullfile(basedir,transstr(trans+1) + "_" + modstr(mod+1),'hf2d.txt'));

            if nargin > 0
                
                cv = gas.cp/gas.gam;
                for nb=1:obj.NB
                    fprintf('nb: %d\n', nb)
                    nib = blk.blockdims(nb,1);
                    njb = blk.blockdims(nb,2);
                    ro = zeros(nib,njb);
                    u = zeros(nib,njb);
                    v = zeros(nib,njb);
                    w = zeros(nib,njb);
                    Et = zeros(nib,njb);
                    k = zeros(nib,njb);
                    omega = zeros(nib,njb);
                    mut = zeros(nib,njb);
                    StR = zeros(nib,njb);
                    if obj.trans == 1
                        gamma = zeros(nib,njb);
                        Reth = zeros(nib,njb);
                    end

                    nodesnow = nodes{nb};
                    for i=1:nib
                        for j=1:njb
                            node = nodesnow(i,j);
                            rotmp = data.density(node);
                            utmp = data.x_velocity(node);
                            vtmp = data.y_velocity(node);
                            Ttmp = data.temperature(node);
                            ro(i,j) = rotmp;
                            u(i,j) = utmp;
                            v(i,j) = vtmp;
                            Et(i,j) = rotmp*(cv*Ttmp + 0.5*(utmp^2+vtmp^2));
                            k(i,j) = data.turb_kinetic_energy(node);
                            omega(i,j) = data.specific_diss_rate(node);
                            mut(i,j) = data.viscosity_turb(node);
                            StR(i,j) = data.strain_rate_mag(node);
                            if obj.trans == 1
                                gamma(i,j) = data.intermittency(node);
                                Reth(i,j) = data.momentum_thickness_re(node);
                            end
                        end
                    end
                    obj.ro{nb} = ro;
                    obj.u{nb} = u;
                    obj.v{nb} = v;
                    obj.w{nb} = w;
                    obj.Et{nb} = Et;
                    obj.k{nb} = k;
                    obj.omega{nb} = omega;
                    obj.mut{nb} = mut;
                    obj.StR{nb} = StR;
                    if obj.trans == 1
                        obj.gamma{nb} = gamma;
                        obj.Reth{nb} = Reth;
                    end
                end
                %clear nodes
            end
            obj.getBCs(blk.inlet_blocks{1});
        end

        function kPlot(obj,prop,ax,lims)
            
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis([-0.2 2 -0.5 0.5])
            axis equal
            cb = colorbar;
            if nargin == 5
                caxis(lims)
            end
        end

        function value = get.Pr(obj)
            disp('Calculating Pr')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = obj.mut{nb}.*obj.StR{nb}.^2;
            end
        end

        
    end
end