function write_input_files(casename,blk,bcs,gas,solver,trip,varargin)
nargin;
p = inputParser;
addOptional(p,'topology',[]);
addOptional(p,'nkproc',[]);
addOptional(p,'casetype','both');
parse(p,varargin{:});
casetype = solver.version;
nkproc = p.Results.nkproc;
topology = p.Results.topology;
if isempty(nkproc)
    nkproc = ceil(solver.nk/solver.npp);
    fprintf('k procs: %d\n', nkproc)
end
if isempty(topology)
    nb = length(blk.x);
    switch nb
        case 1
            topology = 3;
        case 9
            topology = 1;
        case 12
            topology = 2;
    end
end

if ~isfield(gas,'gamma')
    gas.gamma = gas.gam;
end

npp = solver.npp;
fprintf('Using npp = %d\n', npp)
NB = length(blk.x);
ncorner = length(blk.corner);   

dir = fullfile(pwd,casename);
fprintf('Writing input files to directory: %s\n',dir)

if(~exist(dir,'dir'))
mkdir(dir);
end
if strcmp(casetype,'gpu')
    fid = fopen(fullfile(dir,'body.txt'),'a');
    fclose(fid);
end
if ismember(casetype, {'gpu', 'all'})
    % Now write header for new input file
    % GPU input
    fidin = fopen(fullfile(dir,'input_gpu.txt'),'w');
    fprintf(fidin,'%d %d\n', [NB,1]);
       
    nprocs = 0;
    
    for ib=1:NB
        
    x = blk.x{ib};
    y = blk.y{ib};
    [ni,nj]=size(x);   
    
    im_next_block = blk.next_block{ib}.im;
    ip_next_block = blk.next_block{ib}.ip;
    jm_next_block = blk.next_block{ib}.jm;
    jp_next_block = blk.next_block{ib}.jp;
    
    im_next_patch = blk.next_patch{ib}.im;
    ip_next_patch = blk.next_patch{ib}.ip;
    jm_next_patch = blk.next_patch{ib}.jm;
    jp_next_patch = blk.next_patch{ib}.jp;
        
    
    %nkproc = 1;
    %niproc = 1;%max([1,floor(ni/(4*npp))]);
    %njproc = 1;%max([1,floor(nj/(4*npp))]);
    %nprocs = nprocs + niproc*njproc*nkproc;
    
    fprintf(fidin,'%d %d %d\n', [ni nj solver.nk]);
      
    %fprintf(fidin,'%d %d %d\n', [niproc njproc nkproc]);    
       
          
       im = (im_next_block==0)*im_next_patch;
       ip = (ip_next_block==0)*ip_next_patch;
       
       jm = (jm_next_block==0)*jm_next_patch;
       jp = (jp_next_block==0)*jp_next_patch;
       
       
       fprintf(fidin,'%d %d %d %d\n', [im ip jm jp]);
       
       if(im==0)
       fprintf(fidin,'%d %d\n', [im_next_block im_next_patch]);
       end
       
       if(ip==0)
       fprintf(fidin,'%d %d\n', [ip_next_block ip_next_patch]);
       end
       
       if(jm==0)
       fprintf(fidin,'%d %d\n', [jm_next_block jm_next_patch]);
       end
       
       if(jp==0)
       fprintf(fidin,'%d %d\n', [jp_next_block jp_next_patch]);
       end
       
       NI(ib) = ni;
       NJ(ib) = nj;
       
    end
    
    % write corners
    fprintf(fidin,'%d\n', ncorner);
    for n=1:ncorner
        cor_type(n) = 0;
    end
    
    for n=1:ncorner
    fprintf(fidin,'%d %d\n', [blk.corner{n}.Nb, cor_type(n)]);
    for nb=1:blk.corner{n}.Nb
    ib = blk.corner{n}.block{nb};
    ic = blk.corner{n}.i{nb};
    jc = blk.corner{n}.j{nb};
    if(ic>1); ic = NI(ib); end
    if(jc>1); jc = NJ(ib); end
    fprintf(fidin,'%d %d %d\n', [ib ic jc]);
    end
    
    end
    
    fprintf(fidin,'%d\n',blk.nbg); % 1 block group
    for nbg = 1:blk.nbg
        nb_bg = length(blk.block_groups{nbg});
        fprintf(fidin,'%d\n',nb_bg); % NB blocks in group
        for ib=1:nb_bg
            fprintf(fidin,'%d',blk.block_groups{nbg}(ib)); % blocks in group
        end
    end
    
    % 
    % write rest of file
    % check these are ok for your case!
        
        %nsteps nwrite ncut 
        fprintf(fidin,'\n%d %d %d\n', [solver.niter solver.nwrite solver.ncut]);
        
        % cfl, filter coefficient
        fprintf(fidin,'%f %f\n', [solver.cfl solver.sigma]);
        
         % Toin poin pext vref alpha_in pitch_in aturb (not used) ilength (not used) g_z
        fprintf(fidin,'%f %f %f %f %f %f %f %d %d %f\n', [bcs.Toin bcs.Poin bcs.pexit bcs.vin bcs.alpha bcs.gamma bcs.aturb bcs.nturb bcs.iradprof bcs.g_z]);
        
        % gamma cp mu_ref sutherlands constants prandtl no.
        fprintf(fidin,'%f %f %12.5e %f %f %f\n', [gas.gamma gas.cp gas.mu_ref gas.mu_tref gas.mu_cref gas.pr]);
        
        % span, spanwise grid expansion factor
        fprintf(fidin,'%f %f\n', [solver.span solver.fexpan]);
        
        % restart, statistics
        fprintf(fidin,'%d %d\n', [solver.irestart solver.istats]);

        % inlet BL, theta
        fprintf(fidin,'%d %d\n', [solver.ilam bcs.theta]);
    
    fclose(fidin);

    % write bl_recycle.txt
    fblrec = fopen(fullfile(dir,'bl_recycle.txt'),'w');
    fprintf(fblrec,'%f\n', [solver.iblrec]);
    fclose(fblrec);

    % free stream buffer region
    ffreebuff = fopen(fullfile(dir,'free_buffer.txt'),'w');
    fprintf(ffreebuff,'%d \n', [solver.freebuff]);
    fclose(ffreebuff);

    % write in_buffer.txt
    finbuff = fopen(fullfile(dir,'in_buffer.txt'),'w');
    fprintf(finbuff,'%d %f\n', [5, 0.0001]);
    fclose(finbuff);

    % write out_buffer.txt
    foutbuff = fopen(fullfile(dir,'out_buffer.txt'),'w');
    fprintf(foutbuff,'%d %f\n', [10, 0.01]);
    fclose(foutbuff);

    %write freestream.txt
    ffreestream = fopen(fullfile(dir, 'freestream.txt'), 'w');
    fprintf(ffreestream,'%d', [bcs.nfsp]);
    diff = 1/(bcs.nfsp-1);
    x = 0;

    for fsp=1:bcs.nfsp
        u=((1/bcs.vin)-(bcs.k*bcs.Lref*x/gas.mu_ref))^(-1);
        
        vstarsq = (u/sqrt(gas.cp*bcs.Toin))^2;
        gm1 =gas.gamma -1;

        M = sqrt((vstarsq)/(gm1*(1-0.5*vstarsq)));
        p = bcs.Poin*(1+0.5*gm1*M^2)^(-gas.gamma/gm1);
        fprintf(ffreestream, '\n%f %f %f', [x,u,p]);
        x = x+diff;
    end
    fclose(ffreestream);
    % write trip.txt
    ftrip = fopen(fullfile(dir, 'trip.txt'), 'w');
    fprintf(ftrip, '%d', [trip.ifLEtrip]);
    fprintf(ftrip, '\n%f %f %f %f', [trip.x1,trip.y1,trip.x2,trip.y2]);
    fprintf(ftrip, '\n%f %d', [trip.tripscale,trip.kspace]);
    fclose(ftrip);
end


if ismember(casetype, {'cpu','both'})
    % Now write header for new input file
    % CPU input
    fidin = fopen(fullfile(dir,'input_cpu.txt'),'w');
    fprintf(fidin,'%d\n', [NB]);
       
    nprocs = 0;
    
    for ib=1:NB
        
    x = blk.x{ib};
    y = blk.y{ib};
    [ni,nj]=size(x);
    nk = solver.nk;
    
    im_next_block = blk.next_block{ib}.im;
    ip_next_block = blk.next_block{ib}.ip;
    jm_next_block = blk.next_block{ib}.jm;
    jp_next_block = blk.next_block{ib}.jp;
    
    im_next_patch = blk.next_patch{ib}.im;
    ip_next_patch = blk.next_patch{ib}.ip;
    jm_next_patch = blk.next_patch{ib}.jm;
    jp_next_patch = blk.next_patch{ib}.jp;
        
    
    
    niproc = ceil(ni/npp);
    njproc = ceil(nj/npp);
    nprocs = nprocs + niproc*njproc*nkproc;
    
    fprintf(fidin,'%d %d %d\n', [ni nj nk]);
      
    fprintf(fidin,'%d %d %d\n', [niproc njproc nkproc]);    
       
          
       im = (im_next_block==0)*im_next_patch;
       ip = (ip_next_block==0)*ip_next_patch;
       
       jm = (jm_next_block==0)*jm_next_patch;
       jp = (jp_next_block==0)*jp_next_patch;
       
       
       fprintf(fidin,'%d %d %d %d\n', [im ip jm jp]);
       
       if(im==0)
       fprintf(fidin,'%d %d\n', [im_next_block im_next_patch]);
       end
       
       if(ip==0)
       fprintf(fidin,'%d %d\n', [ip_next_block ip_next_patch]);
       end
       
       if(jm==0)
       fprintf(fidin,'%d %d\n', [jm_next_block jm_next_patch]);
       end
       
       if(jp==0)
       fprintf(fidin,'%d %d\n', [jp_next_block jp_next_patch]);
       end
       
       NI(ib) = ni;
       NJ(ib) = nj;
       
    end
    
    % write corners
    fprintf(fidin,'%d\n', ncorner);
    for n=1:ncorner
        cor_type(n) = 0;
    end
    
    for n=1:ncorner
    fprintf(fidin,'%d %d\n', [blk.corner{n}.Nb, cor_type(n)]);
    for nb=1:blk.corner{n}.Nb
    ib = blk.corner{n}.block{nb};
    ic = blk.corner{n}.i{nb};
    jc = blk.corner{n}.j{nb};
    if(ic>1); ic = NI(ib); end
    if(jc>1); jc = NJ(ib); end
    fprintf(fidin,'%d %d %d\n', [ib ic jc]);
    end
    
    end
    % 
    % write rest of file
    % check these are ok for your case!
        
        %nsteps nwrite ncut 
        fprintf(fidin,'%d %d %d\n', [solver.niter solver.nwrite solver.ncut]);
        
        % cfl, filter coefficient, ifsplit, ifsat, if_LES, if
        fprintf(fidin,'%f %f %d %d %d %d %d\n', [solver.cfl, solver.sigma solver.ifsplit solver.ifsat solver.ifLES 0 0]);
        
         % Toin poin pext vref alpha_in pitch_in aturb (not used) ilength (not used) g_z
        fprintf(fidin,'%f %f %f %f %f %f %f %f %d %d\n', [bcs.Toin bcs.Poin bcs.pexit bcs.vin bcs.alpha bcs.cax bcs.aturb bcs.lturb bcs.ilength bcs.radprof]);
        
        % gamma cp mu_ref sutherlands constants prandtl no.
        fprintf(fidin,'%f %f %12.5e %f %f %f\n', [gas.gamma gas.cp gas.mu_ref gas.mu_tref gas.mu_cref gas.pr]);
        
        % span, spanwise grid expansion factor
        fprintf(fidin,'%f %f\n', [solver.span solver.fexpan]);
        
        % restart, statistics
        fprintf(fidin,'%d %d\n', [solver.irestart solver.istats]);
        
        % inlet groups
        fprintf(fidin,'%d\n', [1]);
    
        if topology == 1
            fprintf(fidin,'%d\n', [2]);   
            fprintf(fidin,'%d\n%d\n', [1 2]);
        elseif topology == 2
            fprintf(fidin,'%d\n', [3]);   
            fprintf(fidin,'%d\n%d\n%d\n', [1 2 3]);
        elseif topology == 3 
            inlet_blocks = [1];
            while blk.next_block{inlet_blocks(end)}.jp ~= 0
                inlet_blocks(end+1) = blk.next_block{inlet_blocks(end)}.jp;
            end
            fprintf(fidin,'%d\n', length(inlet_blocks));
            for dum  = 1:length(inlet_blocks)
                fprintf(fidin,'%d\n', inlet_blocks(dum))
            end
        end
         
        % stability stuff   
        fprintf(fidin,'%d %d %d %d %d %f %f\n', [solver.istability 150 177 183 2 5.0 10000.0]);
       
        % inlet boundary layer stuff and adiabatic walls
        fprintf(fidin,'%d %f %f', [0 1e5 bcs.twall]);
    
    
    fclose(fidin);
    fprintf('Total processors: %d\n', nprocs)
end

end