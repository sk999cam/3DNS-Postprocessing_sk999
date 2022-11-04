function solver_process = mesh_and_run_Fluent(geom_params, folder, run_fluent)

    if nargin < 3
        run_fluent = true;
    end

    if ~exist(folder,'dir') && run_fluent
        mkdir(folder);
    end

    templatePath = 'E:\MATLAB\RANS_optimisation\optimisation_setup.jou';
    jouPath = fullfile(folder,'fluent_input.jou');
    
    if ispc
        MATLABpath = 'E:\MATLAB';
    end
    
    profile = fullfile(MATLABpath,'r150_aerofoil.txt')
    
    [xprof,yprof,pitch,stag] = read_profile(profile,true);
    amp = geom_params.amp;% 0.003;
    x1 = geom_params.x1;
    x2 = geom_params.x2;
    xTrip = geom_params.xTrip;
    r = geom_params.r;
    theta = geom_params.theta;
    tmid = geom_params.tmid;
    
    [xprof,yprof] = perturbed_profile(xprof,yprof,amp,x1,x2,xTrip,r,theta,tmid);
    
    msmooths=0;
    npp=6;
    pitch=2;
    bound_angle = 5.0;
    Lup=2.0;
    Ldn=1.0;
    Lo=0.05;
    ywall=Lo/npp;
    
    [blk,next_block,next_patch,corner] = single_blade_topology(xprof,yprof,pitch,bound_angle,Lup,Ldn,Lo,stag,npp,ywall,msmooths,0.1);
    
    NB = length(blk);

%     fig = figure;
%     hold on
%     for i=1:NB
%         xnew = blk{i}.x;
%         ynew = blk{i}.y;
%         plot(xnew,ynew,'k');
%         plot(xnew',ynew','k');
%     end
%     axis equal
    
    % Intermediate mesh imputs
    npp = 30; % points in i,j direction for each block are set to a multiple of npp- this ensures load balancing when running on many cores
    refine_fac = 5.0; % the grid dimensions are scaled by this number in both directions
    ywall = 0.0002; % the final near wall cell height
    msmooths = 30; % number of iterations for the poisson solver
    
    blk=mesh_refinement(blk,refine_fac,npp);
    
    blk=mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,xprof,yprof,ywall);
    
    
    blk = thin_blocks(blk, [1,2,3], [2,2,2], 2, npp);
    blk = thin_blocks(blk, [10,11,12], [1,1,1], 2, npp);
    blk = thin_blocks(blk, [3,8,12], [3,3,3], 2, npp);
    blk = thin_blocks(blk, [2,7,11], [4,4,4], 2, npp);
    
    
%     clear fig;
    fig = figure;
    hold on
    for i=1:NB
        xnew = blk{i}.x;
        ynew = blk{i}.y;
        plot(xnew,ynew,'k');
        plot(xnew',ynew','k');
    end
    axis equal
    
    for i=1:NB
        blk_dims(i,:) = size(blk{i}.x);
    end

    for i=1:NB
        blk2.x{i}=blk{i}.x;
        blk2.y{i}=blk{i}.y;
        blk2.nk{i}=1;
    end
    blk2.oblocks = [4 6 9 5];
    blk2.oblocks_flip = [0 0 1 1];
    blk2.span = 0.1;
    blk2.io_surfaces.blks = [1 2 2 7 11 11 10 12 12 8 3 3];
    blk2.io_surfaces.types = [1 1 3 3 3 2 2 2 4 4 4 1];
    blk2.blockdims = blk_dims;
    
    blk = blk2;
    
    counted_faces(1:NB, 1:4) = false;
    nnodes = 0;
    nfaces = 0;
    ncells = 0;
    
    for ib=1:NB
        [ni, nj] = size(blk.x{ib});
        nii = ni;
        nji = nj-1;
        nij = ni-1;
        njj = nj;
        ncells = ncells + (ni-1)*(nj-1);
        if counted_faces(ib, 1)
            ni = ni-1;
            nii = nii-1;
        end
        if counted_faces(ib, 2)
            ni = ni-1;
            nii = nii-1;
        end
        if counted_faces(ib, 3)
            nj = nj-1;
            njj = njj-1;
        end
        if counted_faces(ib, 4)
            nj = nj-1;
            njj = njj-1;
        end
        nnodes = nnodes + ni*nj;
        nfaces = nfaces + nii*nji + nij*njj;
        
        counted_faces = set_counted(ib, next_block, next_patch, counted_faces);
    end
    
    fprintf('Unstructured mesh stats:\n')
    fprintf('Nodes: %d, 0x%X\n', nnodes, nnodes);
    fprintf('Faces: %d, 0x%X\n', nfaces, nfaces);
    fprintf('Cells: %d, 0x%X\n', ncells, ncells);

    strs2switch = {}; new_strs = {};
    strs2switch{end+1} = 'RUNFOLDER'; new_strs{end+1} = folder;
    strs2switch{end+1} = 'INITDATA'; new_strs{end+1} = 'E:\MATLAB\RANS_optimisation\clean.dat';
    if exist(fullfile(folder,'loading.txt'),"file")
        strs2switch{end+1} = 'WRITELOADING'; new_strs{end+1} = '/plot/plot yes loading.txt yes yes no no absolute-pressure yes 1 0 upper-blade ()';
    else
        strs2switch{end+1} = 'WRITELOADING'; new_strs{end+1} = '/plot/plot yes loading.txt yes no no absolute-pressure yes 1 0 upper-blade ()';
    end


    if run_fluent
        writeFluentMesh(fullfile(folder, 'grid.msh'),blk,next_block,next_patch);
        switch_partial_strings(templatePath,jouPath,strs2switch,new_strs);
        fluentCommand = ['"C:\Program Files\ANSYS Inc\v212\fluent\ntbin\win64\fluent.exe" 2d -t32 -i "' jouPath '"']% ' -initial-file ' init_file]
        solver_process = processManager('command',fluentCommand);
    else
        solver_process = [];
    end
    
    
    function new_counts = set_counted(ib, next_block, next_patch, counts)
        new_counts = counts;
        for dr = ["im" "ip" "jm" "jp"]
            nblk = next_block{ib}.(dr);
            if nblk ~=0
                nptch = next_patch{ib}.(dr);
                new_counts(nblk,nptch) = true;
            end
        end
    end

end