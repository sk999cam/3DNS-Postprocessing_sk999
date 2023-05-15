function blkNodes = writeFluentMeshExtruded(path, blk, next_block, next_patch, nk, span, iWrite)

    if nargin < 7
        iWrite = true;
    end
    
    NB = length(blk.x);
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    
    nodes = [];
    faces = [];
    cells = [];
    nnodes = 1;
    ncells = 1;
    nfaces = 1;

    zv = linspace(0,span,nk);

    for ib = 1:NB
        [ni, nj] = size(blk.x{ib});
        blk.nodes{ib} = zeros(ni, nj, nk);
        blk.ifaces{ib} = zeros(ni, nj-1, nk-1);
        blk.jfaces{ib} = zeros(ni-1, nj, nk-1);
        blk.kfaces{ib} = zeros(ni-1, nj-1, nk);
    end
    
    for ib = 1:NB
    
        is = 1:size(blk.x{ib},1);
        js = 1:size(blk.x{ib},2);
    
        if im_assigned(ib)
            blk.nodes{ib}(1,:,:) = getBoarderNodes(ib, 'im');
            is = is(2:end);
        end
        if ip_assigned(ib)
            blk.nodes{ib}(end,:,:) = getBoarderNodes(ib, 'ip');
            is = is(1:end-1);
        end
        if jm_assigned(ib)
            blk.nodes{ib}(:,1,:) = getBoarderNodes(ib, 'jm');
            js = js(2:end);
        end
        if jp_assigned(ib)
            blk.nodes{ib}(:,end,:) = getBoarderNodes(ib, 'jp');
            js = js(1:end-1);
        end
    
        for k=1:nk
            for i=is
                for j=js
                    blk.nodes{ib}(i,j,k) = nnodes;
                    nodes(nnodes).x = blk.x{ib}(i,j);
                    nodes(nnodes).y = blk.y{ib}(i,j);
                    nodes(nnodes).z = zv(k);
                    nnodes = nnodes+1;
                end
            end
        end

        setAssigned(ib, 'im');
        setAssigned(ib, 'ip');
        setAssigned(ib, 'jm');
        setAssigned(ib, 'jp');
    
    end

    % Set cell nos

    for ib = 1:NB
        for k=1:nk-1
            for i = 1:size(blk.x{ib},1)-1
                for j = 1:size(blk.x{ib},2)-1
                    blk.cells{ib}(i,j,k) = ncells;
                    ncells = ncells+1;
                end
            end
        end
    end

    % Set face nos
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    boundary_zone(1:NB,1:4) = false;

    % Inlet
    inlet_faces(1)=nfaces;
    for ib = [1 2 3]
        for k=1:nk-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(1,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(1,j+1);
                blk.ifaces{ib}(1,j,k) = nfaces;
                faces(nfaces).bc = 4;
                if nfaces == 60
                    disp('')
                end
                nfaces = nfaces+1;
                
            end
        end
        boundary_zone(ib,1) = true;
    end
    inlet_faces(2)=nfaces-1;

    % Upper freestream
    uPer_faces(1)=nfaces;
    for ib = [3 8 12]
        for k=1:nk-1
            for i = 1:size(blk.x{ib},1)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(i,end);
    %             faces(nfaces).n2 = blk.nodes{ib}(i+1,end);
                blk.jfaces{ib}(i,end,k) = nfaces;
                faces(nfaces).bc = 5;
                nfaces = nfaces+1;
            end
        end
        boundary_zone(ib,4) = true;
    end
    uPer_faces(2)=nfaces-1;

    % Lower freestream
    lPer_faces(1)=nfaces;
    for ib = [2 7 11]
        for k=1:nk-1
            for i = 1:size(blk.x{ib},1)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(i,1);
    %             faces(nfaces).n2 = blk.nodes{ib}(i+1,1);
                blk.jfaces{ib}(i,1,k) = nfaces;
                faces(nfaces).bc = 5;
                nfaces = nfaces+1;
            end
        end
        boundary_zone(ib,3) = true;
    end
    lPer_faces(2)=nfaces-1;

    % Outlet
    outlet_faces(1)=nfaces;
    for ib = [10 11 12]
        for k=1:nk-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(end,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
                blk.ifaces{ib}(end,j,k) = nfaces;
                faces(nfaces).bc = 5;
                nfaces = nfaces+1;
            end
        end
        boundary_zone(ib,2) = true;
    end
    outlet_faces(2)=nfaces-1;

    % Blade
    upper_blade_faces(1)=nfaces;
    nLowerFaces = 0;
    lower_faces = [];
    for ib = [4 5 6 9]
        for k=1:nk-1
            for i = 1:size(blk.x{ib},1)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(i,end);
    %             faces(nfaces).n2 = blk.nodes{ib}(i+1,end);
                midpoint = (blk.y{ib}(i,end)+blk.y{ib}(i+1,end))/2;
                if midpoint > 0
                    blk.jfaces{ib}(i,end,k) = nfaces;
                    faces(nfaces).bc = 3;
                    nfaces = nfaces+1;
                else
                    nLowerFaces = nLowerFaces + 1;
                    lower_faces(nLowerFaces,:) = [ib i k];
                end
    
            end
        end
        boundary_zone(ib,4) = true;
    end
    upper_blade_faces(2)=nfaces-1;
    
    lower_blade_faces(1)=nfaces;
    for iface=1:nLowerFaces
        blk.jfaces{lower_faces(iface,1)}(lower_faces(iface,2),end,lower_faces(iface,3)) = nfaces;
        faces(nfaces).bc = 3;
        nfaces = nfaces+1;
    end
    lower_blade_faces(2)=nfaces-1;

    

    % zm peroiodic
    zm_faces(1) = nfaces;
    for ib = 1:12
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(end,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
                blk.kfaces{ib}(i,j,1) = nfaces;
                faces(nfaces).bc = 12;
                nfaces = nfaces+1;
            end
        end
    end
    zm_faces(2)=nfaces-1;

    % zp peroiodic
    zp_faces(1) = nfaces;
    for ib = 1:12
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(end,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
                blk.kfaces{ib}(i,j,end) = nfaces;
                faces(nfaces).bc = 8;
                nfaces = nfaces+1;
            end
        end
    end
    zp_faces(2)=nfaces-1;


    for ib = 1:NB
        fprintf('Calculating block %d/%d\n',ib,NB)
        [ni, nj] = size(blk.x{ib});
        if ib == 1
            ni;
        end

        % Set remaining boarder faces
        if ~boundary_zone(ib,1)
            if ~im_assigned(ib)
                for k=1:nk-1
                    for j = 1:nj-1
                        blk.ifaces{ib}(1,j,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'im')
            else
                blk.ifaces{ib}(1,:,:) = getBoarderFaces(ib,'im');
            end
        end

        if ~boundary_zone(ib,2)
            if ~ip_assigned(ib)
                for k=1:nk-1
                    for j = 1:nj-1
                        blk.ifaces{ib}(end,j,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'ip');
            else
                blk.ifaces{ib}(end,:,:) = getBoarderFaces(ib,'ip');
            end
        end

        if ~boundary_zone(ib,3)
            if ~jm_assigned(ib)
                for k=1:nk-1
                    for i = 1:ni-1
                        blk.jfaces{ib}(i,1,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'jm');
            else
                blk.jfaces{ib}(:,1,:) = getBoarderFaces(ib,'jm');
            end
        end

        if ~boundary_zone(ib,4)
            if ~jp_assigned(ib)
                for k=1:nk-1
                    for i = 1:ni-1
                        blk.jfaces{ib}(i,end,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'jp');
            else
                blk.jfaces{ib}(:,end,:) = getBoarderFaces(ib,'jp');
            end
        end

        % Set interior face nos

        % i faces
        for k=1:nk-1
            for i=2:ni-1
                for j=1:nj-1
                    blk.ifaces{ib}(i,j,k) = nfaces;
                    if nfaces == 605942
                        pause
                    end
                    nfaces = nfaces+1;
                end
            end
        end
        % j faces
        for k=1:nk-1
            for j=2:nj-1
                for i=1:ni-1
                    blk.jfaces{ib}(i,j,k) = nfaces;
                    if nfaces == 605942
                        pause
                    end
                    nfaces = nfaces+1;
                end
            end
        end
        % k faces
        for k=2:nk-1
            for j=1:nj-1
                for i=1:ni-1
                    blk.kfaces{ib}(i,j,k) = nfaces;
                    if nfaces == 605942
                        pause
                    end
                    nfaces = nfaces+1;
                end
            end
        end

        % Now get face nodes and cells
        
        % i faces
        for i=1:ni
            mcells = zeros(1,nj-1,nk-1);
            pcells = zeros(1,nj-1,nk-1);
            if i==1
                mcells = getBoarderCells(ib, 'im');
                pcells = blk.cells{ib}(i,:,:);
            elseif i==ni
                mcells = blk.cells{ib}(i-1,:,:);
                pcells = getBoarderCells(ib, 'ip');
            else
                mcells = blk.cells{ib}(i-1,:,:);
                pcells = blk.cells{ib}(i,:,:);
            end

            if ismember(ib, [4 6])
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for k=1:nk-1
                for j=1:nj-1
                    n = blk.ifaces{ib}(i,j,k);
                    if (i==1) && (j==1) && (k==2)
                        disp('')
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i,j+1,k);
                    faces(n).n3 = blk.nodes{ib}(i,j+1,k+1);
                    faces(n).n4 = blk.nodes{ib}(i,j,k+1);
                    faces(n).cr = pcells(1,j,k);
                    faces(n).cl = mcells(1,j,k);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                end
            end
        end

        % j faces
        for j=1:nj
            mcells = zeros(ni-1,1,nk-1);
            pcells = zeros(ni-1,1,nk-1);
            if j==1
                mcells = blk.cells{ib}(:,j,:);
                pcells = getBoarderCells(ib, 'jm');
            elseif j==nj
                mcells = getBoarderCells(ib, 'jp');
                pcells = blk.cells{ib}(:,j-1,:);
            else
                mcells = blk.cells{ib}(:,j,:);
                pcells = blk.cells{ib}(:,j-1,:);
            end

            if ismember(ib, [4 6])
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for k=1:nk-1
                for i=1:ni-1
                    n = blk.jfaces{ib}(i,j,k);
                    if n == 0
                        n
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i+1,j,k);
                    faces(n).n3 = blk.nodes{ib}(i+1,j,k+1);
                    faces(n).n4 = blk.nodes{ib}(i,j,k+1);
                    faces(n).cr = pcells(i,1,k);
                    faces(n).cl = mcells(i,1,k);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                end
            end
        end

        % k faces
        for k=1:nk
            mcells = zeros(ni-1,nj-1,1);
            pcells = zeros(ni-1,nj-1,1);
            if k==1
                pcells = blk.cells{ib}(:,:,k);
                %pcells = getBoarderCells(ib, 'jm');
            elseif k==nk
                %mcells = getBoarderCells(ib, 'jp');
                mcells = blk.cells{ib}(:,:,k-1);
            else
                pcells = blk.cells{ib}(:,:,k);
                mcells = blk.cells{ib}(:,:,k-1);
            end

            if ismember(ib, [4 6])
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for i=1:ni-1
                for j=1:nj-1
                    n = blk.kfaces{ib}(i,j,k);
                    if n == 0
                        n
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i+1,j,k);
                    faces(n).n3 = blk.nodes{ib}(i+1,j+1,k);
                    faces(n).n4 = blk.nodes{ib}(i,j+1,k);
                    faces(n).cr = pcells(i,j,1);
                    faces(n).cl = mcells(i,j,1);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                    if k==nk
                        faces(n).sf = blk.kfaces{ib}(i,j,1);
                    end
                end
            end
        end
    end

    fprintf('nnodes: %d, should be %d\n', length(nodes), ni*nj*nk)
    fprintf('ncells: %d, should be %d\n', ncells-1, ((ni-1)*(nj-1)*(nk-1)))
    fprintf('nfaces: %d, should be %d\n', length(faces), ((ni-1)*nj*(nk-1) + ni*(nj-1)*(nk-1) + (ni-1)*(nj-1)*nk))

    facesStr = cell(length(faces),1);
    nodesStr = cell(length(nodes),1);
    shadowStr = cell(1+zp_faces(2)-zp_faces(1),1);

    fprintf('Creating ASCII face information to write\n')
    parfor n=1:length(faces)
        facesStr{n} = sprintf('%X %X %X %X %X %X',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
                faces(n).cr faces(n).cl]);
    end

    fprintf('Creating ASCII node information to write\n')
    parfor n=1:length(nodes)
        nodesStr{n} = sprintf('%10.8f %10.8f %10.8f', [nodes(n).x nodes(n).y nodes(n).z]);
    end

    fprintf('Creating ASCII shadow zone information\n')
    parfor i=1:length(shadowStr)
        n = i-1+zp_faces(1);
        shadowStr{i} = sprintf('%X %X',[faces(n).sf n]);
    end

    blkNodes = blk.nodes;

    if iWrite
        
        path
        fid = fopen(path,'w')
    
        fprintf(fid,'(0 "Grid:")\n');
        fprintf(fid,'\n');
    
        fprintf(fid,'(0 "Dimensions:")\n');
        fprintf(fid,'(2 3)\n');
        fprintf(fid,'\n');
        
        fprintf(fid,'(0 "Cells:")\n');
        fprintf(fid,'(12 (0 %X %X 0))\n',[1 (ncells-1)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Faces:")\n');
        fprintf(fid,'(13 (0 %X %X 0))\n',[1 length(faces)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Nodes:")\n');
        fprintf(fid,'(10 (0 %X %X 0 3))\n',[1 length(nodes)]);
        fprintf(fid, '\n');
        
        fprintf(fid,'(12 (b 1 %X 1 4))\n',ncells-1);
        fprintf(fid, '\n');

        fprintf('Writing interior faces\n')
        fprintf(fid,'(0 "Interior:")\n');
        fprintf(fid,'(13 (2 %X %X %X 4)(\n', [zp_faces(2)+1 length(faces) 2]);
%         for n=(zp_faces(2)+1):length(faces)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{zp_faces(2)+1:end});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf('Writing inlet faces\n')
        fprintf(fid,'(0 "Inlet:")\n');
        fprintf(fid,'(13 (3 %X %X %X 4)(\n', [inlet_faces(1) inlet_faces(2) 4]);
%         for n=inlet_faces(1):inlet_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{inlet_faces(1):inlet_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
        
        fprintf('Writing upper freestream faces\n')
        fprintf(fid,'(0 "Upper freestream:")\n');
        fprintf(fid,'(13 (4 %X %X %X 4)(\n', [uPer_faces(1) uPer_faces(2) 5]);
%         for n=uPer_faces(1):uPer_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{uPer_faces(1):uPer_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf('Writing lower freestream faces\n')
        fprintf(fid,'(0 "Lower freestream:")\n');
        fprintf(fid,'(13 (5 %X %X %X 4)(\n', [lPer_faces(1) lPer_faces(2) 5]);
%         for n=lPer_faces(1):lPer_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{lPer_faces(1):lPer_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf('Writing outlet faces\n')
        fprintf(fid,'(0 "Outlet:")\n');
        fprintf(fid,'(13 (6 %X %X %X 4)(\n', [outlet_faces(1) outlet_faces(2) 5]);
%         for n=outlet_faces(1):outlet_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{outlet_faces(1):outlet_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf('Writing blade faces\n')
        fprintf(fid,'(0 "Upper Blade:")\n');
        fprintf(fid,'(13 (7 %X %X %X 4)(\n', [upper_blade_faces(1) upper_blade_faces(2) 3]);
%         for n=upper_blade_faces(1):upper_blade_faces(2)
%             %fprintf('%d, n1: %d %X, n2: %d %X, cr: %d %X, cl: %d %X\n', [n faces(n).n1 faces(n).n1 faces(n).n2 faces(n).n2 faces(n).cr faces(n).cr faces(n).cl faces(n).cl])
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{upper_blade_faces(1):upper_blade_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Lower Blade:")\n');
        fprintf(fid,'(13 (8 %X %X %X 4)(\n', [lower_blade_faces(1) lower_blade_faces(2) 3]);
%         for n=lower_blade_faces(1):lower_blade_faces(2)
%             %fprintf('%d, n1: %d %X, n2: %d %X, cr: %d %X, cl: %d %X\n', [n faces(n).n1 faces(n).n1 faces(n).n2 faces(n).n2 faces(n).cr faces(n).cr faces(n).cl faces(n).cl])
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{lower_blade_faces(1):lower_blade_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');

        fprintf('Writing periodic faces\n')
        fprintf(fid,'(0 "z0 periodic:")\n');
        fprintf(fid,'(13 (9 %X %X %X 4)(\n', [zm_faces(1) zm_faces(2) 12]);
%         for n=zm_faces(1):zm_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{zm_faces(1):zm_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');

        fprintf(fid,'(0 "zp shadow periodic:")\n');
        fprintf(fid,'(13 (a %X %X %X 4)(\n', [zp_faces(1) zp_faces(2) 8]);
%         for n=zp_faces(1):zp_faces(2)
%             fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
%                 faces(n).cr faces(n).cl]);
%         end
        fprintf(fid, '%s\n', facesStr{zp_faces(1):zp_faces(2)});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');

        fprintf('Writing shadow zone\n')
        fprintf(fid,'(0 "Shadow zone:")\n');
        fprintf(fid,'(18 (%X %X %X %X)(\n', [1 (1+zp_faces(2)-zp_faces(1)) 9 10]);
%         for n=zp_faces(1):zp_faces(2)
%             fprintf(fid,'%X %X\n',[faces(n).sf n]);
%         end
        fprintf(fid, '%s\n', shadowStr{:});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
    
        fprintf('Writing nodes\n')
        fprintf(fid,'(10 (1 %X %X 1 3)\n(\n',[1 length(nodes)]);
%         for i=1:length(nodes)
%             fprintf(fid,'%10.8f %10.8f %10.8f\n', [nodes(i).x nodes(i).y nodes(i).z]);
%         end
        fprintf(fid, '%s\n', nodesStr{:});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
        
        fclose(fid);
    end

    function nodesnow = getBoarderNodes(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
        switch nptch
            case 1
            nodesnow = blk.nodes{nblk}(1,:,:);
            if ismember(dr,{'jm','jp'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 2
                nodesnow = blk.nodes{nblk}(end,:,:);
            if ismember(dr,{'jm','jp'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 3
                nodesnow = blk.nodes{nblk}(:,1,:);
            if ismember(dr,{'im','ip'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 4
                nodesnow = blk.nodes{nblk}(:,end,:);
            if ismember(dr,{'im','ip'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
        end
    end

    function facesnow = getBoarderFaces(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    facesnow = blk.ifaces{nblk}(1,:,:);
                    if ismember(dr,{'jm','jp'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 2
                    facesnow = blk.ifaces{nblk}(end,:,:);
                    if ismember(dr,{'jm','jp'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 3
                    facesnow = blk.jfaces{nblk}(:,1,:);
                    if ismember(dr,{'im','ip'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 4
                    facesnow = blk.jfaces{nblk}(:,end,:);
                    if ismember(dr,{'im','ip'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
            end
        end

    function cellsnow = getBoarderCells(ib, dr)
        nblk = next_block{ib}.(dr);
        if nblk > NB
            nblk = 0;
        end
        if nblk == 0
            [nib, njb] = size(blk.x{ib});
            switch dr
                case {'im','ip'}
                    cellsnow(1,1:njb-1,1:nk-1) = 0;
                case {'jm','jp'}
                    cellsnow(1:nib-1,1,1:nk-1) = 0;
            end
        else
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    cellsnow = blk.cells{nblk}(1,:,:);
                    if ismember(dr,{'jm','jp'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 2
                    cellsnow = blk.cells{nblk}(end,:,:);
                    if ismember(dr,{'jm','jp'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 3
                    cellsnow = blk.cells{nblk}(:,1,:);
                    if ismember(dr,{'im','ip'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 4
                    cellsnow = blk.cells{nblk}(:,end,:);
                    if ismember(dr,{'im','ip'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
            end
        end
    end

    function setAssigned(ib, dr)
        nblk = next_block{ib}.(dr);
        if nblk > NB
            nblk = 0;
        end
        if nblk ~= 0
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    im_assigned(nblk) = true;
                case 2
                    ip_assigned(nblk) = true;
                case 3
                    jm_assigned(nblk) = true;
                case 4
                    jp_assigned(nblk) = true;
            end
        end
    end

end
