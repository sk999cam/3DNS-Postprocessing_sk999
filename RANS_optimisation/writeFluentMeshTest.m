function writeFluentMeshTest(path, blk, next_block, next_patch)

    NB = 1 % length(blk.x);
    
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

    for ib = 1:NB
        [ni, nj] = size(blk.x{ib});
        blk.nodes{ib} = zeros(ni, nj);
        blk.ifaces{ib} = zeros(ni, nj-1);
        blk.jfaces{ib} = zeros(ni-1, nj);
    end
    
    for ib = 1:NB
    
        is = 1:size(blk.x{ib},1);
        js = 1:size(blk.x{ib},2);
    
        if im_assigned(ib)
            blk.nodes{ib}(1,:) = getBoarderNodes(next_block{ib}.im, next_patch{ib}.im);
            is = is(2:end);
        end
        if ip_assigned(ib)
            blk.nodes{ib}(end,:) = getBoarderNodes(next_block{ib}.ip, next_patch{ib}.ip);
            is = is(1:end-1);
        end
        if jm_assigned(ib)
            blk.nodes{ib}(:,1) = getBoarderNodes(next_block{ib}.jm, next_patch{ib}.jm);
            js = js(2:end);
        end
        if jp_assigned(ib)
            blk.nodes{ib}(:,end) = getBoarderNodes(next_block{ib}.jp, next_patch{ib}.jp);
            js = js(1:end-1);
        end
    
        for i=is
            for j=js
                blk.nodes{ib}(i,j) = nnodes;
                nodes(nnodes).x = blk.x{ib}(i,j);
                nodes(nnodes).y = blk.y{ib}(i,j);
                nnodes = nnodes+1;
            end
        end

        setAssigned(ib, 'im');
        setAssigned(ib, 'ip');
        setAssigned(ib, 'jm');
        setAssigned(ib, 'jp');
    
    end

    % Set cell nos

    for ib = 1:NB
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
                blk.cells{ib}(i,j) = ncells;
                ncells = ncells+1;
            end
        end
    end

    % Set face nos
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    face_zones(1:NB,1:4) = false;

    % Inlet
    inlet_faces(1)=nfaces;
    for ib = [1]
        for j = 1:size(blk.x{ib},2)-1
%             faces(nfaces).n1 = blk.nodes{ib}(1,j);
%             faces(nfaces).n2 = blk.nodes{ib}(1,j+1);
            blk.ifaces{ib}(1,j) = nfaces;
            faces(nfaces).bc = 4;
            nfaces = nfaces+1;
            
        end
        face_zone(ib,1) = true;
    end
    inlet_faces(2)=nfaces-1;

    % Upper periodic
    uPer_faces(1)=nfaces;
    for ib = [1]
        for i = 1:size(blk.x{ib},1)-1
%             faces(nfaces).n1 = blk.nodes{ib}(i,end);
%             faces(nfaces).n2 = blk.nodes{ib}(i+1,end);
            blk.jfaces{ib}(i,end) = nfaces;
            faces(nfaces).bc = 5;
            nfaces = nfaces+1;
        end
        face_zone(ib,4) = true;
    end
    uPer_faces(2)=nfaces-1;

    % Lower periodic
    lPer_faces(1)=nfaces;
    for ib = [1]
        for i = 1:size(blk.x{ib},1)-1
%             faces(nfaces).n1 = blk.nodes{ib}(i,1);
%             faces(nfaces).n2 = blk.nodes{ib}(i+1,1);
            blk.jfaces{ib}(i,1) = nfaces;
            faces(nfaces).bc = 5;
            nfaces = nfaces+1;
        end
        face_zone(ib,3) = true;
    end
    lPer_faces(2)=nfaces-1;

    % Outlet
    outlet_faces(1)=nfaces;
    for ib = [1]
        for j = 1:size(blk.x{ib},2)-1
%             faces(nfaces).n1 = blk.nodes{ib}(end,j);
%             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
            blk.ifaces{ib}(end,j) = nfaces;
            faces(nfaces).bc = 5;
            nfaces = nfaces+1;
        end
        face_zone(ib,2) = true;
    end
    outlet_faces(2)=nfaces-1;

%     % Blade
%     blade_faces(1)=nfaces;
%     for ib = [4 5 6 9]
%         for i = 1:size(blk.x{ib},1)-1
% %             faces(nfaces).n1 = blk.nodes{ib}(i,end);
% %             faces(nfaces).n2 = blk.nodes{ib}(i+1,end);
%             blk.jfaces{ib}(i,end) = nfaces;
%             faces(nfaces).bc = 3;
%             nfaces = nfaces+1;
%         end
%         face_zone(ib,4) = true;
%     end
%     blade_faces(2)=nfaces-1;

    for ib = 1:NB
        fprintf('Calculating block %d/%d\n',ib,NB)
        [ni, nj] = size(blk.x{ib});
        if ib == 1
            ni;
        end

        % Set remaining boarder faces
        if ~face_zone(ib,1)
            if ~im_assigned(ib)
                for j = 1:nj-1
                    blk.ifaces{ib}(1,j) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'im')
            else
                blk.ifaces{ib}(1,:) = getBoarderFaces(ib,'im');
            end
        end

        if ~face_zone(ib,2)
            if ~ip_assigned(ib)
                for j = 1:nj-1
                    blk.ifaces{ib}(end,j) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'ip');
            else
                blk.ifaces{ib}(end,:) = getBoarderFaces(ib,'ip');
            end
        end

        if ~face_zone(ib,3)
            if ~jm_assigned(ib)
                for i = 1:ni-1
                    blk.jfaces{ib}(i,1) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'jm');
            else
                blk.jfaces{ib}(:,1) = getBoarderFaces(ib,'jm');
            end
        end

        if ~face_zone(ib,4)
            if ~jp_assigned(ib)
                for i = 1:ni-1
                    blk.jfaces{ib}(i,end) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'jp');
            else
                blk.jfaces{ib}(:,end) = getBoarderFaces(ib,'jp');
            end
        end

        % Set interior face nos

        % i faces
        for i=2:ni-1
            for j=1:nj-1
                blk.ifaces{ib}(i,j) = nfaces;
                nfaces = nfaces+1;
            end
        end
        % j faces
        for j=2:nj-1
            for i=1:ni-1
                blk.jfaces{ib}(i,j) = nfaces;
                nfaces = nfaces+1;
            end
        end

        % Now get face nodes and cells
        
        % i faces
        for i=1:ni
            mcells = zeros(1,nj-1);
            pcells = zeros(1,nj-1);
            if i==1
                mcells = getBoarderCells(ib, 'im');
                pcells = blk.cells{ib}(i,:);
            elseif i==ni
                mcells = blk.cells{ib}(i-1,:);
                pcells = getBoarderCells(ib, 'ip');
            else
                mcells = blk.cells{ib}(i-1,:);
                pcells = blk.cells{ib}(i,:);
            end

            for j=1:nj-1
                n = blk.ifaces{ib}(i,j);
                faces(n).n1 = blk.nodes{ib}(i,j);
                faces(n).n2 = blk.nodes{ib}(i,j+1);
                faces(n).cr = mcells(j);
                faces(n).cl = pcells(j);
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end

        % j faces
        for j=1:nj
            mcells = zeros(1,nj-1);
            pcells = zeros(1,nj-1);
            if j==1
                mcells = blk.cells{ib}(:,j);
                pcells = getBoarderCells(ib, 'jm');
            elseif j==nj
                mcells = getBoarderCells(ib, 'jp');
                pcells = blk.cells{ib}(:,j-1);
            else
                mcells = blk.cells{ib}(:,j);
                pcells = blk.cells{ib}(:,j-1);
            end

            for i=1:ni-1
                n = blk.jfaces{ib}(i,j);
                if n == 0
                    n
                end
                faces(n).n1 = blk.nodes{ib}(i,j);
                faces(n).n2 = blk.nodes{ib}(i+1,j);
                faces(n).cr = mcells(i);
                faces(n).cl = pcells(i);
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end
    end

    fprintf('nnodes: %d, should be %d\n', length(nodes), ni*nj)
    fprintf('ncells: %d, should be %d\n', ncells-1, ((ni-1)*(nj-1)))
    fprintf('nfaces: %d, should be %d\n', length(faces), ((ni-1)*nj + ni*(nj-1)))

    fid = fopen(path,'w');

    fprintf(fid,'(0 "Grid:")\n');
    fprintf(fid,'\n');

    fprintf(fid,'(0 "Dimensions:")\n');
    fprintf(fid,'(2 2)\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'(0 "Cells:")\n');
    fprintf(fid,'(12 (0 %X %X 0))\n',[1 (ncells-1)]);
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Faces:")\n');
    fprintf(fid,'(13 (0 %X %X 0))\n',[1 length(faces)]);
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Nodes:")\n');
    fprintf(fid,'(10 (0 %X %X 0 2))\n',[1 length(nodes)]);
    fprintf(fid, '\n');
    
    fprintf(fid,'(12 (8 1 %X 1 3))\n',ncells-1);
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Interior:")\n');
    fprintf(fid,'(13 (2 %X %X %X 2)(\n', [outlet_faces(2)+1 length(faces) 2]);
    for n=(outlet_faces(2)+1):length(faces)
        fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Inlet:")\n');
    fprintf(fid,'(13 (2 %X %X %X 2)(\n', [inlet_faces(1) inlet_faces(2) 4]);
    for n=inlet_faces(1):inlet_faces(2)
        fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Upper periodic:")\n');
    fprintf(fid,'(13 (2 %X %X %X 2)(\n', [uPer_faces(1) uPer_faces(2) 5]);
    for n=uPer_faces(1):uPer_faces(2)
        fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Lower periodic:")\n');
    fprintf(fid,'(13 (3 %X %X %X 2)(\n', [lPer_faces(1) lPer_faces(2) 5]);
    for n=lPer_faces(1):lPer_faces(2)
        fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');

    fprintf(fid,'(0 "Outlet:")\n');
    fprintf(fid,'(13 (4 %X %X %X 2)(\n', [outlet_faces(1) outlet_faces(2) 5]);
    for n=outlet_faces(1):outlet_faces(2)
        fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');

%     fprintf(fid,'(0 "Blade:")\n');
%     fprintf(fid,'(13 (5 %X %X %X 2)(\n', [blade_faces(1) blade_faces(2) 3]);
%     for n=blade_faces(1):blade_faces(2)
%         fprintf('%d, n1: %d %X, n2: %d %X, cr: %d %X, cl: %d %X\n', [n faces(n).n1 faces(n).n1 faces(n).n2 faces(n).n2 faces(n).cr faces(n).cr faces(n).cl faces(n).cl])
%         fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
%     end
%     fprintf(fid, '))\n');
%     fprintf(fid, '\n');

    fprintf(fid,'(10 (1 %X %X 1 2)\n(\n',[1 length(nodes)]);
    for i=1:length(nodes)
        fprintf(fid,'%10.8f %10.8f\n', [nodes(i).x nodes(i).y]);
    end
    fprintf(fid, '))\n');
    fprintf(fid, '\n');
    
    fclose(fid);

    function nodesnow = getBoarderNodes(nblk, nptch)
        switch nptch
            case 1
                nodesnow = blk.nodes{nblk}(1,:);
            case 2
                nodesnow = blk.nodes{nblk}(end,:);
            case 3
                nodesnow = blk.nodes{nblk}(:,1);
            case 4
                nodesnow = blk.nodes{nblk}(:,end);
        end
    end

    function facesnow = getBoarderFaces(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    facesnow = blk.ifaces{nblk}(1,:);
                case 2
                    facesnow = blk.ifaces{nblk}(end,:);
                case 3
                    facesnow = blk.jfaces{nblk}(:,1);
                case 4
                    facesnow = blk.jfaces{nblk}(:,end);
            end
        end

    function cellsnow = getBoarderCells(ib, dr)
        nblk = 0;%next_block{ib}.(dr);
        if nblk == 0
            [nib, njb] = size(blk.x{ib});
            switch dr
                case {'im','ip'}
                    cellsnow(1:njb-1) = 0;
                case {'jm','jp'}
                    cellsnow(1:nib-1) = 0;
            end
        else
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    cellsnow = blk.cells{nblk}(1,:);
                case 2
                    cellsnow = blk.cells{nblk}(end,:);
                case 3
                    cellsnow = blk.cells{nblk}(:,1);
                case 4
                    cellsnow = blk.cells{nblk}(:,end);
            end
        end
    end

    function setAssigned(ib, dr)
        nblk = 0;%next_block{ib}.(dr);
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
