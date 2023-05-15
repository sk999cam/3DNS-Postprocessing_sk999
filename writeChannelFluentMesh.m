function blkNodes = writeChannelFluentMesh(blk, path, nk, span, iWrite)

    if nargin < 5
        iWrite = true;
    end

    if nargin < 3 || isempty(nk)
        span = 0;
        nk = 1;
    end

    nkm1 = max(1,nk-1);

    ib = 1;

    nodelist = [];
    faces = [];
    celllist = [];
    nnodes = 1;
    nfaces = 1;

    [ni, nj] = size(blk.x{1});
    ncells = (ni-1)*(nj-1)*nkm1;
    nodes = zeros(ni, nj, nk);
    ifaces = zeros(ni, nj-1, nkm1);
    jfaces = zeros(ni-1, nj, nkm1);
    kfaces = zeros(ni-1, nj-1, nk);
    cells = [];

    is = 1:ni;
    js = 1:nj;
    ks = 1:nk;

    zv = linspace(0,span,nk);

    % Nodes
    for i=is
        for j=js
            for k = ks
                nodes(i,j,k) = nnodes;
                nodelist(nnodes).x = blk.x{ib}(i,j);
                nodelist(nnodes).y = blk.y{ib}(i,j);
                nodelist(nnodes).z = zv(k);
                nnodes = nnodes+1;
            end
        end
    end

    % Cells

    cells = reshape(1:ncells,ni-1,nj-1,nkm1);

    % Inlet

    inlet_faces(1)=nfaces;
    for k = 1:nkm1
        for j = 1:nj-1
            ifaces(1,j,k) = nfaces;
            faces(nfaces).bc = 4;
            nfaces = nfaces+1;
        end
    end
    inlet_faces(2) = nfaces-1;

    % Inviscid wall

    invc_faces(1) = nfaces;
    for k = 1:nkm1
        for i = 1:ni-1
            jfaces(i,end,k) = nfaces;
            faces(nfaces).bc = 3;
            nfaces = nfaces+1;
        end
    end
    invc_faces(2) = nfaces-1;

    % Outlet

    outlet_faces(1) = nfaces;
    for k = 1:nkm1
        for j = 1:nj-1
            ifaces(end,j,k) = nfaces;
            faces(nfaces).bc = 5;
            nfaces = nfaces+1;
        end
    end
    outlet_faces(2) = nfaces-1;

    % Viscous wall (blade)

    wall_faces(1) = nfaces;
    for k = 1:nkm1
        for i = 1:ni-1
            jfaces(i,1,k) = nfaces;
            faces(nfaces).bc = 3;
            nfaces = nfaces+1;
        end
    end
    wall_faces(2) = nfaces-1;

    if nk > 1

        % Inviscid walls (z=0, z=span)
        
        zm_faces(1) = nfaces;
        for k=[1 nk]
            for j = 1:nj-1
                for i = 1:ni-1
                    kfaces(i,j,k) = nfaces;
                    faces(nfaces).bc = 3;
                    nfaces = nfaces+1;
                end
            end
            if k == 1
                zm_faces(2) = nfaces-1;
                zp_faces(1) = nfaces;
            end
        end
        zp_faces(2) = nfaces-1;
    end

    % Interior faces

    interior_faces(1) = nfaces;

    % i faces
    for k=1:nkm1
        for j=1:nj-1
            for i=2:ni-1
                ifaces(i,j,k) = nfaces;
                nfaces = nfaces+1;
            end
        end
    end
    % j faces
    for k=1:nkm1
        for j=2:nj-1
            for i=1:ni-1
                jfaces(i,j,k) = nfaces;
                nfaces = nfaces+1;
            end
        end
    end
    % k faces
    if nk > 2
        for k=2:nkm1
            for j=1:nj-1
                for i=1:ni-1
                    kfaces(i,j,k) = nfaces;
                    nfaces = nfaces+1;
                end
            end
        end
    end

    interior_faces(2) = nfaces-1;

    % Now get face nodes and cells
        
    % i faces

    for i=1:ni
        mcells = zeros(1,nj-1,nkm1);
        pcells = zeros(1,nj-1,nkm1);
        if i==1
            pcells = cells(i,:,:);
        elseif i == ni
            mcells = cells(end,:,:);
        else
            mcells = cells(i-1,:,:);
            pcells = cells(i,:,:);
        end

        for k=1:nkm1
            for j=1:nj-1
                n = ifaces(i,j,k);
                faces(n).n1 = nodes(i,j,k);
                faces(n).n2 = nodes(i,j+1,k);
                faces(n).cr = mcells(1,j,k);
                faces(n).cl = pcells(1,j,k);
                if nk > 1
                    faces(n).n3 = nodes(i,j+1,k+1);
                    faces(n).n4 = nodes(i,j,k+1);
                end
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end
    end

    % j faces

    for j=1:nj
        mcells = zeros(ni-1, 1, nkm1);
        pcells = zeros(ni-1, 1, nkm1);

        if j==1
            mcells = cells(:,j,:);
        elseif j == nj
            pcells = cells(:,end,:);
        else
            mcells = cells(:, j, :);
            pcells = cells(:, j-1, :);
        end
        for k=1:nkm1
            for i=1:ni-1
                n = jfaces(i,j,k);
                faces(n).n1 = nodes(i,j,k);
                faces(n).n2 = nodes(i+1,j,k);
                faces(n).cr = mcells(i,1,k);
                faces(n).cl = pcells(i,1,k);
                if nk > 1
                    faces(n).n3 = nodes(i+1,j,k+1);
                    faces(n).n4 = nodes(i,j,k+1);
                end
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end
    end

    % k faces
    if nk > 1
    
        for k=1:nk
            mcells = zeros(ni-1, nj-1, 1);
            pcells = zeros(ni-1, nj-1, 1);
    
            if k==1
                mcells = cells(:,:,k);
            elseif k == nk
                pcells = cells(:,:,k-1);
            else
                mcells = cells(:, :, k);
                pcells = cells(:, :, k-1);
            end
            for i=1:ni-1
                for j=1:nj-1
                    n = kfaces(i,j,k);
                    faces(n).n1 = nodes(i,j,k);
                    faces(n).n2 = nodes(i+1,j,k);
                    faces(n).n3 = nodes(i+1,j+1,k);
                    faces(n).n4 = nodes(i,j+1,k);
                    faces(n).cr = mcells(i,j,1);
                    faces(n).cl = pcells(i,j,1);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                end
            end
        end
    end

    if iWrite
        
        path
        fid = fopen(path,'w');

        if nk == 1

            fprintf(fid,'(0 "Grid:")\n');
            fprintf(fid,'\n');
        
            fprintf(fid,'(0 "Dimensions:")\n');
            fprintf(fid,'(2 2)\n');
            fprintf(fid,'\n');
            
            fprintf(fid,'(0 "Cells:")\n');
            fprintf(fid,'(12 (0 %X %X 0))\n',[1 (ncells)]);
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Faces:")\n');
            fprintf(fid,'(13 (0 %X %X 0))\n',[1 length(faces)]);
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Nodes:")\n');
            fprintf(fid,'(10 (0 %X %X 0 2))\n',[1 length(nodelist)]);
            fprintf(fid, '\n');
            
            fprintf(fid,'(12 (7 1 %X 1 3))\n',ncells);
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Interior:")\n');
            fprintf(fid,'(13 (2 %X %X %X 2)(\n', [interior_faces(1) interior_faces(2) 2]);
            for n=(outlet_faces(2)+1):length(faces)
                fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Inlet:")\n');
            fprintf(fid,'(13 (3 %X %X %X 2)(\n', [inlet_faces(1) inlet_faces(2) 4]);
            for n=inlet_faces(1):inlet_faces(2)
                fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Inviscid wall:")\n');
            fprintf(fid,'(13 (4 %X %X %X 2)(\n', [invc_faces(1) invc_faces(2) 3]);
            for n=invc_faces(1):invc_faces(2)
                fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
    
            fprintf(fid,'(0 "Outlet:")\n');
            fprintf(fid,'(13 (5 %X %X %X 2)(\n', [outlet_faces(1) outlet_faces(2) 5]);
            for n=outlet_faces(1):outlet_faces(2)
                fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Wall:")\n');
            fprintf(fid,'(13 (6 %X %X %X 2)(\n', [wall_faces(1) wall_faces(2) 3]);
            for n=wall_faces(1):wall_faces(2)
                fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
        
            fprintf(fid,'(10 (1 %X %X 1 2)\n(\n',[1 length(nodelist)]);
            for i=1:length(nodelist)
                fprintf(fid,'%10.8f %10.8f\n', [nodelist(i).x nodelist(i).y]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');

        else

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
            
            fprintf(fid,'(12 (9 1 %X 1 4))\n',ncells-1);
            fprintf(fid, '\n');
%%

            fprintf('Writing faces\n')
            fprintf(fid,'(0 "Interior:")\n');
            fprintf(fid,'(13 (2 %X %X %X 4)(\n', [interior_faces(1) interior_faces(2) 2]);
            for n=(outlet_faces(2)+1):length(faces)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Inlet:")\n');
            fprintf(fid,'(13 (3 %X %X %X 4)(\n', [inlet_faces(1) inlet_faces(2) 4]);
            for n=inlet_faces(1):inlet_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Inviscid wall:")\n');
            fprintf(fid,'(13 (4 %X %X %X 4)(\n', [invc_faces(1) invc_faces(2) 3]);
            for n=invc_faces(1):invc_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
    
            fprintf(fid,'(0 "Outlet:")\n');
            fprintf(fid,'(13 (5 %X %X %X 4)(\n', [outlet_faces(1) outlet_faces(2) 5]);
            for n=outlet_faces(1):outlet_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf(fid,'(0 "Wall:")\n');
            fprintf(fid,'(13 (6 %X %X %X 4)(\n', [wall_faces(1) wall_faces(2) 3]);
            for n=wall_faces(1):wall_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');

            fprintf(fid,'(0 "z=0:")\n');
            fprintf(fid,'(13 (7 %X %X %X 4)(\n', [zm_faces(1) zm_faces(2) 3]);
            for n=wall_faces(1):wall_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');

            fprintf(fid,'(0 "z=span:")\n');
            fprintf(fid,'(13 (8 %X %X %X 4)(\n', [zp_faces(1) zp_faces(2) 3]);
            for n=wall_faces(1):wall_faces(2)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 ...
                    faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        
            fprintf('Writing nodes\n')
            fprintf(fid,'(10 (1 %X %X 1 3)\n(\n',[1 length(nodelist)]);
            for i=1:length(nodelist)
                fprintf(fid,'%10.8f %10.8f %10.8f\n', [nodelist(i).x nodelist(i).y, nodelist(i).z]);
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');

        end
        
        fclose(fid);
    end

