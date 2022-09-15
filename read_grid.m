function [blk, blocks] = read_grid(case_name)
    blk = {};

    base_folder = cd;
    temp_slash = '/'; if ispc, temp_slash = '\'; end
    path = [base_folder temp_slash case_name];
    blocks = readmatrix([path temp_slash 'blockdims.txt']);
    NB = size(blocks,1);
    blk.blockdims = zeros(NB,3);
    for i=1:NB
        [ni,nj,nk] = feval(@(x) x{:}, num2cell(blocks(i,:)));
        blk.blockdims(i,:) = [ni nj nk];
        grid = readmatrix([path temp_slash 'grid_' num2str(i) '.txt']);
        blk.x{i} = reshape(grid(:,1),[ni,nj]);
        blk.y{i} = reshape(grid(:,2),[ni,nj]);
        blk.nk{i} = nk;
    end
end