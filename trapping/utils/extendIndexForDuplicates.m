function [c_sub] = extendIndexForDuplicates(c, c_u, c_sub)
    % Extend array for virtual cells connected to multiple semi-permeable
    % layers by duplicating associated cell indices
    % PARAMETERS:
    %   c     - cell indices from hybrid grid, includes duplicates
    %   c_u   - unique cell indices from hybrid grid
    %   c_sub - associated cell indices for c from subgrid
    c_ctf = [c_u, c_sub]; % make matrix mapping from hybrid cell indices to associated local index of subgrid    
    c_ctf1 = c_ctf(:,1);
    c_ctf2 = c_ctf(:,2);

    edges = min(c):max(c)+1; % +1 to include last index
    [counts, values] = histcounts(c, edges);
    duplicates = values(counts > 1);
    
    c_sub = zeros(size(c));
    idx_cH = ~ismember(c, duplicates); % no-duplicate-indexes for cells cH
    cH_no_duplicates = c(idx_cH);
    idx_ctf = ismember(c_ctf1, cH_no_duplicates); % no duplicate indexes for cells ctf
    c_sub(idx_cH) = c_ctf2(idx_ctf); % add no-duplicate elements by vectorization
    
    for i=1:numel(duplicates) % extend cH_sub sequentially with repeated cells
        d = duplicates(i);       
        d_idx_glob = ismember(c, d); % logical index for duplicate
        d_sub = c_ctf2(ismember(c_ctf1, d)); % mapped value
               
        c_sub(d_idx_glob) = d_sub; % insert mapped value into associated position in cH_sub
    end
end