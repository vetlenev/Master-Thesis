function [cB_idx, cB_sub, ...
          cH_idx, cH_sub, ...
          cBH_idx, cBH_sub] = getRelaxedVEIndices(c_ve, cB, cH, cBH)
% Extract indices from all relazed VE columns, where _idx has one
% element per relaxed VE cell and _sub has one element per cell in
% subgrid.
% PARAMETERS:
%   c_ve - RVE cells (from global grid) to get indices from
%   cB   - bottom sealing layers
%   cH   - horizontal sealing layers
%   cBH  - bottom + horizontal sealing layers

    cB_idx = ismember(cB, c_ve); % index for arrays with one element per bottom cell cB (from hybrid indexes for full grid)
    cB_subidx = ismember(c_ve, cB); % index for arrays with one element per cell in subgrid (from hybrid indexes for full grid)    
    cB_sub = find(cB_subidx); % global cell index
           
    cH_idx = ismember(cH, c_ve); % with duplicates    
    cH_subidx = ismember(c_ve, cH);
    cH = cH(cH_idx);
    cH_sub = find(cH_subidx); % indices for cHorz cells in subgrid
    [cH_u, ~, uHidx] = unique(cH); 
    
    if ~isempty(cH)
        % For cells with multiple rising plumes we need to duplicate the
        % index of the cell in the global list (c_ve)
        cH_sub = extendIndexForDuplicates(cH, cH_u, cH_sub);
    end
        
    cBH_idx = ismember(cBH, c_ve);
    cBH_subidx = ismember(c_ve, cBH);
    cBH = cBH(cBH_idx); 
    cBH_sub = find(cBH_subidx);
    [cBH_u, ~, uBHidx] = unique(cBH);
    
    if ~isempty(cBH)
        cBH_sub = extendIndexForDuplicates(cBH, cBH_u, cBH_sub);
    end
end