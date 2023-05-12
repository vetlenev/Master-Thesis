function [extrac, extrac_f] = distributeFineRegions(G, asc, asc_f)
    % Distribute additional fine cells around sealing barriers.
    % Purpose: reduce errors from RVE assumption
    % INPUTS:
    %   G: full-dimensional grid
    %   asc: cell array of sealing barriers cells
    %   asc_f: cell array of faces bounding the sealing barriers
    % RETURNS:
    %   extrac: cell array of additional fine cells
    %   extrac_f: cell array of additional bounding faces
    
    [ii, ~, kk] = gridLogicalIndices(G);
    %N = G.cartDims;
    sealing_barriers = asc;
    extrac = {}; % extra fine cells
    extrac_f = {}; % faces of extra fine cells

    for i=1:numel(sealing_barriers)
        sb = sealing_barriers{i};
        if isempty(find(sb,1))
            continue;
        end
        %diff_i = max(ii(sb)) - min(ii(sb));       
        diff_i = max(ii);
        %diff_k = max(kk(sb)) - min(kk(sb));
        diff_k = max(kk);
        i_stop = min(ii(sb)) + ceil(0.05*diff_i);
        i_start = min(ii(sb)) - ceil(0.05*diff_i);%(i_stop - min(ii(sb)));
        k_stop = max(kk(sb)) + 0.02*diff_k;% + (diff_k == 0);
        k_start = min(kk(sb)) - 0.02*diff_k;% - (diff_k == 0);

        [sealing_faces, sealing_bottom, sealingCells] = addConfiningLayers(G, 'type', 'cells', 'i_range', [i_start, i_stop], 'k_range', [k_start, k_stop+1]);       
        extrac_f = cat(2, extrac_f, sealing_faces);
        extrac = cat(2, extrac, sealingCells);

        i_stop = max(ii(sb)) + ceil(0.05*diff_i);
        i_start = max(ii(sb)) - ceil(0.05*diff_i);%(i_stop - min(ii(sb)));

        [sealing_faces, sealing_bottom, sealingCells] = addConfiningLayers(G, 'type', 'cells', 'i_range', [i_start, i_stop], 'k_range', [k_start, k_stop+1]);       
        extrac_f = cat(2, extrac_f, sealing_faces);
        extrac = cat(2, extrac, sealingCells);
       
    end
end