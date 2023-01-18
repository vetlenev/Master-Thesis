function [Gs, cmap, fmap, nmap] = ExtractLayerSubgrid(G, Gh, bottom_faces, sealingCells, sealingFaces)
      % Extract subgrid for semi-permeable layer, as part of hybrid
      % model.
      % PARAMS:
      %     G: underlying 3D grid (full domain).
      %     Gh: hybrid grid (NB: not needed!)
      %     bottom_faces: array of faces bounding the bottom of desired
      %                     layer.
      %     sealingCells    - all cell representations of semi-perm layers
      %     sealingFaces    - all face representations of semi-perm layers
      % RETURNS:
      %     Gs: 3D subgrid upper bounded by given layer
      %     cmap: mapping from cells in subgrid to cells in full grid
      %     fmap: mappnig from faces in subgrid to faces in full grid
      %     nmap: mappnig from nodes in subgrid to nodes in full grid
        [ii, jj, kk] = gridLogicalIndices(G);
        p = Gh.partition;
        subcells = G.faces.neighbors(bottom_faces,:);

        if ~any(all(subcells, 2)) % remove neighbors outside domain
           subcells = unique(subcells(subcells ~= 0));
        end

        [~, col_idx] = max(kk(subcells), [], 2); % choose bottom neighbor, not top
        subcells = diag(subcells(:, col_idx));

        iis = ii(subcells);
        jjs = jj(subcells);
        kks = kk(subcells);

        subcells_top = min(kks); % NB: assuming all layering is parallel to formation top and bottoms
        subcells_left = min(iis);
        subcells_right = max(iis);
        subcells_west = min(jjs);
        subcells_east = max(jjs);                      

        subcells = G.cells.indexMap(ii >= subcells_left & ii <= subcells_right & ...
                                    jj >= subcells_west & jj <= subcells_east & ...
                                    kk >= subcells_top);


        remove_cells = {};
        for i=1:numel(sealingCells) % loop through semi-perm layers and find overlaps and "subsubgrids" to remove from subgrid
            sealingC = G.cells.indexMap(sealingCells{i});
            sealing_overlap = sealingC(ismember(sealingC, subcells));
            sealing_top = min(kk(sealing_overlap));
            sealing_left = min(ii(sealing_overlap));
            sealing_right = max(ii(sealing_overlap));
            sealing_west = min(jj(sealing_overlap));
            sealing_east = max(jj(sealing_overlap));

            if ~isempty(sealing_overlap)
                remove_c = G.cells.indexMap(ii >= sealing_left & ii <= sealing_right & ...
                                                jj >= sealing_west & jj <= sealing_east & ...
                                                kk >= sealing_top);
            else
                remove_c = [];
            end
            remove_cells = cat(1, remove_cells, remove_c);
        end

        for i=1:numel(sealingFaces) % handle overlap for sealing faces
            if isequal(sealingFaces{i}, bottom_faces)
               continue; % avoid checking the same sealing layer 
            end
            sealingF = sealingFaces{i};
            sealingC = G.faces.neighbors(sealingF,:);
            [~, col_idx] = max(kk(sealingC), [], 2); % choose bottom neighbor, not top
            sealingC = diag(sealingC(:, col_idx)); % choose from col-index corresponding to bottom neighbor

            sealing_overlap = sealingC(ismember(sealingC, subcells));
            sealing_top = min(kk(sealing_overlap)); % last term necessary in case of no overlap
            sealing_left = min(ii(sealing_overlap));
            sealing_right = max(ii(sealing_overlap));
            sealing_west = min(jj(sealing_overlap)); 
            sealing_east = max(jj(sealing_overlap));

            if ~isempty(sealing_overlap)
                remove_c = G.cells.indexMap(ii >= sealing_left & ii <= sealing_right & ...
                                                jj >= sealing_west & jj <= sealing_east & ...
                                                kk >= sealing_top);
            else
                remove_c = [];
            end
            remove_cells = cat(1, remove_cells, remove_c);
        end
        % remove other semi-perm layers inside subgrid
        remove_cells = unique(vertcat(remove_cells{:})); % blocks of cells to remove may overlap -> duplicates removed by selecting unique            

        subcells = setdiff(subcells, remove_cells);
        ve_mask = Gh.cells.discretization(p(subcells)) > 1;
        vecells = subcells(ve_mask); % only ve
        finecells = subcells(~ve_mask); % only fine            

%             % extract separate grids
%             [Gss, cmaps, fmaps, nmaps] = extractSubgrid(G, vecells); % subgrid for VE regions under sealing layer
%             [Gsf, cmapf, fmapf, nmapf] = extractSubgrid(G, finecells); % subgrid for fine regions under sealing layer
%             
%             Gs = {Gss, Gsf, Gs_tot};
%             cmap = {cmaps, cmapf};
%             fmap = {fmaps, fmapf};
%             nmap = {nmaps, nmapf};
        % extract combined grid
        [Gs, cmap, fmap, nmap] = extractSubgrid(G, subcells);  
        test = 0;
   end