function [Gs, cmap, fmap, nmap] = ExtractLayerSubgrid_FF(Gf, Gh, bottom_faces, pidx, sealingCells, sealingFaces)
      % Extract subgrid for a layer of the FluidFlower rig, 
      % as part of hybrid model.
      % PARAMS:
      %     Gf: full-dimensional grid
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
        Gs = struct;
        cmap = struct; fmap = struct; nmap = struct;
        p_num = 1;

        ii = Gf.i;
        jj = Gf.j;        

        if isempty(bottom_faces)
            p_idx = strcat(pidx, '_', string(p_num));
            Gs.(p_idx) = [];
            cmap.(p_idx) = []; fmap.(p_idx) = []; nmap.(p_idx) = [];
            return;
        end
        subcells = Gf.faces.neighbors(bottom_faces,:);
        subcells_sep = {};
        k = 1;

        if ~all(subcells(:)) % remove neighbors outside domain
           ob_cells = ~all(subcells, 2); % out of bound cells
           ob_cells_idx = find(ob_cells);
           subcells_sep{k} = subcells(1:ob_cells_idx(1)-1, :);
           subcells_sep{k+1} = subcells(ob_cells_idx(end)+1:end, :);
           k = k+2;
        else
            subcells_sep{k} = subcells;
            k = k+1;
        end

        nan_cells = any(isnan(jj(subcells_sep{k-1})), 2); % fault cells
        if any(nan_cells)
            nan_cells_idx = find(nan_cells);
            subcells_dummy = subcells_sep{k-1};
            subcells_sep{k-1} = subcells_dummy(1:nan_cells_idx(1)-1, :);
            subcells_sep{k} = subcells_dummy(nan_cells_idx(end)+1:end, :);   
        end
      
        % Split top surface into multiple sub-surfaces
        for p_num=1:numel(subcells_sep)
            subcells = subcells_sep{p_num};
            p_idx = strcat(pidx, '_', string(p_num));

            [~, col_idx] = max(jj(subcells), [], 2); % choose bottom neighbor, not top
            subcells = diag(subcells(:, col_idx));
    
            iis = ii(subcells);
            jjs = jj(subcells);
    
            subcells_top = min(jjs); % NB: assuming all layering is parallel to formation top and bottoms    
            subcells_west = min(iis);
            subcells_east = max(iis);                      
    
            subcells = Gf.cells.indexMap(ii >= subcells_west & ii <= subcells_east & ...
                                        jj >= subcells_top);
    
    
            remove_cells = {};
            for i=1:numel(sealingCells) % loop through semi-perm layers and find overlaps and "subsubgrids" to remove from subgrid
                sealingC = Gf.cells.indexMap(sealingCells{i});
                sealing_overlap = sealingC(ismember(sealingC, subcells));
                sealing_top = min(jj(sealing_overlap));
                sealing_west = min(ii(sealing_overlap));
                sealing_east = max(ii(sealing_overlap));
    
                if ~isempty(sealing_overlap)
                    remove_c = Gf.cells.indexMap(ii >= sealing_west & ii <= sealing_east & ... 
                                                    jj >= sealing_top);
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
                sealingC = Gf.faces.neighbors(sealingF,:);
                if i > 26
                    test = 0;
                end
                % remove cells out of bounds
                c_int = ~any(sealingC == 0, 2);
                sealingC = sealingC(c_int, :);
                jjC = jj(sealingC);
                % remove nan cells
                c_nonan = ~any(isnan(jjC), 2);
                jjC = jjC(c_nonan, :);
                sealingC_nonan = sealingC(c_nonan, :);
                
                [~, col_idx] = max(jjC, [], 2); % choose bottom neighbor, not top
                sealingC = diag(sealingC_nonan(:, col_idx)); % choose from col-index corresponding to bottom neighbor
    
                sealing_overlap = sealingC(ismember(sealingC, subcells));
                sealing_top = min(jj(sealing_overlap)); % last term necessary in case of no overlap
                sealing_west = min(ii(sealing_overlap)); 
                sealing_east = max(ii(sealing_overlap));
    
                if ~isempty(sealing_overlap)
                    remove_c = Gf.cells.indexMap(ii >= sealing_west & ii <= sealing_east & ...
                                                    jj >= sealing_top);
                else
                    remove_c = [];
                end
                remove_cells = cat(1, remove_cells, remove_c);
            end
            % remove other semi-perm layers inside subgrid
            remove_cells = unique(vertcat(remove_cells{:})); % blocks of cells to remove may overlap -> duplicates removed by selecting unique            
    
            subcells = setdiff(subcells, remove_cells);
            
            % extract combined grid
            [Gs_i, cmap_i, fmap_i, nmap_i] = extractSubgrid(Gf, subcells);  
            Gs_i.i = ii; % ii(Gs_i.cells.global);
            Gs_i.j = jj; % jj(Gs_i.cells.global);
            Gs.(p_idx) = Gs_i;            
            cmap.(p_idx) = cmap_i;
            fmap.(p_idx) = fmap_i;
            nmap.(p_idx) = nmap_i;
        end
   end