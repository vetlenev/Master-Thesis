classdef UtilFunctions
   methods (Static)
       
       function seed = setSeed(data_dir, my_seed)
            stochastic_files = dir(data_dir);
            for i=1:length(stochastic_files)
                seed_found = regexp(stochastic_files(i).name, string(my_seed), 'match');
                if ~isempty(seed_found)
                    seed = load(strcat(data_dir, '/', strcat('seed_', string(my_seed), '.mat')), 'seed');                
                    seed = seed.seed;
                    disp('Seed already used. Replicating state.')
                    break;
                elseif i == length(stochastic_files) % no runs for this seed yet --> initialize rng
                    seed = rng();
                    seed.Seed = my_seed;
                    save(sprintf(strcat(data_dir, '/', 'seed_%d'), my_seed), 'seed');
                    disp('New seed. Setting new state.')                    
                end
            end
       end
       
       function f = fullsizeFig(k)          
           f = figure(k);
           f.Position = get(0, 'Screensize');           
       end
       
       function sorted_struct = sortStructByField(my_struct, fieldname)
           tab = struct2table(my_struct);
           sorted_tab = sortrows(tab, fieldname);
           sorted_struct = table2struct(sorted_tab);
       end       
       
       function [cutted_cells] = cutLayer(grid, cells, hole_min, hole_max, layer_type, layer_idx, n_cuts)
            holes = [];
            for i=1:n_cuts
                hole_start = min(grid.X) + rand(1)*((1-hole_min)*max(grid.X)-min(grid.X));
                hole_stop = hole_start + max(grid.X)*(hole_min + rand(1)*(hole_max - hole_min)); 

                hole = grid.G.cells.indexMap(grid.X >= hole_start & ...
                                            grid.X <= hole_stop & ...
                                            grid.Z > grid.zStart.(layer_type)(layer_idx) & ...
                                            grid.Z < grid.zStop.(layer_type)(layer_idx));
                holes = cat(1, holes, hole);                
            end

            cutted_cells = grid.G.cells.indexMap(setdiff(cells, holes));
       end
       
       function [boundary_faces, bottom] = localBoundaryFaces(G, cells, varargin)
           % cells: list of indices of cells to extract boundary faces of
            opt = struct('full_dim', false); % defaults: sharp interface, negligable horizontal flux for NVEHorz cells (r=0) 
            opt = merge_options(opt, varargin{:});
            
            all_faces = zeros([numel(cells) 6]); % each cell has 6 faces in 3D
            
            for i=1:size(cells)
                % Only in pocessing of grid before simulation, so despite
                % potentially many loops, not too expensive overall
                icell_start = G.cells.facePos(cells(i));
                icell_stop = G.cells.facePos(cells(i)+1)-1;
                all_faces(i,:) = G.cells.faces(icell_start:icell_stop, 1);
            end
            
            [ii, jj, kk] = gridLogicalIndices(G);
            nx = numel(unique(ii(cells)));                                  
            
            if opt.full_dim
                ny = numel(unique(jj(cells)));
                
                bottom = UtilFunctions.boundingFaces(G, 'bottom', all_faces, kk, @max);
                top = UtilFunctions.boundingFaces(G, 'top', all_faces, kk, @min);
                left = UtilFunctions.boundingFaces(G, 'left', all_faces, ii, @min);
                right = UtilFunctions.boundingFaces(G, 'right', all_faces, ii, @max);
                west = UtilFunctions.boundingFaces(G, 'west', all_faces, jj, @min);
                east = UtilFunctions.boundingFaces(G, 'east', all_faces, jj, @max);
                                
                boundary_faces = cat(1, bottom, top, left, right, west, east);                
            else        
                bottom = all_faces(:,6);       
                bottom = bottom(numel(bottom)-(nx-1):numel(bottom));
                left = all_faces(:,1);
                left = left(1:nx:numel(left));
                right = all_faces(:,2);
                right = right(nx:nx:numel(right));
                top = all_faces(:,5);
                top = top(1:nx);
                boundary_faces = cat(1, bottom, top, left, right); 
            end                                            
       end
       
       function [faces] = boundingFaces(G, side, all_faces, xx, max_or_min)
           side_idx = struct('left', 1, 'right', 2, 'west', 3, ...
                                'east', 4, 'top', 5, 'bottom', 6);
           
            faces = all_faces(:, side_idx.(side));
            n = G.faces.neighbors(faces, :);
            [row, col] = find(n == 0);
            for r=1:numel(row)
               n(row(r), col(r)) = n(row(r), 3-col(r)); % if boundary face, select cell adjacent to boundary but inside domain 
            end
           
            bound_n = any(xx(n) == max_or_min(xx(n)), 2);
            faces = faces(bound_n);
       end
       
       function [Gs, cmap, fmap, nmap] = extractLayerSubgrid(G, Gh, bottom_faces, sealingCells, sealingFaces)
          % Extract subgrid for semi-permeable layer, as part of hybrid
          % model.
          % PARAMS:
          %     G: underlying 3D grid (full domain).
          %     Gh: hybrid grid
          %     bottom_faces: array of faces bounding the bottom of desired
          %                     layer.
          %     ?all_top_faces?: array of faces bounding the top of
          %                   all semi-perm layers
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
            
            for i=1:numel(sealingFaces) % handle overlap for sealing faces
                sealingF = sealingFaces{i};
                sealingC = G.faces.neighbors(sealingF,:);
                [~, col_idx] = max(kk(subcells), [], 2); % choose bottom neighbor, not top
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
                          
            % extract separate grids
            [Gss, cmaps, fmaps, nmaps] = extractSubgrid(G, vecells); % subgrid for VE regions under sealing layer
            [Gsf, cmapf, fmapf, nmapf] = extractSubgrid(G, finecells); % subgrid for fine regions under sealing layer
            
            % extract combined grid
            [Gs_tot, ~, ~, ~] = extractSubgrid(G, subcells);
            
            Gs = {Gss, Gsf, Gs_tot};
            cmap = {cmaps, cmapf};
            fmap = {fmaps, fmapf};
            nmap = {nmaps, nmapf};          
       end
       
       function [var] = initNanADI(ad_var)
           var = ADI(nan(size(ad_var.val)), {nan(size(ad_var.jac{1})), nan(size(ad_var.jac{2})), nan(size(ad_var.jac{3})), ...
                                                        nan(size(ad_var.jac{4})), nan(size(ad_var.jac{5}))});
       end
   end
end