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
       
       function [Gs, cmap, fmap, nmap] = extractLayerSubgrid(G, bottom_faces)
          % Extract subgrid for semi-permeable layer, as part of hybrid
          % model.
          % PARAMS:
          %     G: underlying 3D grid (full domain).
          %     bottom_faces: array of faces bounding the bottom of desired
          %                     layer.
          % RETURNS:
          %     Gs: 3D subgrid upper bounded by given layer
          %     cmap: mapping from cells in subgrid to cells in full grid
          %     fmap: mappnig from faces in subgrid to faces in full grid
          %     nmap: mappnig from nodes in subgrid to nodes in full grid
          [ii, jj, kk] = gridLogicalIndices(G);
            subcells = G.faces.neighbors(bottom_faces,:);
            [~, col_idx] = max(kk(subcells), [], 2); % choose bottom neighbor, not top
            subcells = diag(subcells(:, col_idx));

            subcells_top = min(kk(subcells)); % NB: assuming all layering is parallel to formation top and bottoms
            subcells_left = min(ii(subcells));
            subcells_right = max(ii(subcells));
            subcells_west = min(jj(subcells));
            subcells_east = max(jj(subcells));
            
            subcells = G.cells.indexMap(ii >= subcells_left & ii <= subcells_right & ...
                                        jj >= subcells_west & jj <= subcells_east & ...
                                        kk >= subcells_top);
                                  
            [Gs, cmap, fmap, nmap] = extractSubgrid(G, subcells);
          
       end
       
       function [var] = initNanADI(ad_var)
           var = ADI(nan(size(ad_var.val)), {nan(size(ad_var.jac{1})), nan(size(ad_var.jac{2})), nan(size(ad_var.jac{3})), ...
                                                        nan(size(ad_var.jac{4})), nan(size(ad_var.jac{5}))});
       end
   end
end