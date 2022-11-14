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
       
       function [boundary_faces] = localBoundaryFaces(G, cells, varargin)
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
                
%                 bottom = all_faces(:,6);                
%                 nb = G.faces.neighbors(bottom, :);
%                 nb_i = all(nb > 0, 2);
%                 nb = nb(nb_i);
%                 bottom = bottom(nb_i);
%                 max_nb = any(kk(nb) == max(kk(nb)), 2); % indices of faces neighbors to bottommost cells in sealing layer
%                 bottom = bottom(max_nb);              
%                 %bottom = bottom(numel(bottom)-(nx-1):numel(bottom));
%                 
%                 left = all_faces(:,1);
%                 nl = G.faces.neighbors(left, :);
%                 nl_i = all(nl > 0, 2); % interior neighbors
%                 nl = nl(nl_i);
%                 left = left(nl_i);
%                 min_nl = any(ii(nl) == min(ii(nl)), 2);
%                 left = left(min_nl);
%                 %left = left(1:nx:numel(left));
%                 
%                 right = all_faces(:,2);
%                 nr = G.faces.neighbors(right, :);
%                 nr_i = all(nr > 0, 2);
%                 nr = nr(nr_i);
%                 right = right(nr_i);
%                 max_nr = any(ii(nr) == max(ii(nr)), 2);
%                 right = right(max_nr);
%                 %right = right(nx:nx:numel(right));                
%                 
%                 top = all_faces(:,5);
%                 nt = G.faces.neighbors(top, :);
%                 nt_i = all(nt > 0, 2);
%                 nt = nt(nt_i);
%                 top = top(nt_i);
%                 min_nt = any(kk(nt) == min(kk(nt)), 2);
%                 top = top(min_nt);
%                 %top = top(1:nx);
%                 
%                 west = all_faces(:,3);
%                 nw = G.faces.neighbors(west, :);
%                 nw_i = all(nw > 0, 2);
%                 nw = nw(nw_i);
%                 west = west(nw_i);
%                 min_nw = any(jj(nw) == min(jj(nw)), 2);
%                 west = west(min_nw);
%                 
%                 east = all_faces(:,4);
%                 ne = G.faces.neighbors(east, :);
%                 ne_i = all(ne > 0, 2);
%                 ne = ne(ne_i);
%                 east = east(ne_i);
%                 max_ne = any(jj(ne) == max(jj(ne)), 2);
%                 east = east(max_ne);
                
%                 west_old = all_faces(:,3); 
%                 west_single_size = numel(1:(nx*(ny-1)):numel(west_old));
%                 west = zeros(nx, west_single_size);
%                 for i=1:nx
%                     west(i,:) = west_old(i:(nx*(ny-1)):numel(west_old));
%                 end                            
%                 
%                 east_old = all_faces(:,4);
%                 east_single_size = numel((nx*(ny-1)+1):(nx*(ny-1)):numel(east_old));
%                 east = zeros(nx, east_single_size);
%                 for i=1:nx
%                     east(i,:) = east_old((nx*(ny-1)+i):(nx*(ny-1)):numel(east_old));
%                 end                
                
                % Flatten back to original ordering
%                 west = west(:);
%                 east = east(:); 
                
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
               n(row(r), col(r)) = n(row(r), 3-col(r));  
            end
            %n_i = all(n > 0, 2); % interior faces 
            %n = n(n_i);
            %faces = faces(n_i);
            bound_n = any(xx(n) == max_or_min(xx(n)), 2);
            faces = faces(bound_n);
       end
       
       function [var] = initNanADI(ad_var)
           var = ADI(nan(size(ad_var.val)), {nan(size(ad_var.jac{1})), nan(size(ad_var.jac{2})), nan(size(ad_var.jac{3})), ...
                                                        nan(size(ad_var.jac{4})), nan(size(ad_var.jac{5}))});
       end
   end
end