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
            opt = struct('full_dim', false, 'FF', false); % defaults: sharp interface, negligable horizontal flux for NVEHorz cells (r=0) 
            opt = merge_options(opt, varargin{:});
            
            if opt.FF
                all_faces = zeros([numel(cells) 4]);
            else
                all_faces = zeros([numel(cells) 6]); % each cell has 6 faces in 3D
                for i=1:size(cells)
                    % Only in pocessing of grid before simulation, so despite
                    % potentially many loops, not too expensive overall
                    icell_start = G.cells.facePos(cells(i));
                    icell_stop = G.cells.facePos(cells(i)+1)-1;
                    all_faces(i,:) = G.cells.faces(icell_start:icell_stop, 1);
                end
            end                       
            
            if opt.FF
                ii = G.i;
                kk = G.j;
            else
                [ii, jj, kk] = gridLogicalIndices(G);
            end
            nx = numel(unique(ii(cells)));                                  
           
            if opt.FF
                [G_sub, gc, gf] = extractSubgrid(G, cells);
                G_sub.i = ii(gc);
                G_sub.j = kk(gc);
                boundary_faces_sub = boundaryFaces(G_sub);
                boundary_faces = gf(boundary_faces_sub); 
                
                bottom_cells_sub = G_sub.cells.indexMap(G_sub.j == min(G_sub.j));
                %bottom_cells = gc(bottom_cells_sub);
                bottom_cells = bottom_cells_sub;
                bottom = [];
                for i=1:numel(bottom_cells)
                    bface = G_sub.cells.faces(G_sub.cells.facePos(i):G_sub.cells.facePos(i+1)-1,:);
                    bface = intersect(gf(bface), boundary_faces);
                    bottom = [bottom; bface];
                end
                bottom_centroids = G.faces.centroids(bottom, 1);
                [~, min_idx] = min(bottom_centroids);
                [~, max_idx] = max(bottom_centroids);
                bottom([min_idx, max_idx]) = []; % remove face on west+east side



            elseif opt.full_dim
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
       
       function [faces] = boundingFaces(G, side, all_faces, xx, max_or_min, varargin)
           %Extract faces from a list of 'all_faces' that are directed
           %towards specified 'side' and satisfy 'max_or_min'.     
           opt = struct('FF', false); % defaults: sharp interface, negligable horizontal flux for NVEHorz cells (r=0) 
           opt = merge_options(opt, varargin{:});

           if opt.FF
               side_idx = struct('left', 1, 'right', 2, 'west', 3, ...
                                'east', 4, 'top', 5, 'bottom', 6);
           else
               side_idx = struct('left', 1, 'right', 2, 'west', 3, ...
                                'east', 4, 'top', 5, 'bottom', 6);
           end
           
            faces = all_faces(:, side_idx.(side));
            n = G.faces.neighbors(faces, :);
            [row, col] = find(n == 0);
            for r=1:numel(row)
               n(row(r), col(r)) = n(row(r), 3-col(r)); % if boundary face, select cell adjacent to boundary but inside domain 
            end
           
            bound_n = any(xx(n) == max_or_min(xx(n)), 2);
            faces = faces(bound_n);
       end              
       
       function [var] = initNanADI(ad_var)
           var = ADI(nan(size(ad_var.val)), {nan(size(ad_var.jac{1})), nan(size(ad_var.jac{2})), nan(size(ad_var.jac{3})), ...
                                                        nan(size(ad_var.jac{4})), nan(size(ad_var.jac{5}))});
       end
       
       function storeProblemSettings(data, tot_time_, median_time_, dirpath, sim_id, varargin)           
           if nargin > 5
               info = {varargin{1}};
           else
               info = {''};
           end

           if ~isfile(dirpath) % file does not exist -> create table names
               id = {sim_id};
               pe_sealing = data{1};
               pe_rest = data{2};
               perm_sealing = data{3};
               perm_rest = data{4};              
               trans_mult = data{5};               
               swr = data{6};
               snr = data{7};               
               poro = data{8};
               num_years = data{9};
               inj_stop = data{10};
               pv_rate = data{11};
               nx = data{12};
               ny = data{13};
               nz = data{14};
               n_sealing = data{15};
               n_rel = data{16};
               my_seed = data{17};
               max_trap_diff = data{18};
               max_RVE_diff = data{19};
               tot_time = tot_time_;
               median_time = median_time_;               

               T_new = table(id,pe_sealing,pe_rest,perm_sealing,perm_rest,trans_mult, ...
                            swr,snr,poro,num_years,inj_stop,pv_rate, ...
                            nx,ny,nz,n_sealing,n_rel,my_seed,max_trap_diff,max_RVE_diff, ...
                            tot_time,median_time,info);
               writetable(T_new, dirpath);               
           else % file exists -> append to bottom
               T_new = table({sim_id},data{1},data{2},data{3},data{4},data{5}, ...
                            data{6},data{7},data{8},data{9},data{10},data{11}, ...
                            data{12},data{13},data{14},data{15},data{16},data{17}, ...
                            data{18},data{19},tot_time_,median_time_,info);
               T = readtable(dirpath);
               if ~any(strcmp(T.id,string(sim_id))) % only append if unique simulation not stored yet
                   writetable(T_new, dirpath, 'WriteMode', 'Append', 'WriteVariableNames', false);
               end
           end
            
       end

       function storeTrappingData(diff, CO2, dirpath, sim_id, varargin)           
           if nargin > 3
               info = {varargin{1}};
           else
               info = {''};
           end

           if ~isfile(dirpath) % file does not exist -> create table names
               id = {sim_id};
               trap_diff = round(diff.tot/CO2.tot * 100, 2);               
               sr_diff = round(diff.struct_res/CO2.struct_res * 100, 2);
               sr_diff(isinf(sr_diff)) = nan;
               r_diff = round(diff.res/CO2.res * 100, 2);
               r_diff(isinf(r_diff)) = nan;
               fr_diff = round(diff.free_res/CO2.free_res * 100, 2);
               fr_diff(isinf(fr_diff)) = nan;
               sm_diff = round(diff.struct_mob/CO2.struct_mob * 100, 2);
               sm_diff(isinf(sm_diff)) = nan;
               fm_diff = round(diff.free_mob/CO2.free_mob * 100, 2);
               fm_diff(isinf(fm_diff)) = nan;
               exit_diff = round(diff.exit/CO2.exit * 100, 2);
               exit_diff(isinf(exit_diff)) = nan;
                           
               T_new = table(id,trap_diff, sr_diff, r_diff, ...
                                fr_diff, sm_diff, fm_diff, ...
                                exit_diff, info);
               writetable(T_new, dirpath);               
           else % file exists -> append to bottom               
               data = cell(7,1);
               data{1} = round(diff.tot/CO2.tot * 100, 2);
               data{1}(isinf(data{1})) = nan;
               data{2} = round(diff.struct_res/CO2.struct_res * 100, 2);
               data{2}(isinf(data{2})) = nan;
               data{3} = round(diff.res/CO2.res * 100, 2);
               data{3}(isinf(data{3})) = nan;
               data{4} = round(diff.free_res/CO2.free_res * 100, 2);
               data{4}(isinf(data{4})) = nan;
               data{5} = round(diff.struct_mob/CO2.struct_mob * 100, 2);
               data{5}(isinf(data{5})) = nan;
               data{6} = round(diff.free_mob/CO2.free_mob * 100, 2);
               data{6}(isinf(data{6})) = nan;
               data{7} = round(diff.exit/CO2.exit * 100, 2);
               data{7}(isinf(data{7})) = nan;

               T_new = table({sim_id},data{1},data{2},data{3},data{4},data{5}, ...
                              data{6},data{7},info);
               T = readtable(dirpath);
               if ~any(strcmp(T.id,string(sim_id))) % only append if unique simulation not stored yet
                   writetable(T_new, dirpath, 'WriteMode', 'Append', 'WriteVariableNames', false);
               end
           end
            
       end
   end
end