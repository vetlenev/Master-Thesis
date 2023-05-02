classdef UtilFunctionsFF
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

       function rand_pt = randomPointOnLine(p1, p2, rng_type)
            % Create a random point on line segment from p1 to p2
            %   Detailed explanation goes here
            r = max(0, min(1, 0.2*rand(1) + 0.5));
            if strcmp(rng_type, 'uniform') 
                rand_pt = r*p1 + (1-r)*p2;
            elseif strcmp(rng_type, 'nonuniform')    
                rand_pt = sqrt(r)*p1 + (1-sqrt(r))*p2;
            end
            
       end

       function array_new = removeOverlappingElements(array)
        % Remove any elements occuring more than once from array.
            array_new = sort(array);
            idx_nonunique = diff(array_new) == 0;
            array_remove = array_new(idx_nonunique);
            array_new(ismember(array_new, array_remove)) = [];            
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
                
                bottom = UtilFunctionsFF.boundingFaces(G, 'bottom', all_faces, kk, @max);
                top = UtilFunctionsFF.boundingFaces(G, 'top', all_faces, kk, @min);
                left = UtilFunctionsFF.boundingFaces(G, 'left', all_faces, ii, @min);
                right = UtilFunctionsFF.boundingFaces(G, 'right', all_faces, ii, @max);
                west = UtilFunctionsFF.boundingFaces(G, 'west', all_faces, jj, @min);
                east = UtilFunctionsFF.boundingFaces(G, 'east', all_faces, jj, @max);
                                
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
       
      
   end
end