classdef AnticlineSetupPEBI < handle
    %AnticlineSetupPEBI: Set up grid of straight and anticline forms
    %conforming to a PEBI structure.
    properties
        G
        X
        Z
        nLowperm
        nImperm
        ii
        kk
        dX
        dZ
        angle
        xStart
        xStop
        zStart
        zStop
        corr_len_x
        corr_len_z
        Perm
        lowperm
        baseperm
    end
    
    methods
        function o = AnticlineSetupPEBI(G, n_lowperm_layers, n_imperm_layers, perms)
            [ii, ~, kk] = gridLogicalIndices(G);            
            
            if nargin == 4
                o.G = G;
                o.X = G.cells.centroids(:,1);
                o.Z = G.cells.centroids(:,3);
                o.nLowperm = n_lowperm_layers;
                o.nImperm = n_imperm_layers;
                o.ii = ii;
                o.kk = kk;
                o.dX = mode(diff(o.X));
                o.dZ = mode(diff(o.Z));
                
                o.angle = pi/4;
                o.xStart = struct; o.xStop = struct;
                o.zStart = struct; o.zStop = struct;               
                o.corr_len_x = 1; o.corr_len_z = 1;
                           
                o.lowperm = perms{1};
                o.baseperm = perms{2};
                o.Perm = repmat(o.baseperm, [o.G.cells.num 1]);
            end
        end
        
        function [line_idx, anticline_idx] = DefineLayers(o)
            % Randomly distribute bunch of low/imperm layers           
            
            perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);         

            anticline_idx = 1:3:(o.nLowperm+o.nImperm); % every third layer is anticline
            line_idx = setdiff(1:(o.nLowperm+o.nImperm), anticline_idx); % indices for line layers
            
            rand_x_mid = linspace(min([o.X])-o.dX/2, max(o.X)+o.dX/2, o.nLowperm+o.nImperm); % o.dX/2 to move from centroid to endpoint
            rand_x_mid = rand_x_mid(randperm(length(rand_x_mid))); % random permute

            o.xStart.lowperm = rand_x_mid(1:o.nLowperm);
            o.xStart.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm);

            x_length = randi([round(perc_val(o.X, 0.05)), ...
                                                round(perc_val(o.X, 0.35))], ...
                                                [o.nLowperm+o.nImperm, 1]).';

            o.xStart.lowperm = rand_x_mid(1:o.nLowperm) - round(x_length(1:o.nLowperm)/2);
            o.xStop.lowperm = rand_x_mid(1:o.nLowperm) + round(x_length(1:o.nLowperm)/2);

            o.xStart.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm) - round(x_length(o.nLowperm+1:o.nLowperm+o.nImperm)/2);
            o.xStop.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm) + round(x_length(o.nLowperm+1:o.nLowperm+o.nImperm)/2);           

            rand_z_start = linspace(perc_val(o.Z, 0.1), perc_val(o.Z, 0.95), o.nLowperm+o.nImperm);
            o.zStart.lowperm = sort(randsample(rand_z_start, o.nLowperm));
            rand_z_start(ismember(rand_z_start, o.zStart.lowperm)) = [];
            o.zStart.imperm = sort(rand_z_start);

            o.zStop.lowperm = o.zStart.lowperm + perc_val(o.Z, 0.02);
            o.zStop.imperm = o.zStart.imperm + perc_val(o.Z, 0.02); 
            
            o.corr_len_x = mean(o.xStop.lowperm - o.xStart.lowperm) / 100; % default value
            o.corr_len_z = mean(o.zStop.lowperm - o.zStart.lowperm) / 10; % default value
            
        end
        
        function [PG, constraints, added_layers, trapped_cells] = ...
                GenerateLayers(o, line_idx, anticline_idx) 
            % Include defined layers in grid and set their physical
            % properties
            perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);
            
            added_layers = {[], []};
            layer_types = {'imperm', 'lowperm'};

            anticline_idx_lowperm = anticline_idx(1:round(o.nLowperm/(o.nLowperm+o.nImperm)*length(anticline_idx)));
            anticline_idx_lowperm(anticline_idx_lowperm > o.nLowperm) = []; % remove invalid values
            anticline_idx_imperm = anticline_idx(find(anticline_idx == anticline_idx_lowperm(end))+1:end) - o.nLowperm;
            anticline_idx_imperm(anticline_idx_imperm > o.nImperm) = []; % remove invalid values

            anticline_idxs = {anticline_idx_imperm, anticline_idx_lowperm};

            line_idx_lowperm = line_idx(1:round(o.nLowperm/(o.nLowperm+o.nImperm)*length(line_idx)));
            line_idx_lowperm(line_idx_lowperm > o.nLowperm) = []; % remove invalid values
            line_idx_imperm = line_idx(find(line_idx == line_idx_lowperm(end))+1:end) - o.nLowperm;
            line_idx_imperm(line_idx_imperm > o.nImperm) = []; % remove invalid values

            line_idxs = {line_idx_imperm, line_idx_lowperm};            
            
            o.xStart.imperm(anticline_idxs{1}) = o.xStart.imperm(anticline_idxs{1}) - perc_val(o.X, 0.05);
            o.xStop.imperm(anticline_idxs{1}) = o.xStop.imperm(anticline_idxs{1}) + perc_val(o.X, 0.05);
            o.xStart.lowperm(anticline_idxs{2}) = o.xStart.lowperm(anticline_idxs{2}) - perc_val(o.X, 0.05);
            o.xStop.lowperm(anticline_idxs{2}) = o.xStop.lowperm(anticline_idxs{2}) + perc_val(o.X, 0.05);           
            
            o.zStop.imperm(anticline_idxs{1}) = o.zStop.imperm(anticline_idxs{1}) + perc_val(o.Z, 0.01);
            o.zStop.lowperm(anticline_idxs{2}) = o.zStop.lowperm(anticline_idxs{2}) + perc_val(o.Z, 0.01);            
                       
            perm = repmat(o.baseperm, [o.G.cells.num 1]);
            perm_vals = {1e-5*milli*darcy, o.lowperm};                        

            trapped_cells = {{}, {}}; % {structural imperm, structural free}
            sin_arg = @(x, x_start, x_stop) 0.8*(x - x_start)/(x_stop - x_start) + ...
                                             0.2*(x_stop - x)/(x_stop - x_start);
            constraints = {};

            for k=1:numel(layer_types)
                layer_type = layer_types{k};               
                % Make straight forms
                for i=line_idxs{k}                    
                                                        
                    line_cells = o.G.cells.indexMap(o.X >= o.xStart.(layer_type)(i) & ...
                                                    o.X <= o.xStop.(layer_type)(i) & ...
                                                    o.Z > o.zStart.(layer_type)(i) & ...
                                                    o.Z < o.zStop.(layer_type)(i));

                    max_z = max(o.kk(line_cells));
                    z_bottom = unique(o.Z(o.kk == max_z + 1));

                    [ix, ~, kz] = gridLogicalIndices(o.G, line_cells);
                    nxi = numel(unique(ix)); nzi = numel(unique(kz)); % dimensions of particular low-perm layer

                    line_cells_G = intersect(o.G.cells.indexMap, line_cells);
                    added_layers_G = intersect(o.G.cells.indexMap, vertcat(added_layers{:})); 
                    overlapped_cells = o.G.cells.indexMap(intersect(line_cells_G, ...
                                                                    added_layers_G));
                    store_perm = perm(overlapped_cells);
                    
                    perm(line_cells) = perm_vals{k} + perm_vals{k}*FastGaussian([nxi nzi], 0.8, [o.corr_len_x o.corr_len_z]);
                    perm(overlapped_cells) = store_perm;  
                    
                    added_layers{k} = cat(1, added_layers{k}, line_cells);
                end

                % Make anticline forms
                for i=anticline_idxs{k}
                    num_x_cells = numel(o.G.cells.indexMap(o.X > o.xStart.(layer_type)(i) ... 
                                        & o.X < o.xStop.(layer_type)(i) ...
                                        & o.Z == min(o.Z))); % only select ONE horizontal patch
                    theta0 = o.angle;
                    theta = linspace(theta0, pi - theta0, num_x_cells);
                    r = (o.xStop.(layer_type)(i) - o.xStart.(layer_type)(i))/2;
                    z_anticline_start = -r*sin(theta) + o.zStart.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % -r*sin to get anticline (since z positive downwards)
                    z_anticline_stop = -r*sin(theta) + o.zStop.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % + r*cos(pi/4) to shift to bottom of spherical cap

                    anticline_cells = [];
                    anticline_cells_dummy = o.G.cells.indexMap(o.X > o.xStart.(layer_type)(i) & ...
                                                    o.X < o.xStop.(layer_type)(i) & ...
                                                    o.Z > min(z_anticline_start) & ...
                                                    o.Z < max(z_anticline_stop));

                    max_z = max(o.kk(anticline_cells_dummy));
                    min_z = min(o.kk(anticline_cells_dummy)); 
                    z_bottom = unique(o.Z(o.kk == max_z)); % z-coord bottom part of layer                    

                    
                    x_anticline = linspace(o.xStart.(layer_type)(i), o.xStop.(layer_type)(i), num_x_cells);
                    amp = sin_arg(x_anticline, o.xStart.(layer_type)(i), o.xStop.(layer_type)(i));
                    z_anticline = o.zStart.(layer_type)(i) + r*sin(pi*amp);
                    constraint_top = [x_anticline(:), z_anticline(:)];
                    
                    z_anticline = o.zStop.(layer_type)(i) + r*sin(pi*amp);
                    constraint_bottom = {x_anticline(:), z_anticline(:)};
                    %disp(constraint_top)
                    
                    constraints = cat(2, constraints, constraint_top);                                     
                    
                    trapped_slice = o.G.cells.indexMap(o.X > o.xStart.(layer_type)(i) & ...
                                                     o.X < o.xStop.(layer_type)(i) & ...
                                                     o.Z > z_anticline & ...
                                                     o.Z <= z_bottom).';  
                                                     
                    trapped_cells{k} = cat(1, trapped_cells{k}, trapped_slice);

                    num_z_cells = max_z-min_z+1; % +1 to get correct nr vertical cells in anticline form
                    
                    added_layers_G = intersect(intersect(o.G.cells.indexMap, anticline_cells_dummy), ...
                                                intersect(o.G.cells.indexMap, vertcat(added_layers{:}))); 
              
                    overlapped_cells = o.G.cells.indexMap(added_layers_G);
                    
                    store_perm = perm(overlapped_cells);
            
                    perm(anticline_cells_dummy) = perm_vals{k} + perm_vals{k}*FastGaussian([num_x_cells num_z_cells], 0.8, [o.corr_len_x o.corr_len_z]); % rectangular box around anticline
                    perm(overlapped_cells) = store_perm;
                    background_cells = o.G.cells.indexMap(setdiff(anticline_cells_dummy, ...
                                                        vertcat(anticline_cells, vertcat(added_layers{:}))));
                    perm(background_cells) = o.baseperm; % reset perm for parts of rectangle NOT part of anticline
                    added_layers{k} = cat(1, added_layers{k}, anticline_cells);
                end
                
                trapped_cells{k} = horzcat(trapped_cells{k}{:});
            end
            
            %F = surfaceSites2D(constraints, r/10, 'interpolateFC', [true, true])
            nx = o.G.cartDims(1); nz = o.G.cartDims(3);
            disp(constraints{1})
            
            PG = pebiGrid2D(1, [max(o.X)+o.dX, max(o.Z)+o.dZ], ...
                                       'faceConstraints', constraints);          
            %all_added_layers = vertcat(added_layers{:});            
            perm = perm.*sign(perm); % set all positive values
            perm = max(perm, perm_vals{1}); % perm must be positive -> cap at non-zero value to avoid singular matrix
            o.Perm = perm;
        end
    end
end

