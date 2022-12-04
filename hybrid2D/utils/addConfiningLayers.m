function [faces, bottom_faces, sealingCells] = addConfiningLayers(G, varargin)
% Add a single confining layer, represented as faces or cells.
% All cells constituting the confining layer must be direct neighbors.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('type', 'faces', ...
                 'full_dim', false, ...
                 'i_range', [-inf, inf], ...
                 'j_range', [-inf, inf], ...
                 'k_range', [-inf, inf], ...
                 'x_range', [-inf, inf], ...
                 'y_range', [-inf, inf], ...
                 'z_range', [-inf, inf], ...
                 'fraction', 1 ...
                 );

    opt = merge_options(opt, varargin{:});
    [ii, jj, kk] = gridLogicalIndices(G);
    neighbors = G.faces.neighbors;
    act = all(neighbors > 0, 2);
    N = neighbors(act, :);
    all_faces = find(act);  
    
    sealingCells = zeros(G.cells.num, 1); 
    
    x = G.faces.centroids(act, 1);
    y = G.faces.centroids(act, 2);
    z = G.faces.centroids(act, 3);   

    log_mask = all(ii(N) >= opt.i_range(1), 2) & all(ii(N) <= opt.i_range(2), 2) & ...
               all(jj(N) >= opt.j_range(1), 2) & all(jj(N) <= opt.j_range(2), 2) & ...
               all(kk(N) >= opt.k_range(1), 2) & all(kk(N) <= opt.k_range(2), 2);

    coord_mask =  all(x >= opt.x_range(1), 2) & all(x <= opt.x_range(2), 2) & ...
                  all(y >= opt.y_range(1), 2) & all(y <= opt.y_range(2), 2) & ...
                  all(z >= opt.z_range(1), 2) & all(z <= opt.z_range(2), 2);
              
    mask = log_mask & coord_mask;   
    
    %if strcmp(opt.type, 'faces')
        mask = mask & kk(N(:, 1)) ~= kk(N(:, 2)); % omit faces in horizontal direction  
    %end

    faces = all_faces(mask);
    if opt.fraction < 1 % only include a fraction of the sealing faces
        tmp = rand(size(faces));
        faces = faces(tmp < opt.fraction);
    end
    bottom_faces = [];
    
    if strcmp(opt.type, 'cells')
%         z_faces = G.faces.centroids(faces, 3);
%         top_confining_faces = faces(z_faces == min(z_faces));
%         bottom_confining_faces = faces(z_faces == max(z_faces));                      
        
        sealing_N = G.faces.neighbors(faces,:);
        
        sealing_N1 = sealing_N(:,1); % upper neighbor
        sealing_N2 = sealing_N(:,2); % lower neighbor
        
        top_boundary_idx = kk(sealing_N1) == min(kk);
        bottom_boundary_idx = kk(sealing_N2) == max(kk);
        
        top_neighbors = (~ismember(sealing_N1, sealing_N2) + top_boundary_idx) == 1; % top layer of sealing cells, excluding top boundary (confining layer)          
        top_neighbors_idx = sealing_N1(top_neighbors);
        bottom_neighbors = (~ismember(sealing_N2, sealing_N1) + bottom_boundary_idx) == 1; % bottom layer of sealing cells, excluding bottom boundary (confining layer)
        bottom_neighbors_idx = sealing_N2(bottom_neighbors);
        
        if all(abs(kk(sealing_N2) - kk(sealing_N1)) == 1) % unit-thickness face -> make cell one-layered
            sealing_N(bottom_neighbors, 2) = nan; % only remove bottom neighbor, retain layer right above sealing face
        else                
            sealing_N(top_neighbors, 1) = nan; % get correct sealing faces         
            sealing_N(bottom_neighbors, 2) = nan; % remove cell layer below bottom sealing face         
        end
        
        sealing_N = sealing_N(~isnan(sealing_N));
        sealingCells(sealing_N) = 1;  
        sealingCells = logical(sealingCells);  

        % find bounding sealing faces
        sealingCellsIdx = G.cells.indexMap(sealingCells);    

        [faces, bottom_faces] = UtilFunctions.localBoundaryFaces(G, sealingCellsIdx, 'full_dim', opt.full_dim);
    end
        
end
