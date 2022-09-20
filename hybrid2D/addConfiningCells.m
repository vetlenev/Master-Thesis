function [faces_boundary, sealingCells] = addConfiningCells(G, varargin)
% Add sealing layers to grid, assuming non-zero thickness.
% Modelled as fine cells in a hybrid VE model.

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

    opt = struct('type', 'horizontal', ...
                 'i_range', [-inf, inf], ...
                 'j_range', [-inf, inf], ...
                 'k_range', [-inf, inf], ...
                 'x_range', [-inf, inf], ...
                 'y_range', [-inf, inf], ...
                 'z_range', [-inf, inf], ...
                 'thickness', 1, ...
                 'fraction', 1 ...
                 );

    opt = merge_options(opt, varargin{:});
    
    [ii, jj, kk] = gridLogicalIndices(G);
    neighbors = G.faces.neighbors;
    act = all(neighbors > 0, 2);
    N = neighbors(act, :);
    all_faces = find(act);
    
    sealingCells = zeros(G.cells.num, 1);   
    
    if strcmp(opt.type, 'sloped')    
        m = (opt.k_range(2) - opt.k_range(1)) / (opt.i_range(2) - opt.i_range(1)); % slope

        zk = @(xi, t) m*(xi-opt.i_range(1)) + opt.k_range(1) + t;

        sgn = sign(opt.k_range(2) - opt.k_range(1)); % sign of slope
        delta_z = opt.thickness;

        mask = all(ii(N) >= opt.i_range(1), 2) & all(ii(N) <= opt.i_range(2), 2) & ...
                   all(jj(N) >= opt.j_range(1), 2) & all(jj(N) <= opt.j_range(2), 2) & ...
                   all(sgn*kk(N) >= sgn*(opt.k_range(1)-delta_z/2), 2) & all(sgn*kk(N) <= sgn*(opt.k_range(2)+delta_z/2), 2); % within minimum rectangle bounding sloping layer

        mask = mask & all(kk(N) >= zk(ii(N), -delta_z/2), 2) ...
                    & all(kk(N) <= zk(ii(N), delta_z/2), 2);
                
        N1 = N(:,1);
        N2 = N(:,2);
        n_idx = ii(N1) ~= min(ii(N1)) & ii(N2) ~= min(ii(N2));                
        n1 = N1(n_idx);
        n2 = N2(n_idx);
        
        kkk = kk;
        kkk(N1(~n_idx)) = kkk(N2(~n_idx));
        %n12 = ismember(N1, n1) & ismember(N2, n2);
        mask = mask & kk(N1) ~= kk(N2); % omit faces in horizontal direction
        
    elseif strcmp(opt.type, 'horizontal')
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
        %mask = mask & kk(N(:, 1)) ~= kk(N(:, 2)); % omit faces in horizontal direction

    else
        error('%s is not a valid type of sealing face. Valid options are "sloped" and "horizontal".', opt.type);
    end

    faces = all_faces(mask);
    if opt.fraction < 1 % only include a fraction of the sealing faces
        tmp = rand(size(faces));
        faces = faces(tmp < opt.fraction);
    end
    
    % find cells bounded by sealing faces
    sealing_N = G.faces.neighbors(faces, :);
     
    sealing_N1 = sealing_N(:,1); % upper neighbor
    sealing_N2 = sealing_N(:,2); % lower neighbor
     
    top_neighbors = ~ismember(sealing_N1, sealing_N2); % top layer of sealing cells          
    top_neighbors_idx = sealing_N1(top_neighbors);
    bottom_neighbors = ~ismember(sealing_N2, sealing_N1); % bottom layer of sealing cells
    bottom_neighbors_idx = sealing_N2(bottom_neighbors);
    sealing_N(top_neighbors | bottom_neighbors, :) = []; % get correct sealing faces         
    
    sealingCells(sealing_N) = 1;  
    sealingCells = logical(sealingCells);  
    
    % find bounding sealing faces
    sealingCellsIdx = G.cells.indexMap(sealingCells);    
    
    top_c2f = gridCellFaces(G, top_neighbors_idx);
    bottom_c2f = gridCellFaces(G, bottom_neighbors_idx);
    sealing_c2f = gridCellFaces(G, sealingCellsIdx);
    
    % faces at boundary between sealing cells and VE cells
    faces_boundary = intersect(sealing_c2f, [top_c2f; bottom_c2f]);  
end
