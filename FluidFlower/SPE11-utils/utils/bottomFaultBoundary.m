function [poly_obj, external_nodes] = bottomFaultBoundary(all_polys, poly_obj, poly_num, poly_neighbors, ...
                                                                  nodes_overlap, pts_overlap, G_glob, varargin)  
   
% Extract external boundaries of bottom fault
    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,1));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);
    
    poly_glob = Faults(all_polys, 1, 1, true); % just assign a dummy poly number
    poly_glob.p_idx = 'pBFPEBI';
    poly_obj.pBFPEBI = poly_glob;   

    % First: make Fault object for internal faults
    for p=1:numel(poly_num)
        poly_n = poly_num(p);
        poly = Faults(all_polys, poly_n, 1, false); 
        p_idx_pebi = strcat('p', string(poly_n), 'PEBI');
        poly.p_idx = p_idx_pebi;                
    
       poly_obj.(p_idx_pebi) = poly;
    end
   
    all_nodes = [];
    all_pts = [];
    % 1. Assign external boundary NODES
    for i=1:size(poly_neighbors, 1)
        pn = poly_neighbors{i,1}; % neighboring polygon instance        
        pn_face = poly_neighbors{i,2}; % intersecting face of neighbor
        p_idx = strcat(poly_neighbors{i,3}, 'PEBI'); % index of associated polygon in fault
        poly = poly_obj.(p_idx);

        [nodes, pts] = findOverlappingNodes(poly, pn, pn_face);

        if strcmp(pn.p_idx, 'p27')
            %nodes = [nodes(1,:); nodes(3,:); nodes(end,:); nodes(end-1,:); nodes(2,:)];
            nodes = flip(sortrows(nodes, 2));
        elseif strcmp(pn.p_idx, 'p7small') || strcmp(pn.p_idx, 'p11left')
            nodes = sortrows(nodes, 2);
        elseif strcmp(pn.p_idx, 'p7left')
            p7small_exclude = poly_obj.p7small.G.nodes.coords(~poly_obj.p7small.top_mask, :);
            nodes = setdiff(nodes, p7small_exclude, 'rows');
        elseif any(strcmp(pn_face, 'west')) || any(strcmp(pn_face, 'top')) ...
                || strcmp(pn.p_idx, 'p29B') || strcmp(pn.p_idx, 'p19B')
            nodes = flip(nodes); % necessary to get in counter-clockwise order
        end
        all_nodes = [all_nodes; nodes];
        all_pts = [all_pts; pts];
        %facies = [facies; repmat(pn.facies, numel(nodes), 1)];
    end
    [~, unique_idx] = uniquetol(all_nodes, 'ByRows', true); % stable to give in same order as list of neighbors
    external_nodes = all_nodes(sort(unique_idx), :); 
    [~, unique_idx] = uniquetol(all_pts, 'ByRows', true); % stable to give in same order as list of neighbors
    external_pts = all_pts(sort(unique_idx), :); 

    poly_obj.pBFPEBI.bnodes = external_nodes;
    poly_obj.pBFPEBI.p = external_pts;
   
    internal_nodes = [];
    % 2. Assign internal boundary POINTS (no nodes defined here yet)
    for p=1:numel(poly_num)
        poly_n = poly_num(p);
        poly = Faults(all_polys, poly_n, 1, false); 
        p_idx_pebi = strcat('p', string(poly_n), 'PEBI');
        poly.p_idx = p_idx_pebi;                
    
       poly_obj.(p_idx_pebi) = poly;

       internal_nodes = [internal_nodes; poly.p];

    end
end