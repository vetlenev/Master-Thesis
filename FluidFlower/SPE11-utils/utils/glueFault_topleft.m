function [poly_obj, nodes_overlap, pts_overlap] = glueFault_topleft(all_polys, poly_obj, poly_num, poly_neighbors, ...
                                                                  nodes_overlap, pts_overlap, G_glob, varargin)  
   
% NB: poly_neighbors must be defined in counter-clockwise order around
% fault!
    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,1));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);
   
    poly = Faults(all_polys, poly_num, 1, false); 
    p_idx_pebi = strcat('p', string(poly_num), 'PEBI');
    poly.p_idx = p_idx_pebi;   

    % Find all overlapping nodes with neighbors -> this gives boundary
    % nodes of fault
    all_nodes = [];
    all_pts = [];
    facies = [];
    for i=1:size(poly_neighbors,1)
        pn = poly_neighbors{i,1}; % neighboring polygon instance        
        pn_face = poly_neighbors{i,2}; % intersecting face of neighbor
        [nodes, pts] = findOverlappingNodes(poly, pn, pn_face);
        if strcmp(pn_face, 'east') || strcmp(pn_face, 'top')
            nodes = flip(nodes); % necessary to get in counter-clockwise order
            %pts = flip(pts);
        end
        all_nodes = [all_nodes; nodes];
        all_pts = [all_pts; pts];
        facies = [facies; repmat(pn.facies, numel(nodes), 1)];
    end
    [~, unique_idx] = uniquetol(all_nodes, 'ByRows', true); % stable to give in same order as list of neighbors
    poly.bnodes = all_nodes(sort(unique_idx), :);
    poly.bfacies = facies(sort(unique_idx), :); % Do we need this ??
%     [~, unique_idx] = uniquetol(all_pts, 'ByRows', true);   
%     poly.bnodes = all_pts(sort(unique_idx), :);

   poly_obj.(p_idx_pebi) = poly;
   nodes_overlap(p_idx_pebi) = poly.bnodes;
end
