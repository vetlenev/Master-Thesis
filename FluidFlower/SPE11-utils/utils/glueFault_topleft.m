function [all_nodes, nodes_overlap, pts_overlap] = glueFault_topleft(all_polys, poly_obj, poly_num, poly_neighbors, ...
                                                                  nodes_overlap, pts_overlap, G_glob, varargin)  
   Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,1));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);
   
    poly = PolygonGrid(all_polys, poly_num); 
    p_idx_pebi = strcat('p', string(poly_num), 'PEBI');
    poly.p_idx = p_idx_pebi;   

    % Find all overlapping nodes with neighbors -> this gives boundary
    % nodes of fault
    all_nodes = [];
    for i=1:size(poly_neighbors,1)
        pn = poly_neighbors{i,1}; % neighboring polygon instance
        pn_face = poly_neighbors{i,2}; % intersecting face of neighbor
        [nodes, pts] = findOverlappingNodes(poly, pn, pn_face);
        all_nodes = [all_nodes; nodes];
    end
    all_nodes = unique(all_nodes, 'rows');
       
end
