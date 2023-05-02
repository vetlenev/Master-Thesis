function [poly_obj, nodes_overlap, pts_overlap] = triangulateFault(all_polys, poly_obj, poly_num, poly_neighbors, ...
                                                                  nodes_overlap, pts_overlap, G_glob, varargin)  
   
% NB: poly_neighbors must be defined in counter-clockwise order around
% fault!
    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,1));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);
   
    if nargin > 7
        fac_scale = varargin{1};
    else
        fac_scale = 1;
    end

    poly = Faults(all_polys, poly_num, fac_scale, false); 
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

        is_west = any(strcmp(pn_face, 'west'));
        is_top = any(strcmp(pn_face, 'top'));
        is_east = any(strcmp(pn_face, 'east'));
        is_bottom = any(strcmp(pn_face, 'bottom'));

        if strcmp(pn.p_idx, 'p29C')
            nodes = [nodes(1:end-4,:); flip(nodes(end-3:end,:))];
        elseif strcmp(pn.p_idx, 'p19C')
            nodes = [nodes(1:end-6,:); nodes(end-4,:); flip(nodes(end-3:end,:)); nodes(end-5,:)];
        elseif is_top || (is_west && poly_num == 25) ...
                      || (is_east && poly_num == 31) ...
                      || strcmp(pn.p_idx, 'p29B') ...
                      || strcmp(pn.p_idx, 'p19B')
            nodes = flip(nodes); % necessary to get in counter-clockwise order               
        end
        all_nodes = [all_nodes; nodes];
        all_pts = [all_pts; pts];
        facies = [facies; repmat(pn.facies, size(nodes,1), 1)];
    end

    % interpolate east side of p12 fault
    if poly_num == 12
        p = poly.p;
        p_east = p(p(:,1) == max(p(:,1)), :);
        z_min = min(p_east(:,2));
        z_max = max(p_east(:,2));
        Lz = z_max - z_min;
        fac = Lz/Lz_glob;
        Nz = ceil(poly.fac_scale*fac*Nz_glob);
        dz = Lz/Nz;
        z = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;
        x = repmat(max(p(:,1)), Nz+1, 1);
        east_nodes = [x,z];
        all_nodes = [all_nodes; east_nodes];
        facies = [facies; zeros(size(east_nodes,1), 1)]; % no neighbor at east side -> set facie to zero
    end

    [~, unique_idx] = uniquetol(all_nodes, 'ByRows', true); % stable to give in same order as list of neighbors
    poly.bnodes = all_nodes(sort(unique_idx), :);
    poly.bfacies = facies(sort(unique_idx), :); % Do we need this ??
 
    poly_obj.(p_idx_pebi) = poly;
    nodes_overlap.(p_idx_pebi) = poly.bnodes;
end
