function [poly_obj, nodes_overlap, pts_overlap] = glueFault_topright(all_polys, poly_obj, poly_num, poly_neighbors, ...
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
        if strcmp(pn.p_idx, 'p29C')
            nodes = [nodes(1:end-4,:); flip(nodes(end-3:end,:))];
        elseif strcmp(pn.p_idx, 'p19C')
            nodes = [nodes(1:end-6,:); nodes(end-4,:); flip(nodes(end-3:end,:)); nodes(end-5,:)];
        elseif any(strcmp(pn_face, 'west')) || any(strcmp(pn_face, 'top')) ...
                || strcmp(pn.p_idx, 'p29B') || strcmp(pn.p_idx, 'p19B')
            nodes = flip(nodes); % necessary to get in counter-clockwise order
        end
        all_nodes = [all_nodes; nodes];
        all_pts = [all_pts; pts];
        facies = [facies; repmat(pn.facies, numel(nodes), 1)];
    end
    [~, unique_idx] = uniquetol(all_nodes, 'ByRows', true); % stable to give in same order as list of neighbors
    poly.bnodes = all_nodes(sort(unique_idx), :);     


   %% Make PEBI subgrid
    poly = Faults.globalSubgrid(poly, Lx, Lz, Nx, Nz, []);
    
    %% In polygon
    poly = Faults.cellsInsideFault(poly);
    
    %% Plot grid
    figure()
    pG = poly.G;
    bneighbors = unique(pG.faces.neighbors(pG.faces.tag, :));
    plotGrid(pG)
    bcells = zeros(pG.cells.num, 1);
    bcells(bneighbors) = 1;
    plotCellData(pG, bcells)
    
    figure()
    pGF = poly.G_fault;
    plotGrid(pGF)
    
    %% Set facies
    [poly, bneighbors] = Faults.assignFacies25(poly, poly_neighbors);
    
    %% Find discontinuous face transitions
    pGF = poly.G_fault;
    tip_faces = find(pGF.faces.centroids(:,2) < 0.673);
    tip_faces = setdiff(tip_faces, boundaryFaces(pGF));
    
    edge_case_faces = find(pGF.faces.centroids(:,1) > 1.433 & ...
                            pGF.faces.centroids(:,1) < 1.44 & ...
                            pGF.faces.centroids(:,2) > 0.988 & ...
                            pGF.faces.centroids(:,2) < 0.9897);  % somehow removeRedundantFaces does not work on this face ..., add it manually
    edge_case_faces = [edge_case_faces; find(pGF.faces.centroids(:,1) > 1.245 & ...
                                            pGF.faces.centroids(:,1) < 1.25 & ...
                                            pGF.faces.centroids(:,2) > 0.904 & ...
                                            pGF.faces.centroids(:,2) < 0.906)];  % somehow removeRedundantFaces does not work on this face ..., add it manually
    edge_case_faces = [edge_case_faces; find(pGF.faces.centroids(:,1) > 1.22 & ...
                                            pGF.faces.centroids(:,1) < 1.225 & ...
                                            pGF.faces.centroids(:,2) > 0.894 & ...
                                            pGF.faces.centroids(:,2) < 0.896)];  % somehow removeRedundantFaces does not work on this face ..., add it manually
    
    
    tip_faces = [tip_faces; edge_case_faces];
    faces2keep = find(pGF.faces.centroids(:,1) > 1.2085 & ...
                        pGF.faces.centroids(:,1) < 1.2211 & ...
                        pGF.faces.centroids(:,2) > 0.895 & ...
                        pGF.faces.centroids(:,2) < 0.8975);
    
    poly = Faults.removeRedundantFaces(poly, bneighbors, tip_faces, faces2keep);
    
    %% Fix faces unmatched with global grid
    point_set = {[1,2], [100,109], [109,115]};
    poly = Faults.fixUnmatchedFaces25(poly, point_set);
    pGF = poly.G_fault; % update
    
    %% Change coordinate of degenerate face
    bf = boundaryFaces(pGF);
    bf_mask = zeros(pGF.faces.num, 1);
    bf_mask(bf) = 1;
    
    old_face = find(pGF.faces.centroids(:,2) == min(pGF.faces.centroids(~bf_mask,2)));
    true_coord_left = poly.bnodes(108,:);
    true_coord_right = poly.bnodes(110,:);
    
    poly = Faults.changeNodesForFaces(poly, old_face, true_coord_left, true_coord_right);

    poly_obj.(p_idx_pebi) = poly;
    nodes_overlap(p_idx_pebi) = poly.bnodes;
end
