function [poly_obj, nodes_overlap, pts_overlap, ...
            top_nodes, bottom_nodes, west_nodes] = gluePolygon1New(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                  nodes_overlap, pts_overlap, G_glob, node_density, varargin)  
   Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,2));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);

    p_idx_upper = {strcat('p', string(poly_num_upper(1))), ...
                        strcat('p', string(poly_num_upper(2)), 'left'), ...
                        strcat('p', string(poly_num_upper(2)), 'mid'), ...
                        strcat('p', string(poly_num_upper(2)), 'right')};
    poly_upper = {poly_obj.(p_idx_upper{1}), ...
                   poly_obj.(p_idx_upper{2}), ...
                   poly_obj.(p_idx_upper{3}), ...
                   poly_obj.(p_idx_upper{4})}; 

    %% TRIANGULAR END-PIECE   
    poly_tri = PolygonGrid(all_polys, poly_num); 
    p_idx_small = strcat('p', string(poly_num), 'small');
    poly_tri.p_idx = p_idx_small;
    p_triangular_left = [poly_tri.p(20,1), 0]; % bottom-left of small triangle
    poly_tri.p(end+1,:) = p_triangular_left; % add point projected to bottom boundary to make triangular shape
    poly_tri.p = poly_tri.p([1,20,end], :);
    poly_tri.p_added = [poly_tri.p_added; poly_tri.p(end,:)];

    poly_tri.p_we = [];
    poly_tri.we = [];
    poly_tri.p_bt = poly_tri.p;
    poly_tri.top_side = poly_tri.p_bt([1,2],:);
    poly_tri.bottom_side = poly_tri.p_bt([1,3],:);    
   

    %% PARENT POLYGON
    poly = PolygonGrid(all_polys, poly_num);
    p_idx = strcat('p', string(poly_num));
    poly.p(end+1,:) = [poly.p(20,1), 0]; % add point projected to bottom boundary to make triangular shape                            
   
    % Remove right endpoint
    poly.p(1, :) = []; % new indexing: 1:20 -> 1:19             
    for i=1:numel(poly_upper)
        if ~isempty(poly_upper{i}.p_added)
            n_added = size(poly_upper{i}.p_added, 1);
            poly.p(end+1:end+n_added, :) = poly_upper{i}.p_added;
        end
    end

    [nodes_overlap_all, pts_overlap_all] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');    
        
    top_nodes = nodes_overlap_all;    
    num_x_overlap = size(top_nodes, 1);
        
    poly.top_side = pts_overlap_all; % poly points, not nodes   
    poly.bottom_side = poly.p([1,end], :);

    % Create (dummy) subgrid for full polygon
    poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap);
    % Create bottom nodes for parent polygon (this is queried later)
    [poly, closest_bottom] = correct_dx_poly_new(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
    poly = interpolateZ_remaining_new(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'linear');   

    %% Handle LEFT part
    p_idx_left = strcat(p_idx, 'left');
    poly_left = PolygonGrid(all_polys, poly_num);    
    poly_left.p_idx = p_idx_left;

    poly_left.p = poly.p; 
    poly_left.top_side = poly.top_side;
    poly_left.bottom_side = poly.bottom_side;
    poly_left.p_we = [];
    poly_left.we = [];
   
    [nodes_overlap.(p_idx_left), pts_overlap.(p_idx_left)] = findOverlappingNodes(poly_left, poly_upper{1}, 'top'); % left part only overlaps with first upper polygon
    nodes_overlap_left = nodes_overlap.(p_idx_left);
    num_overlap_left = size(nodes_overlap_left, 1);

    % Query associated part of bottom side from parent polygon
    sep_point_left = nodes_overlap_left(nodes_overlap_left(:,1) == max(nodes_overlap_left(:,1)), :);
    bottom_left_end = [sep_point_left(:,1), 0];
    poly_left.top_side = poly_left.top_side(poly_left.top_side(:,1) <= sep_point_left(1), :);
    poly_left.bottom_side = poly_left.bottom_side(poly_left.bottom_side(:,1) <= sep_point_left(1), :);
    poly_left.bottom_side = [poly_left.bottom_side; bottom_left_end]; % add endpoint to be used for interpolation on bottom side
    poly_left.p_added = [poly_left.p_added; bottom_left_end];

    Lz_left = max(poly_left.top_side(:,2));
    Nz_left = ceil(Lz_left/Lz_glob * Nz_glob);

    poly_left = cartesianSubgrid(poly_left, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_left, Nz_left); 
  
    poly_left.G.nodes.coords(poly_left.top_mask, :) = nodes_overlap_left; % glue to overlapping nodes on top side       
       
    % Set bottom nodes (by interpolation up to east end of left polygon)
    xmin = min(poly_left.bottom_side(:,1));
    xmax = max(poly_left.bottom_side(:,1));
    Nx_left = poly_left.G.cartDims(1);
    Nz_left = poly_left.G.cartDims(2);
    dx = (xmax - xmin)/Nx_left;
    poly_left.G.nodes.coords(poly_left.bottom_mask, 1) = xmin + cumsum(repmat(dx, Nx_left+1, 1)) - dx;
    % z-value at bottom is zero by default -> no need for interpolateZ_rem
  
    poly_left = interpolateInternal(poly_left, poly_left.top_mask, poly_left.bottom_mask, []);    
    
    poly_left = logicalIndicesUpperOverlap(poly_left, poly_upper{1});
    % UNCOMMENT THIS IF INCLUDING LEFT PART:
    poly_obj.(p_idx_left) = poly_left;    

    %% Handle MIDDLE part
    p_idx_mid = strcat(p_idx, 'mid');
    poly_obj.(p_idx_mid) = PolygonGrid(all_polys, poly_num);           
    poly_mid = poly_obj.(p_idx_mid); % dummy variable for readability    
    poly_mid.p_idx = p_idx_mid;

    poly_mid.p = [poly.p(6:7, :); poly_left.p_added]; 
    poly_mid.top_side = poly.p(6:7, :);
    poly_mid.bottom_side = poly_left.p_added;
    poly_mid.p_we = [];
    poly_mid.we = [];
    
    poly_mid = cartesianSubgrid(poly_mid, Lx_glob, Lz_glob, Nx_glob, Nz_glob, [], Nz_left); 

    % Fix top nodes by interpolation
    xmin = min(poly_mid.top_side(:,1));
    zmin = min(poly_mid.top_side(:,2));
    xmax = max(poly_mid.top_side(:,1));
    zmax = max(poly_mid.top_side(:,2));
    Nx_mid = poly_mid.G.cartDims(1);
    Nz_mid = poly_mid.G.cartDims(2);
    dx = (xmax - xmin)/Nx_mid;
    dz = (zmax - zmin)/Nx_mid;
    x_interp = xmin + cumsum(repmat(dx, Nx_mid+1, 1)) - dx;
    z_interp = zmin + cumsum(repmat(dz, Nx_mid+1, 1)) - dz;
    poly_mid.G.nodes.coords(poly_mid.top_mask, 1) = x_interp; 
    poly_mid.G.nodes.coords(poly_mid.top_mask, 2) = z_interp; 

    % Add point to represent bottom right corner
    sep_point_right = poly_mid.p(poly_mid.p(:,1) == max(poly_mid.p(:,1)));
    bottom_right_start = [sep_point_right(:,1), 0];
    poly_mid.bottom_side = [poly_mid.bottom_side; bottom_right_start]; % add endpoint to be used for interpolation on bottom side
    poly_mid.p_added = [poly_mid.p_added; bottom_right_start];
  
    % Set bottom nodes (by interpolation up to pinch point)
    xmin = min(poly_mid.bottom_side(:,1));
    xmax = max(poly_mid.bottom_side(:,1));
    dx = (xmax - xmin)/Nx_mid;
    poly_mid.G.nodes.coords(poly_mid.bottom_mask, 1) = xmin + cumsum(repmat(dx, Nx_mid+1, 1)) - dx;

    poly_mid = interpolateInternal(poly_mid, poly_mid.top_mask, poly_mid.bottom_mask, []); 

    %poly_mid = logicalIndicesWestOverlap(poly_mid, poly_left);
    % No logical indices for this subgrid because it is right below fault!
    poly_mid.G.i = nan(poly_mid.G.cells.num,1);
    poly_mid.G.j = nan(poly_mid.G.cells.num,1);

    poly_obj.(p_idx_mid) = poly_mid;       

    %% Handle RIGHT part
    % Part A
    [poly_obj, pts_overlap, nodes_overlap] = constructSubgridRight(poly_obj, all_polys, poly_num, p_idx, poly, ...
                                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...   
                                                                    'A', poly_upper{2}, sep_point_right(1), ...
                                                                    pts_overlap, nodes_overlap);


    % Part B
    [poly_obj, pts_overlap, nodes_overlap] = constructSubgridRight(poly_obj, all_polys, poly_num, p_idx, poly, ...
                                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                                    'B', poly_upper{3}, max(poly_obj.p1rightA.top_side(:,1)), ...
                                                                    pts_overlap, nodes_overlap);

    % Part C
    [poly_obj, pts_overlap, nodes_overlap] = constructSubgridRight(poly_obj, all_polys, poly_num, p_idx, poly, ...
                                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                                    'C', poly_upper{4}, max(poly_obj.p1rightB.top_side(:,1)), ...
                                                                    pts_overlap, nodes_overlap);

    poly_rightA = poly_obj.p1rightA;
    poly_rightB = poly_obj.p1rightB;
    poly_rightC = poly_obj.p1rightC;

    %% Finish triangular part
    [nodes_overlap_all, pts_overlap_all] = findOverlappingNodes(poly_tri, poly_upper{4}, 'top'); 
    top_nodes = nodes_overlap_all;
    num_x_overlap = size(top_nodes,1);
    Nz_right = poly_rightC.G.cartDims(2);

    poly_tri = cartesianSubgrid(poly_tri, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap, Nz_right);

    poly_tri.G.nodes.coords(poly_tri.top_mask,:) = top_nodes;    

    west_nodes = poly_rightC.G.nodes.coords(poly_rightC.east_mask, :);
    poly_tri.G.nodes.coords(poly_tri.west_mask,:) = west_nodes;
    east_nodes = repmat(top_nodes(top_nodes(:,1) == max(top_nodes(:,1)), :), Nz_right+1, 1);
    poly_tri.G.nodes.coords(poly_tri.east_mask,:) = east_nodes;

    xmax = max(top_nodes(:,1));
    xmin = min(top_nodes(:,1));
    Nx_tri = poly_tri.G.cartDims(1);
    dx = (xmax - xmin)/Nx_tri;
    x = xmin + cumsum(repmat(dx, Nx_tri+1, 1)) - dx;
    z = zeros(numel(x), 1);

    bottom_nodes = [x,z];   
    poly_tri.G.nodes.coords(poly_tri.bottom_mask,:) = bottom_nodes;

    poly_tri = interpolateInternal(poly_tri, poly_tri.top_mask, poly_tri.bottom_mask, []); 

    poly_tri = logicalIndicesUpperOverlap(poly_tri, poly_upper{4}); 

    tri_bnodes = [flip(top_nodes); flip(west_nodes); bottom_nodes];       
    [~, unique_idx] = uniquetol(tri_bnodes, 'ByRows', true); % stable to give in same order as list of neighbors
    tri_bnodes = tri_bnodes(sort(unique_idx), :);
    poly_tri.bnodes = tri_bnodes;
  
    poly_obj.(p_idx_small) = poly_tri;
    pts_overlap.(p_idx_small) = pts_overlap_all;
    nodes_overlap.(p_idx_small) = top_nodes;
end

function [poly_obj, pts_overlap, nodes_overlap] = constructSubgridRight(poly_obj, all_polys, poly_num, p_idx, poly_parent, ...
                                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                                            right_part, poly_upper_neighbor, xstart_left, ...
                                                                            pts_overlap, nodes_overlap)
    p_idx_right = strcat(p_idx, 'right', right_part);
    poly_obj.(p_idx_right) = PolygonGrid(all_polys, poly_num);           
    poly_right = poly_obj.(p_idx_right); % dummy variable for readability    
    poly_right.p_idx = p_idx_right;

    poly_right.p = poly_parent.p; 
    poly_right.top_side = poly_parent.top_side;
    poly_right.bottom_side = poly_parent.bottom_side;
    poly_right.p_we = [];
    poly_right.we = [];
    
    [nodes_overlap.(p_idx_right), pts_overlap.(p_idx_right)] = findOverlappingNodes(poly_right, poly_upper_neighbor, 'top'); % part A only overlaps with second upper polygon
    nodes_overlap_right = nodes_overlap.(p_idx_right);
    num_overlap_right = size(nodes_overlap_right, 1);
    Nz_mid = poly_obj.p1mid.G.cartDims(2);

    poly_right = cartesianSubgrid(poly_right, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_right, Nz_mid); 

    poly_right.G.nodes.coords(poly_right.top_mask, :) = nodes_overlap_right; % glue to overlapping nodes on top side       
   
    poly_right.top_side = poly_right.top_side(poly_right.top_side(:,1) >= xstart_left & ...
                                                poly_right.top_side(:,1) <= max(nodes_overlap_right(:,1)), :);

    poly_right.bottom_side = poly_right.bottom_side(poly_right.bottom_side(:,1) >= xstart_left & ...
                                                     poly_right.bottom_side(:,1) <= max(nodes_overlap_right(:,1)), :);
    %poly_right.bottom_side = [bottom_right_start; poly_right.bottom_side]; % add endpoint to be used for interpolation on bottom side
    %poly_right.p_added = [poly_right.p_added; bottom_right_start];
    
    poly_right.G.nodes.coords(poly_right.bottom_mask, 1) = poly_right.G.nodes.coords(poly_right.top_mask, 1);
    poly_right.G.nodes.coords(poly_right.bottom_mask, 2) = 0;

    poly_right = interpolateInternal(poly_right, poly_right.top_mask, poly_right.bottom_mask, []); 

    poly_right = logicalIndicesUpperOverlap(poly_right, poly_upper_neighbor);  
    poly_obj.(p_idx_right) = poly_right;
end
