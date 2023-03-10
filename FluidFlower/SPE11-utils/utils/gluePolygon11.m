function [poly_obj, nodes_overlap, pts_overlap] = gluePolygon11(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                  nodes_overlap, pts_overlap, G_glob, varargin)  
    p_idx = strcat('p', string(poly_num));
    
    p_idx_upper = {strcat('p', string(poly_num_upper(1))), ...
                        strcat('p', string(poly_num_upper(2))), ...
                        strcat('p', string(poly_num_upper(3)), 'A'), ...
                        strcat('p', string(poly_num_upper(3)), 'C')};
    poly_upper = {poly_obj.(p_idx_upper{1}), ...
                   poly_obj.(p_idx_upper{2}), ...
                   poly_obj.(p_idx_upper{3}), ...
                   poly_obj.(p_idx_upper{4})};

    poly = PolygonGrid(all_polys, poly_num);               

    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,2));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);   

    [nodes_overlap_left, pts_overlap_left] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');    
        
    top_nodes = nodes_overlap_left;
    top_pts = pts_overlap_left;
    num_x_overlap = size(top_nodes, 1);
        
    poly.bottom_side = poly.p(3:17, :);
    poly.top_side = top_pts;
    % Create (dummy) subgrid for full polygon
    poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap);           
    % CONTINUE HERE:
    pinch_left = poly.p(30, :);
    pinch_right = poly.p(28, :);

    % Handle LEFT part of polygon 11
    p_idx_left = strcat(p_idx, 'left');
    poly_obj.(p_idx_left) = PolygonGrid(all_polys, poly_num);           
    poly_left = poly_obj.(p_idx_left); % dummy variable for readability    
    poly_left.p_idx = p_idx_left;

    poly_left.p = poly.p([1:8, 30:end], :); 
    poly_left.top_side = poly.p(33:end, :);
    poly_left.bottom_side = poly.p(3:8, :);
    poly_left.p_we = poly.p([1:2, 30:32], :);
    poly_left.we = [1,1,0,0,0]';

    [nodes_overlap.(p_idx_left), pts_overlap.(p_idx_left)] = findOverlappingNodes(poly_left, poly_upper{1}, 'top'); % left part only overlaps with first upper polygon
    num_overlap_left = size(nodes_overlap.(p_idx_left), 1);
    poly_left = cartesianSubgrid(poly_left, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_left); 
       
    % Fix top nodes
    poly_left.G.nodes.coords(poly_left.top_mask, :) = nodes_overlap.(p_idx_left); % glue to overlapping nodes on top side       
    % Create bottom nodes for parent polygon (this is queried later)
    [poly, closest_bottom] = correct_dx_poly(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
    poly = interpolateZ_remaining(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'spline');
    % Query associated part of bottom side from parent polygon  
    bottom_nodes = poly.G.nodes.coords(poly.bottom_mask, :);      
    [~, pinch_left_idx] = min(abs(bottom_nodes(:,1) - pinch_left(1)));
    bottom_left_end = bottom_nodes(pinch_left_idx, :);
    poly_left.bottom_side = [poly_left.bottom_side; bottom_left_end];
    poly_left.p_added = [poly_left.p_added; bottom_left_end];
    % Set bottom nodes (by interpolation up to pinch point)
    [poly_left, closest_bottom] = correct_dx_poly(poly_left, G_glob, poly_left.bottom_mask, 'bottom');
    poly_left = interpolateZ_remaining(poly_left, closest_bottom, ...
                                        poly_left.bottom_mask, 'bottom', 'spline');

    % Interpolate internal nodes (not accounting for side-points yet)
    poly_left = interpolateInternal(poly_left, poly_left.top_mask, poly_left.bottom_mask, []);    
    % Now we can fix side nodes (assigned to closest internal node)
    poly_left = interpolateSide(poly_left);
    % Reinterpolate horizontally after shifting boundary nodes
    poly_left = interpolateHorizontal(poly_left);

    poly_left = logicalIndicesUpperOverlap(poly_left, poly_upper{1});
    poly_left.p = [poly_left.p; poly_left.p_added];
    poly_obj.(p_idx_left) = poly_left;

    % Handle MIDDLE part
    p_idx_mid = strcat(p_idx, 'mid');
    poly_obj.(p_idx_mid) = PolygonGrid(all_polys, poly_num);           
    poly_mid = poly_obj.(p_idx_mid); % dummy variable for readability    
    poly_mid.p_idx = p_idx_mid;

    poly_mid.p = poly.p([9:10, 28:30], :); 
    poly_mid.top_side = poly.p(28:30, :);
    poly_mid.bottom_side = poly.p(9:10, :);
    poly_mid.p_we = [];
    poly_mid.we = [];

    G_left = poly_left.G;
    Nx_left = G_left.cartDims(1)+1;
    Nz_left = G_left.cartDims(2)+1;
    G_left_east_nodes = G_left.nodes.coords(Nx_left:Nx_left:Nx_left*Nz_left, 2);
    Nz_mid = nnz(G_left_east_nodes <= pinch_left(2)) - 1; % -1 since one less cell than nodes

    [nodes_overlap.(p_idx_mid), pts_overlap.(p_idx_mid)] = findOverlappingNodes(poly_mid, poly_upper{2}, 'top'); % middle part only overlaps with second upper polygon
    num_overlap_mid = size(nodes_overlap.(p_idx_mid), 1);
    poly_mid = cartesianSubgrid(poly_mid, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_mid, Nz_mid); 
      
    % Fix top nodes
    poly_mid.G.nodes.coords(poly_mid.top_mask, :) = nodes_overlap.(p_idx_mid); % glue to overlapping nodes on top side       
   
    [~, pinch_right_idx] = min(abs(bottom_nodes(:,1) - pinch_right(1)));
    bottom_mid_start = bottom_left_end;
    bottom_mid_end = bottom_nodes(pinch_right_idx, :);
    poly_mid.bottom_side = [bottom_mid_start; poly_mid.bottom_side; bottom_mid_end];
    poly_mid.p_added = [poly_mid.p_added; bottom_mid_start; bottom_mid_end];

    [poly_mid, closest_bottom] = correct_dx_poly(poly_mid, G_glob, poly_mid.bottom_mask, 'bottom');
    poly_mid = interpolateZ_remaining(poly_mid, closest_bottom, ...
                                        poly_mid.bottom_mask, 'bottom', 'spline');

    poly_mid = interpolateInternal(poly_mid, poly_mid.top_mask, poly_mid.bottom_mask, []);    
    poly_mid = interpolateSide(poly_mid);
    poly_mid = interpolateHorizontal(poly_mid);

    poly_mid = logicalIndicesUpperOverlap(poly_mid, poly_upper{2});
    poly_mid.p = [poly_mid.p; poly_mid.p_added];
    poly_obj.(p_idx_mid) = poly_mid;

    % Handle RIGHT part
    p_idx_right = strcat(p_idx, 'right');
    poly_obj.(p_idx_right) = PolygonGrid(all_polys, poly_num);           
    poly_right = poly_obj.(p_idx_right); % dummy variable for readability    
    poly_right.p_idx = p_idx_right;

    poly_right.p = poly.p(11:28, :); 
    poly_right.top_side = poly.p(18:26, :);
    poly_right.bottom_side = poly.p(11:17, :);
    poly_right.p_we = poly.p(27:28, :);
    poly_right.we = [1,1]';
    
    z_left = poly_left.top_side(1,2) - poly_left.bottom_side(end,2);
    z_right = poly_right.top_side(end,2) - poly_right.bottom_side(1,2);
    Nz_right = ceil(Nz_left*(z_right/z_left));

    [nodes_overlap.(p_idx_right), pts_overlap.(p_idx_right)] = PinchOuts.findOverlappingNodesMultiple(poly_right, {poly_upper{3},poly_upper{4}}, 'top'); % right part overlaps with third AND fourth upper polygon
    num_overlap_right = size(nodes_overlap.(p_idx_right), 1);
    poly_right = cartesianSubgrid(poly_right, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_right); 
    
    % Fix top nodes
    poly_right.G.nodes.coords(poly_right.top_mask, :) = nodes_overlap.(p_idx_right); % glue to overlapping nodes on top side       
   
    bottom_right_start = bottom_mid_end;
    poly_right.bottom_side = [bottom_right_start; poly_right.bottom_side];
    poly_right.p_added = [poly_right.p_added; bottom_right_start];

    [poly_right, closest_bottom] = correct_dx_poly(poly_right, G_glob, poly_right.bottom_mask, 'bottom');
    poly_right = interpolateZ_remaining(poly_right, closest_bottom, ...
                                        poly_right.bottom_mask, 'bottom', 'spline');

    poly_right = interpolateInternal(poly_right, poly_right.top_mask, poly_right.bottom_mask, []);    
    poly_right = interpolateSide(poly_right);
    poly_right = interpolateHorizontal(poly_right);

    poly_right = logicalIndicesUpperOverlap(poly_right, {poly_upper{3},poly_upper{4}});
    poly_right.p = [poly_right.p; poly_right.p_added];
    poly_obj.(p_idx_right) = poly_right;
end
