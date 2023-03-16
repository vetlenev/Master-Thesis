function [poly_obj, nodes_overlap, pts_overlap] = glueUpperPolygonMultiple(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                           nodes_overlap, pts_overlap, G_glob, varargin)  

    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,2));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);

    %% SMALL, LEFT POLYGON
    if poly_num == 7
        poly = PolygonGrid(all_polys, poly_num); 
        p_idx_small = strcat('p', string(poly_num), 'small');
        poly.p_idx = p_idx_small;
        poly.p = poly.p(21:25, :);
    
        poly.p_we = poly.p(3,:);
        poly.we = 1;
        poly.p_bt = poly.p([1:2, 4:5], :);
        poly.p = [poly.p_bt; poly.p_we];
        poly.top_side = poly.p_bt(3:4,:);
        poly.bottom_side = poly.p_bt(1:2,:);
    
        poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, fix(Nx_glob/100)); % no overlap with upper polygon, fix a very small horizontal dimension (to not get too tight cells)            
    
        [poly, bottom_idx] = correct_dx_poly_new(poly, G_glob, poly.bottom_mask, 'bottom');
        [poly, top_idx] = correct_dx_poly_new(poly, G_glob, poly.top_mask, 'top');
    
        poly = interpolateZ_remaining_new(poly, bottom_idx, poly.bottom_mask, 'bottom', 'linear');
        poly = interpolateZ_remaining_new(poly, top_idx, poly.top_mask, 'top', 'linear');
    
        poly = interpolateInternal(poly, poly.top_mask, poly.bottom_mask, []);    
        poly = interpolateSide(poly);
        poly = interpolateHorizontal(poly);
        %poly = interpolatePartlyHorizontal(poly, 0.2);
    
        Nx_small = poly.G.cartDims(1)+1;
        Nz_small = poly.G.cartDims(2)+1;
        poly_obj.(p_idx_small) = poly;
        % NB: Here the overlaps are to the east, not top!
        nodes_overlap.(p_idx_small) = poly.G.nodes.coords(poly.east_mask, :); % these nodes are used to glue left part of polygon 7 to the small part
        east_nodes_small = nodes_overlap.(p_idx_small);
        pts_overlap.(p_idx_small) = [poly.bottom_side(1,:); poly.top_side(end,:)];
    end
    %% PARENT POLYGON
    poly = PolygonGrid(all_polys, poly_num);
    p_idx = strcat('p', string(poly_num));
    
    p_idx_upper = {strcat('p', string(poly_num_upper), 'left'), ...
                        strcat('p', string(poly_num_upper), 'mid'), ...
                        strcat('p', string(poly_num_upper), 'right')};
    poly_upper = {poly_obj.(p_idx_upper{1}), ...
                   poly_obj.(p_idx_upper{2}), ...
                   poly_obj.(p_idx_upper{3})};                      
    
    if poly_num == 7
        % Remove small left part
        poly.p(22:24, :) = []; % new indexing: 25:end -> 22:end            
    end
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
    if poly_num == 7
        poly.bottom_side = poly.p(1:21, :);    
    elseif poly_num == 2
        poly.bottom_side = poly.p(1:15, :);
    end
    % Create (dummy) subgrid for full polygon
    poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap);           
    poly.bottom_mask = poly.G.nodes.coords(:,2) == min(poly.G.nodes.coords(:,2));
    % Create bottom nodes for parent polygon (this is queried later)
    [poly, closest_bottom] = correct_dx_poly_new(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
    poly = interpolateZ_remaining_new(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'linear');   

    %% Handle LEFT part
    p_idx_left = strcat(p_idx, 'left');
    poly_obj.(p_idx_left) = PolygonGrid(all_polys, poly_num);           
    poly_left = poly_obj.(p_idx_left); % dummy variable for readability    
    poly_left.p_idx = p_idx_left;

    poly_left.p = poly.p; 
    poly_left.top_side = poly.top_side;
    poly_left.bottom_side = poly.bottom_side;
    if poly_num == 7
        poly_left.p_we = poly.p(22, :);
        poly_left.we = 1;
    elseif poly_num == 2
        poly_left.p_we = poly.p(16:20, :);
        poly_left.we = [1,1,1,1,1]';
    end
   
    [nodes_overlap.(p_idx_left), pts_overlap.(p_idx_left)] = findOverlappingNodes(poly_left, poly_upper{1}, 'top'); % left part only overlaps with first upper polygon
    nodes_overlap_left = nodes_overlap.(p_idx_left);
    num_overlap_left = size(nodes_overlap_left, 1);
    poly_left = cartesianSubgrid(poly_left, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_left); 

    % Fix top nodes
    poly_left.G.nodes.coords(poly_left.top_mask, :) = nodes_overlap.(p_idx_left); % glue to overlapping nodes on top side       
    
    % Query associated part of bottom side from parent polygon  
    bottom_nodes = poly.G.nodes.coords(poly.bottom_mask, :);

    sep_point_left = nodes_overlap_left(nodes_overlap_left(:,1) == max(nodes_overlap_left(:,1)), :);
    [~, sep_left_idx] = min(abs(bottom_nodes(:,1) - sep_point_left(1)));
    bottom_left_end = bottom_nodes(sep_left_idx, :);

    poly_left.top_side = poly_left.top_side(poly_left.top_side(:,1) <= sep_point_left(1), :);
    poly_left.bottom_side = poly_left.bottom_side(poly_left.bottom_side(:,1) <= sep_point_left(1), :);
    poly_left.bottom_side = [poly_left.bottom_side; bottom_left_end]; % add endpoint to be used for interpolation on bottom side
    poly_left.p_added = [poly_left.p_added; bottom_left_end];
    % Set bottom nodes (by interpolation up to east end of left polygon)
    [poly_left, closest_bottom] = correct_dx_poly_new(poly_left, G_glob, poly_left.bottom_mask, 'bottom');
    poly_left = interpolateZ_remaining_new(poly_left, closest_bottom, ...
                                        poly_left.bottom_mask, 'bottom', 'linear');
    
    poly_left = interpolateInternal(poly_left, poly_left.top_mask, poly_left.bottom_mask, []);    

    if poly_num == 7
        % Fix west side points to conform with left, small polygon
        Nx_left = poly_left.G.cartDims(1)+1;
        Nz_left = poly_left.G.cartDims(2)+1;    
        poly_left.G.nodes.coords(1:Nx_left:Nx_left*Nz_small, :) = east_nodes_small;
    
        west_interp_bottom = east_nodes_small(east_nodes_small(:,2) == max(east_nodes_small(:,2)), :);
        poly_left_top = poly_left.G.nodes.coords(poly_left.top_mask, :);
        west_interp_top = poly_left_top(poly_left_top(:,1) == min(poly_left_top(:,1)), :);
    
        Nz_left_rem = Nz_left - Nz_small + 1; % point in vertical ABOVE small part
        xt = west_interp_top(1); zt = west_interp_top(2);
        xb = west_interp_bottom(1); zb = west_interp_bottom(2);
        dx = (xt - xb)/(Nz_left_rem - 1);
        dz = (zt - zb)/(Nz_left_rem - 1);
        x_interp = xb + cumsum(repmat(dx, Nz_left_rem, 1)) - dx;
        z_interp = zb + cumsum(repmat(dz, Nz_left_rem, 1)) - dz;
    
        poly_left.G.nodes.coords(1+Nx_left*(Nz_small-1):Nx_left:Nx_left*Nz_left, 1) = x_interp;
        poly_left.G.nodes.coords(1+Nx_left*(Nz_small-1):Nx_left:Nx_left*Nz_left, 2) = z_interp;
    elseif poly_num == 2
        poly_left = interpolateSide(poly_left);
        %poly_left = interpolateHorizontal(poly_left);
    end    
    poly_left = logicalIndicesUpperOverlap(poly_left, poly_upper{1});
    poly_obj.(p_idx_left) = poly_left;

    %% Handle MIDDLE part
    p_idx_mid = strcat(p_idx, 'mid');
    poly_obj.(p_idx_mid) = PolygonGrid(all_polys, poly_num);           
    poly_mid = poly_obj.(p_idx_mid); % dummy variable for readability    
    poly_mid.p_idx = p_idx_mid;

    poly_mid.p = poly.p; 
    poly_mid.top_side = poly.top_side;
    poly_mid.bottom_side = poly.bottom_side;
    poly_mid.p_we = [];
    poly_mid.we = [];

    [nodes_overlap.(p_idx_mid), pts_overlap.(p_idx_mid)] = findOverlappingNodes(poly_mid, poly_upper{2}, 'top'); % mid part only overlaps with second upper polygon
    nodes_overlap_mid = nodes_overlap.(p_idx_mid);
    num_overlap_mid = size(nodes_overlap_mid, 1);
    poly_mid = cartesianSubgrid(poly_mid, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_mid); 

    % Fix top nodes
    poly_mid.G.nodes.coords(poly_mid.top_mask, :) = nodes_overlap.(p_idx_mid); % glue to overlapping nodes on top side       
    
    % Query associated part of bottom side from parent polygon   
    sep_point_right = nodes_overlap_mid(nodes_overlap_mid(:,1) == max(nodes_overlap_mid(:,1)), :);

    [~, sep_right_idx] = min(abs(bottom_nodes(:,1) - sep_point_right(1)));
    bottom_right_start = bottom_nodes(sep_right_idx, :);
    poly_mid.top_side = poly_mid.top_side(poly_mid.top_side(:,1) >= sep_point_left(1) & ...
                                           poly_mid.top_side(:,1) <= sep_point_right(1), :);

    poly_mid.bottom_side = poly_mid.bottom_side(poly_mid.bottom_side(:,1) >= sep_point_left(1) & ...
                                                  poly_mid.bottom_side(:,1) <= sep_point_right(1), :);
    poly_mid.bottom_side = [bottom_left_end; poly_mid.bottom_side; bottom_right_start]; % add endpoint to be used for interpolation on bottom side
    poly_mid.p_added = [poly_mid.p_added; bottom_left_end; bottom_right_start];
    % Set bottom nodes (by interpolation up to pinch point)
    [poly_mid, closest_bottom] = correct_dx_poly_new(poly_mid, G_glob, poly_mid.bottom_mask, 'bottom');
    poly_mid = interpolateZ_remaining_new(poly_mid, closest_bottom, ...
                                        poly_mid.bottom_mask, 'bottom', 'linear');

    poly_mid = interpolateInternal(poly_mid, poly_mid.top_mask, poly_mid.bottom_mask, []); 

    %poly_mid = interpolateHorizontal(poly_mid);

    poly_mid = logicalIndicesUpperOverlap(poly_mid, poly_upper{2});
    poly_obj.(p_idx_mid) = poly_mid;

    %% Handle RIGHT part
    p_idx_right = strcat(p_idx, 'right');
    poly_obj.(p_idx_right) = PolygonGrid(all_polys, poly_num);           
    poly_right = poly_obj.(p_idx_right); % dummy variable for readability    
    poly_right.p_idx = p_idx_right;

    poly_right.p = poly.p; 
    poly_right.top_side = poly.top_side;
    poly_right.bottom_side = poly.bottom_side;
    poly_right.p_we = [];
    poly_right.we = [];
    
    [nodes_overlap.(p_idx_right), pts_overlap.(p_idx_right)] = findOverlappingNodes(poly_right, poly_upper{3}, 'top'); % right part only overlaps with second upper polygon
    nodes_overlap_right = nodes_overlap.(p_idx_right);
    num_overlap_right = size(nodes_overlap_right, 1);
    poly_right = cartesianSubgrid(poly_right, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_overlap_right); 

    poly_right.G.nodes.coords(poly_right.top_mask, :) = nodes_overlap.(p_idx_right); % glue to overlapping nodes on top side       
   
    poly_right.top_side = poly_right.top_side(poly_right.top_side(:,1) >= sep_point_right(1), :);

    poly_right.bottom_side = poly_right.bottom_side(poly_right.bottom_side(:,1) >= sep_point_right(1), :);
    poly_right.bottom_side = [bottom_right_start; poly_right.bottom_side]; % add endpoint to be used for interpolation on bottom side
    poly_right.p_added = [poly_right.p_added; bottom_right_start];
    
    [poly_right, closest_bottom] = correct_dx_poly_new(poly_right, G_glob, poly_right.bottom_mask, 'bottom');
    poly_right = interpolateZ_remaining_new(poly_right, closest_bottom, ...
                                        poly_right.bottom_mask, 'bottom', 'linear');

    poly_right = interpolateInternal(poly_right, poly_right.top_mask, poly_right.bottom_mask, []); 

    %poly_right = interpolateHorizontal(poly_right);

    poly_right = logicalIndicesUpperOverlap(poly_right, poly_upper{3});
    poly_obj.(p_idx_right) = poly_right;

end
