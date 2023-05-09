function [poly_obj, nodes_overlap, pts_overlap] = gluePinchOuts(all_polys, poly_obj, poly_num, poly_num_upper, poly_num_pinch, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)   
    % Create a structured, Cartesian grid for a layer with a pinch-out,
    % Layer is separated into three parts (A, B and C), where A conforms to
    % the tip of the pinch-out, B conforms to the top surface of the
    % pinch-out, and C conforms to the bottom surface of the pinch-out.
    % Resulting grids are compatible for gluing.
    % 
    % INPUTS:
    %   all_polys: polygonal points for each polygon of FluidFlower
    %   poly_obj: struct with properties of discretized polygon
    %   poly_num: index of polygon to discretize
    %   poly_num_upper: index of the polygon's upper neighbor
    %   poly_num_pinch: index of polygon representing the pinch-out
    %   nodes_overlap: struct of overlapping nodes used for gluing grids
    %   pts_overlap: struct of overlapping original polygonal points
    %   G_glob: global, virtual background grid for FluidFlower
    %
    % RETURNS:
    %   poly_obj: updated poly_obj struct
    %   nodes_overlap: updated nodes_overlap struct
    %   pts_overlap: updated pts_overlap struct
    %
    p_idx = strcat('p', string(poly_num));   

    p_name = fieldnames(poly_obj); 
    is_pinch = 0;
    for i=1:numel(p_name)
        is_pinch = is_pinch + ~isempty(regexp(p_name{i}, strcat('p', string(poly_num_upper), '[ABC]'), 'match'));
    end

    if is_pinch
        p_idx_upper = {strcat('p', string(poly_num_upper), 'A'), ...
                        strcat('p', string(poly_num_upper), 'C')}; % subgrid A and C together comprise the bottom side of upper polygon
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2})};
    else
        p_idx_upper = strcat('p', string(poly_num_upper));
        poly_upper = poly_obj.(p_idx_upper);
    end
    p_idx_pinch = strcat('p', string(poly_num_pinch));

    poly = PolygonGrid(all_polys, poly_num); % poly_obj.(p_idx) = ...
    poly_pinch = PolygonGrid(all_polys, poly_num_pinch);
    %poly = poly_obj.(p_idx); % dummy variable for readability    
    p = poly.p;
    pp = poly_pinch.p;

    % Create separation points at end of pinch-out   
    switch poly_num
        case 5
            pin_idx = 13;
        case 29
            pin_idx = 12;
        case 19
            pin_idx = 15;
    end
    pinch = p(pin_idx, :);
    z_pinch = p(pin_idx, 2);
    x_pinch = p(pin_idx, 1);

    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,2));
    N = G_glob.cartDims;
    Nx_glob = N(1); Nz_glob = N(2);          

    if numel(poly_upper) > 1 % Upper neighbor contains a pinch-out
        [nodes_overlap_pX, pts_overlap_pX] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');
    else
        [nodes_overlap_pX, pts_overlap_pX] = findOverlappingNodes(poly, poly_upper, 'top');
    end
    top_nodes = nodes_overlap_pX;
    top_pts = pts_overlap_pX;
    num_x_overlap = size(top_nodes, 1);
        
    % Find corresponding node in bottom side     
    switch poly_num
        case 5
            bottom_pts = [p(26:end,:); p(1:8,:)];
        case 29
            bottom_pts = p(2:6, :);
        case 19
            bottom_pts = p(8:11, :);
    end
    poly.bottom_side = bottom_pts;
    poly.top_side = top_pts;
    % Create (dummy) subgrid for full polygon
    poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap);           

    [bottom_nodes, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(poly, poly.bottom_side, G_glob);
    poly.bottom_side_new = poly.bottom_side;
    poly.top_side_new = poly.top_side;
    poly.bottom_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(poly.bottom_side_new, bottom_nodes, closest_mask, 'linear'); % changed first argument from poly.bottom_side to poly.bottom_side_new
    bottom_nodes(rem_idx,2) = z_new;

    % Find top/bottom-node closest to tip of pinch --> determines separation point for subgrid A and subgrids B and C.      
    dist_pinch = @(node) sqrt((node(:,1)-x_pinch).^2 + (node(:,2)-z_pinch).^2);
    top2pinch = dist_pinch(top_nodes);
    bottom2pinch = dist_pinch(bottom_nodes);
    net_dist = top2pinch + bottom2pinch;
    [~, pin_idx] = min(net_dist);
  
    sep_point_top = top_nodes(pin_idx,:);
    sep_point_bottom = bottom_nodes(pin_idx,:);    
    
    x_top = top_nodes(:,1);
    x_bottom = bottom_nodes(:,1);    

    % --- A: Make catesian subgrid of pXA based on nodes and points to left
    % of separation point ---
    top_nodesA = top_nodes(x_top <= sep_point_top(:,1), :);
    bottom_nodesA = bottom_nodes(x_bottom <= sep_point_bottom(:,1), :);        
    num_x_pXA = size(top_nodesA, 1);

    top_ptsA = [top_pts(top_pts(:,1) <= sep_point_top(:,1), :); sep_point_top]; % polygon data points on top side LEFT of separation point, including top separation point
    bottom_ptsA = [bottom_pts(bottom_pts(:,1) <= sep_point_bottom(:,1), :); sep_point_bottom]; % (:,1) -> compare x-value
    side_ptsA = [top_ptsA; bottom_ptsA];    

    polyA = PinchOuts(all_polys, poly_num, side_ptsA);
    polyA.bottom_side = bottom_ptsA;
    polyA.top_side = top_ptsA;
    polyA = cartesianSubgrid(polyA, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXA);     

    % CHANGE COORDS (pXA)
    [polyA, z_sep_idx] = PinchOuts.coordCorrectionSubgridA(polyA, top_nodesA, bottom_nodesA, pinch);   
   
    % Get logical indices for new polygon
    polyA = logicalIndicesUpperOverlap(polyA, poly_upper);     

    p_idxA = strcat(p_idx, 'A');
    polyA.p_idx = p_idxA;
    poly_obj.(p_idxA) = polyA;
    nodes_overlap.(p_idxA) = top_nodesA; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxA) = top_ptsA;

    % Get nodes of east side (needed when gluing pXB and pXC)
    NxA = polyA.G.cartDims(1)+1; NzA = polyA.G.cartDims(2)+1; % +1 since we are indexing NODES not CELLS
    east_nodesA = polyA.G.nodes.coords(NxA:NxA:NxA*NzA, :);
    east_nodesAC = east_nodesA(1:z_sep_idx, :); % from bottom to pinch
    east_nodesAB = east_nodesA(z_sep_idx:end, :); % from pinch to top


    % --- B: Make cartesian subgrid of pXB. ---
    top_nodesB = top_nodes(x_top >= sep_point_top(:,1), :);       
    [bottom_ptsB, top_ptsC] = PinchOuts.getSides_pinchout(poly_pinch, poly_num_pinch, poly_num); % bottom points of polygon pXB are top points of fault, and top points of pXC are bottom point of fault
    num_x_pXB = size(top_nodesB, 1);  

    top_ptsB = [sep_point_top; top_pts(top_pts(:,1) >= sep_point_top(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsB = [top_ptsB; bottom_ptsB];

    polyB = PinchOuts(all_polys, poly_num, side_ptsB);       
    polyB.bottom_side = bottom_ptsB;
    polyB.top_side = top_ptsB;
    NzB = NzA - z_sep_idx; % subtract num cells below tip of pinch to get num cells for subgrid B
    polyB = cartesianSubgrid(polyB, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXB, NzB);
   
    [bottom_nodesB, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(polyB, bottom_ptsB, G_glob);
    polyB.bottom_side_new = polyB.bottom_side;
    polyB.top_side_new = polyB.top_side;
    polyB.bottom_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(polyB.bottom_side_new, bottom_nodesB, closest_mask, 'linear');
    bottom_nodesB(rem_idx,2) = z_new;

    polyB = PinchOuts.coordCorrectionSubgridBC(polyB, top_nodesB, bottom_nodesB, east_nodesAB);        

    polyB = logicalIndicesUpperOverlap(polyB, poly_upper);      

    p_idxB = strcat(p_idx, 'B');
    polyB.p_idx = p_idxB;
    poly_obj.(p_idxB) = polyB;
    nodes_overlap.(p_idxB) = top_nodesB; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxB) = top_ptsB;


    % --- C: Make cartesian subgrid of pXC. ---
    bottom_nodesC = bottom_nodes(x_bottom >= sep_point_bottom(:,1), :);
    num_x_pXC = num_x_pXB; % enforce same number of columns (not strictly necessary, but ensures same number of nodes for bottom and top side of full polygon subgrid
    NzC = z_sep_idx - 1; % -1 to get num CELLS

    bottom_ptsC = [sep_point_bottom; bottom_pts(bottom_pts(:,1) >= sep_point_bottom(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsC = [top_ptsC; bottom_ptsC];

    polyC = PinchOuts(all_polys, poly_num, side_ptsC);            
    polyC.bottom_side = bottom_ptsC;
    polyC.top_side = top_ptsC;
    polyC = cartesianSubgrid(polyC, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXC, NzC);

    [top_nodesC, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(polyC, top_ptsC, G_glob);
    polyC.bottom_side_new = polyC.bottom_side;
    polyC.top_side_new = polyC.top_side;
    polyC.top_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(polyC.top_side_new, top_nodesC, closest_mask, 'linear');
    top_nodesC(rem_idx,2) = z_new;

    polyC = PinchOuts.coordCorrectionSubgridBC(polyC, top_nodesC, bottom_nodesC, east_nodesAC);       

    polyC = logicalIndicesWestOverlap(polyC, polyA);           

    p_idxC = strcat(p_idx, 'C');
    polyC.p_idx = p_idxC;
    poly_obj.(p_idxC) = polyC;
    nodes_overlap.(p_idxC) = top_nodesC; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxC) = top_ptsC;
end