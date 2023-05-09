function [poly_obj, nodes_overlap, pts_overlap] = glueToUpperPolygon(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)  
    % Create a structured, Cartesian grid for a given polygon layer,
    % whose top surface conforms to bottom surface of its upper neigbhor
    % (assuming grid for upper neighbor has already been generated).
    % Resulting grids are compatible for gluing.
    % 
    % INPUTS:
    %   all_polys: polygonal points for each polygon of FluidFlower
    %   poly_obj: struct with properties of discretized polygon
    %   poly_num: index of polygon to discretize
    %   poly_num_upper: index of the polygon's upper neighbor
    %   nodes_overlap: struct of overlapping nodes used for gluing grids
    %   pts_overlap: struct of overlapping original polygonal points
    %   G_glob: global, virtual background grid for FluidFlower
    %
    % RETURNS:
    %   poly_obj: updated poly_obj struct
    %   nodes_overlap: updated nodes_overlap struct
    %   pts_overlap: updated pts_overlap struct
    %
    fac_scale = 1;
    if nargin == 8
        inter_horz = varargin{1};
    elseif nargin > 8
        inter_horz = varargin{1};
        fac_scale = varargin{2};
    else        
        inter_horz = true;
    end

    p_idx = strcat('p', string(poly_num));
    
    p_name = fieldnames(poly_obj); 
    is_pinch = 0;
    for i=1:numel(p_name)
        is_pinch = is_pinch + ~isempty(regexp(p_name{i}, strcat('p', string(poly_num_upper), '[ABC]'), 'match'));
    end

    if poly_num == 2 % handle this separately
        p_idx_upper = {strcat('p', string(poly_num_upper), 'left'), ...
                        strcat('p', string(poly_num_upper), 'mid'), ...
                        strcat('p', string(poly_num_upper), 'right')};
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2}), ...
                       poly_obj.(p_idx_upper{3})};

    elseif is_pinch
        p_idx_upper = {strcat('p', string(poly_num_upper), 'A'), ...
                        strcat('p', string(poly_num_upper), 'C')}; % subgrid A and C together comprise the bottom side of upper polygon
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2})};
    else
        p_idx_upper = strcat('p', string(poly_num_upper));
        poly_upper = poly_obj.(p_idx_upper);
    end

    poly_obj.(p_idx) = PolygonGrid(all_polys, poly_num);           
    poly = poly_obj.(p_idx); % dummy variable for readability    
    poly.p_idx = p_idx;    

    for i=1:numel(poly_upper)
        % if points have been added to neighboring polygon, make sure to
        % add these to current polygon
        if ~iscell(poly_upper)
            poly_i = poly_upper;
        else
            poly_i = poly_upper{i};
        end
        if ~isempty(poly_i.p_added)
            n_added = size(poly_i.p_added, 1);
            poly.p(end+1:end+n_added, :) = poly_i.p_added;
        end
    end

    Lx = max(G_glob.faces.centroids(:,1));
    Lz = max(G_glob.faces.centroids(:,2));
    N = G_glob.cartDims;
    Nx = N(1); Nz = N(2);   
    
    if numel(poly_upper) > 1 % Upper neighbor contains a pinch-out
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');
    else
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = findOverlappingNodes(poly, poly_upper, 'top'); 
    end
    
    num_x_overlap = size(nodes_overlap.(p_idx), 1);

    [poly, split_pts] = reorderPts(poly, poly_num);
    
    [top_side, bottom_side] = topAndBottomSurfaces(poly, split_pts, poly_num, pts_overlap); % poly_num required for edge-cases
    
    poly.top_side = unique(top_side, 'rows');
    poly.bottom_side = unique(bottom_side, 'rows');
    
    poly = cartesianSubgrid(poly, Lx, Lz, Nx, Nz, num_x_overlap, [], fac_scale);                
    
    poly.G.nodes.coords(poly.top_mask, :) = nodes_overlap.(p_idx); % glue to overlapping nodes on top side
    
    [poly, closest_bottom] = correct_dx_poly(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
       
    poly = interpolateZ_remaining(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'linear');
           
    poly = interpolateInternal(poly, poly.top_mask, poly.bottom_mask, []);

    if ~isempty(poly.p_we)
        poly = interpolateSide(poly);
    end
    % Interpolate x-points in internal to conform with shifted boundary                         
    if inter_horz
        poly = interpolateHorizontal(poly);
    else
        poly = interpolatePartlyHorizontal(poly, 0.2); % on west and east side, interpolate 20% of horizontal extent
    end

    % Get logical indices for new polygon   
    poly = logicalIndicesUpperOverlap(poly, poly_upper);  

    poly_obj.(p_idx) = poly;
end
