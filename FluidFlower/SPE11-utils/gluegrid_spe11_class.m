%% Load geometry
filename = 'draft_spe11_with_facies_markers.geo';
[pts, loops, facies] = parse_spe11_geo(filename);
% pts (287 x 3): 3D coords of each Point
% loops (32 x variable): ordered array of points forming the surface of loop i
% facies (7 x variable): cell array of surfaces (set of loops) connected to each facie 

%% Other definitions
poly = struct; % to hold each polynomial
split_pts = struct; % to store split points of each poly
pts_overlap = struct;
nodes_overlap = struct;

logicalBottom = @(G) G.nodes.coords(:,2) == min(G.nodes.coords(:,2));
logicalTop = @(G) G.nodes.coords(:,2) == max(G.nodes.coords(:,2));
facesBottom = @(G) G.faces.centroids(:,2) == min(G.faces.centroids(:,2));
facesTop = @(G) G.faces.centroids(:,2) == max(G.faces.centroids(:,2));

gravity reset on

%% Global grid
[G_glob, x_glob, z_glob] = PolygonGrid.globalCartGrid(pts, 200, 150); % nx=200, nz=150

Lx = max(G_glob.faces.centroids(:,1));
Lz = max(G_glob.faces.centroids(:,2));
N = G_glob.cartDims;
Nx = N(1); Nz = N(2);

%% Organize polygons
polys = PolygonGrid.polyPoints(pts, loops, facies);
% Upper three surfaces: 28, 30, 32
poly.top = PolygonGrid(polys, 32);

ptop_orig = poly.top.p_orig;
ptop = poly.top.p;

%% Find top and bottom surfaces
[poly.top, split_pts.top] = reorderPts(poly.top, 32);

[top_side, bottom_side] = topAndBottomSurfaces(poly.top, split_pts.top, 32, pts_overlap);
poly.top.top_side = top_side;
poly.top.bottom_side = bottom_side;

poly.top = cartesianSubgrid(poly.top, Lx, Lz, Nx, Nz, []); % []: no conforming of horizontal grid size

%% dx-correction
% For each surface point, find closest x-coord in subgrid and change it to
% equal the surface coordinate
poly.top.bottom_mask = logicalBottom(poly.top.G); % only select nodes at bottom of grid
poly.top.top_mask = logicalTop(poly.top.G);
poly.top.G = computeGeometry(poly.top.G);
poly.top.faces_bottom = facesBottom(poly.top.G);
poly.top.faces_top = facesTop(poly.top.G);
% Bottom:
[poly.top, closest_bottom] = correct_dx_poly(poly.top, G_glob, poly.top.bottom_mask, 'bottom');
% Top:
[poly.top, closest_top] = correct_dx_poly(poly.top, G_glob, poly.top.top_mask, 'top');

%% Interpolation
% Interpolate z-values at surface of polygon
poly.top = interpolateZ_remaining(poly.top, closest_bottom, ...
                                poly.top.bottom_mask, 'bottom', 'spline'); % top face is the "remaining" parts of the polygon we want bottom surface
poly.top = interpolateZ_remaining(poly.top, closest_top, ...
                                poly.top.top_mask, 'top', 'spline');

% Interpolate x+z values in interior
poly.top = interpolateInternal(poly.top, poly.top.top_mask, poly.top.bottom_mask, []);

%% Finally, update grid
[ii, jj] = gridLogicalIndices(poly.top.G);
poly.top.G.i = ii;
poly.top.G.j = max(jj) - jj + 1; % To make lowest index start from top and increase downwards

Gtop = poly.top.G;
ptop_bottom = poly.top.bottom_side;
ptop_top = poly.top.top_side;

% Rename to polygon idx
poly = cell2struct(struct2cell(poly), {'p32'});

%% Create subgrid for neighboring polygon below
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 30, 32, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 18, 30, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 13, 18, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 28, 32, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 20, 28, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 15, 20, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 6, 15, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Create subgrids around pinch-outs
[poly, nodes_overlap, pts_overlap] = gluePinchOut(polys, poly, 5, 13, 12, ...   
                                                     nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = gluePinchOut(polys, poly, 29, 32, 25, ...
                                                     nodes_overlap, pts_overlap, G_glob);
% For polygon 19, upper neighbor contains pinch-out
[poly, nodes_overlap, pts_overlap] = gluePinchOut(polys, poly, 19, 29, 25, ...
                                                     nodes_overlap, pts_overlap, G_glob);

%% Continue downwards until polygon 11 -> this one must be implemented manually
%For polygon 14, upper neighbor contains pinch-out -> handled internally in
%glueToUpperPolygon
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 14, 19, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 4, 14, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 27, 6, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 17, 6, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Handle polygon 11
[poly, nodes_overlap, pts_overlap] = gluePolygon11(polys, poly, 11, [17,4,5], ...
                                                     nodes_overlap, pts_overlap, G_glob);

%% Continue downwards on both sides of bottom fault
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 16, 27, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 10, 16, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 8, 10, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 3, 8, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Remaining polygons are glued to multiple sub-polygons
[poly, nodes_overlap, pts_overlap] = glueUpperPolygonMultiple(polys, poly, 7, 11, ...
                                                    nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueUpperPolygonMultiple(polys, poly, 2, 7, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Bottom polygon glued to structured polygon and fault
[poly, nodes_overlap, pts_overlap] = gluePolygon1(polys, poly, 1, [3,2], ...
                                                   nodes_overlap, pts_overlap, G_glob);

%% PEBI - top left fault
poly31_neighbors = cell(11, 2);
poly31_neighbors{1,1} = poly.p32; poly31_neighbors{1,2} = 'top';
poly31_neighbors{2,1} = poly.p28; poly31_neighbors{2,2} = 'west';
poly31_neighbors{3,1} = poly.p20; poly31_neighbors{3,2} = 'west';
poly31_neighbors{4,1} = poly.p15; poly31_neighbors{4,2} = 'west';
poly31_neighbors{5,1} = poly.p6; poly31_neighbors{5,2} = 'west';
poly31_neighbors{6,1} = poly.p17; poly31_neighbors{6,2} = 'west';
poly31_neighbors{7,1} = poly.p11left; poly31_neighbors{7,2} = 'west';
poly31_neighbors{8,1} = poly.p4; poly31_neighbors{8,2} = 'east';
poly31_neighbors{9,1} = poly.p14; poly31_neighbors{9,2} = 'east';
poly31_neighbors{10,1} = poly.p19A; poly31_neighbors{10,2} = 'east';
poly31_neighbors{11,1} = poly.p29A; poly31_neighbors{11,2} = 'east';

[all_nodes, nodes_overlap, pts_overlap] = glueFault_topleft(polys, poly, 31, poly31_neighbors, ...
                                                        nodes_overlap, pts_overlap, G_glob);


%% Plots / tests
figure()
poly_idxs = fieldnames(poly);

overlap_idxs = fieldnames(nodes_overlap);
for i = 1:numel(poly_idxs) % plot subgrids first
    px = poly_idxs{i};
    if ~isempty(poly.(px).G)
        plotGrid(poly.(px).G)
        hold on
    end
end
for i = 1:numel(overlap_idxs) % then plot overlapping nodes
    px = overlap_idxs{i};
    plot(nodes_overlap.(px)(:,1), nodes_overlap.(px)(:,2), 'r.', 'markersize', 10)
    hold on
    plot(pts_overlap.(px)(:,1), pts_overlap.(px)(:,2), 'b.', 'markersize', 15)
    hold on
%     plot(poly.(px).bottom_side(:,1), poly.(px).bottom_side(:,2), 'b.', 'markersize', 15)
%     hold on
end

plot(all_nodes(:,1), all_nodes(:,2), 'g.', 'markersize', 10)

xlim([0, 2.8])
ylim([0, 1.2])

%% Make top polygon basis for gluing
assert(strcmp(poly_idxs{1}, 'p32')) % make sure top subgrid is first element
poly_idxs = poly_idxs(2:end);

%% Nested gluing
poly.glued = poly.p32; % top grid is basis for gluing

% NB: Gluing doesn't work until faults are set up!
for i=1:numel(poly_idxs)    
    px = poly_idxs{i};
    disp(px)
    G_new = poly.(px).G;
    G_glued = glue2DGrid(poly.glued.G, G_new);

    G_glued.cells.indexMap = [poly.glued.G.cells.indexMap; G_new.cells.indexMap + poly.glued.G.cells.num];
    G_glued.i = [poly.glued.G.i; G_new.i];
    G_glued.j = [poly.glued.G.j; G_new.j];

    poly.glued.G = G_glued;  
end


%% FUNTCION DEFINTIONS
function [poly_obj, nodes_overlap, pts_overlap] = glueToUpperPolygon(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)
    if nargin > 7
        plot_overlap = varargin{1};
    else
        plot_overlap = false;
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
        %[nodes_overlap.(p_idx), pts_overlap.(p_idx)] = PinchOuts.findOverlappingNodesPinch(poly, poly_upper{1}, poly_upper{2}, 'top');
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');
    else
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = findOverlappingNodes(poly, poly_upper, 'top'); 
    end
    
    num_x_overlap = size(nodes_overlap.(p_idx), 1);

    [poly, split_pts] = reorderPts(poly, poly_num); % two last args required for edge-cases
    
    [top_side, bottom_side] = topAndBottomSurfaces(poly, split_pts, poly_num, pts_overlap); % poly_num required for edge-cases
    
    poly.top_side = unique(top_side, 'rows');
    poly.bottom_side = unique(bottom_side, 'rows');
    
    poly = cartesianSubgrid(poly, Lx, Lz, Nx, Nz, num_x_overlap);   

    if plot_overlap
        figure()
        plotGrid(poly.G)
        hold on
        for i=1:numel(poly_upper)
            plotGrid(poly_upper{i}.G)           
            hold on
        end
        plot(nodes_overlap.(p_idx)(:,1), nodes_overlap.(p_idx)(:,2), 'r.', 'markersize', 10)
        hold on
        plot(pts_overlap.(p_idx)(:,1), pts_overlap.(p_idx)(:,2), 'b.', 'markersize', 15)
        title('Overlapping geometry points for neighboring polygons')
    end                       
    
    poly.G.nodes.coords(poly.top_mask, :) = nodes_overlap.(p_idx); % glue to overlapping nodes on top side
    
    [poly, closest_bottom] = correct_dx_poly(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
    
    poly = interpolateZ_remaining(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'spline');
           
    poly = interpolateInternal(poly, poly.top_mask, poly.bottom_mask, []);

    if ~isempty(poly.p_we)
        poly = interpolateSide(poly);
    end
    % Interpolate x-points in internal to conform with shifted boundary                         
    poly = interpolateHorizontal(poly);

    % Get logical indices for new polygon
    if poly_num == 2
        poly = logicalIndicesUpperOverlapNew(poly, poly_upper);
    else
        poly = logicalIndicesUpperOverlap(poly, poly_upper);
    end

    % --- Is this needed? ---
    poly.faces_bottom = poly.G.faces.centroids(:,2) == min(poly.G.faces.centroids(:,2)); 
    poly.faces_top = poly.G.faces.centroids(:,2) == max(poly.G.faces.centroids(:,2));
    % ---
    
    poly_obj.(p_idx) = poly;
end


function [poly_obj, nodes_overlap, pts_overlap] = gluePinchOut(all_polys, poly_obj, poly_num, poly_num_upper, poly_num_pinch, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)   
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
    top_side = nodes_overlap_pX;
    top_pts = pts_overlap_pX;
    num_x_overlap = size(top_side, 1);
        
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

    % Find node in top side with x-coord closest to x_pin
    [~, pin_idx] = min(abs(top_side(:,1) - x_pinch));
    sep_point_top = top_side(pin_idx,:);
    pin_idx = nnz(top_side(top_side(:,1) <= sep_point_top(:,1), 1)); % this is index for separation point at bottom    
    
    [bottom_side, closest_mask] = PinchOuts.separationPoint_pinchout(poly, bottom_pts, G_glob);

    [z_new, rem_idx] = PinchOuts.interpolateZSide(bottom_pts, bottom_side, closest_mask, 'spline');
    bottom_side(rem_idx,2) = z_new;
    sep_point_bottom = bottom_side(pin_idx,:);

    % Using separation points, divide polygon into pXA, pXB and pXC (for
    % polygon of number X):
    x_top = top_side(:,1);
    x_bottom = bottom_side(:,1);    


    % --- A: Make catesian subgrid of pXA based on nodes and points to left
    % of separation point ---
    top_nodesA = top_side(x_top <= sep_point_top(:,1), :);
    bottom_nodesA = bottom_side(x_bottom <= sep_point_bottom(:,1), :);        
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
    polyA = interpolateHorizontal(polyA);
   
    % Get logical indices for new polygon
    polyA = logicalIndicesUpperOverlap(polyA, poly_upper);     

    p_idxA = strcat(p_idx, 'A');
    polyA.p_idx = p_idxA;
    poly_obj.(p_idxA) = polyA;
    nodes_overlap.(p_idxA) = top_nodesA; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxA) = top_ptsA;

    % Get nodes of east side (needed when gluing pXB and pXC)
    NxA = polyA.G.cartDims(1)+1; NzA = polyA.G.cartDims(2)+1; % +1 since we are indexing NODES not CELLS
    east_nodesA = polyA.G.nodes.coords(NxA:NxA:NxA*NzA, :); % correct indexing ???
    east_nodesAC = east_nodesA(1:z_sep_idx, :); % from bottom to pinch
    east_nodesAB = east_nodesA(z_sep_idx:end, :); % from pinch to top


    % --- B: Make cartesian subgrid of pXB. ---
    top_nodesB = top_side(x_top >= sep_point_top(:,1), :);       
    [bottom_ptsB, top_ptsC] = PinchOuts.getSides_pinchout(poly_pinch, poly_num_pinch, poly_num); % bottom points of polygon pXB are top points of fault, and top points of pXC are bottom point of fault
    num_x_pXB = size(top_nodesB, 1);  

    top_ptsB = [sep_point_top; top_pts(top_pts(:,1) >= sep_point_top(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsB = [top_ptsB; bottom_ptsB];

    polyB = PinchOuts(all_polys, poly_num, side_ptsB);       
    polyB.bottom_side = bottom_ptsB;
    polyB.top_side = top_ptsB;
    NzB = NzA - z_sep_idx; % NB: number of CELLS not NODES
    polyB = cartesianSubgrid(polyB, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXB, NzB);
  
    [bottom_nodesB, closest_mask] = PinchOuts.separationPoint_pinchout(polyB, bottom_ptsB, G_glob);

    [z_new, rem_idx] = PinchOuts.interpolateZSide(bottom_ptsB, bottom_nodesB, closest_mask, 'spline');
    bottom_nodesB(rem_idx,2) = z_new;

    polyB = PinchOuts.coordCorrectionSubgridBC(polyB, top_nodesB, bottom_nodesB, east_nodesAB);    
    polyB = interpolateHorizontal(polyB);

    polyB = logicalIndicesUpperOverlap(polyB, poly_upper);      

    p_idxB = strcat(p_idx, 'B');
    polyB.p_idx = p_idxB;
    poly_obj.(p_idxB) = polyB;
    nodes_overlap.(p_idxB) = top_nodesB; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxB) = top_ptsB;


    % --- C: Make cartesian subgrid of pXC. ---
    bottom_nodesC = bottom_side(x_bottom >= sep_point_bottom(:,1), :);
    num_x_pXC = num_x_pXB; % enforce same number of columns (not strictly necessary, but ensures same number of nodes for bottom and top side of full polygon subgrid
    NzC = z_sep_idx - 1; % -1 to remove end-node

    bottom_ptsC = [sep_point_bottom; bottom_pts(bottom_pts(:,1) >= sep_point_bottom(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsC = [top_ptsC; bottom_ptsC];

    polyC = PinchOuts(all_polys, poly_num, side_ptsC);            
    polyC.bottom_side = bottom_ptsC;
    polyC.top_side = top_ptsC;
    polyC = cartesianSubgrid(polyC, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXC, NzC);

    [top_nodesC, closest_mask] = PinchOuts.separationPoint_pinchout(polyC, top_ptsC, G_glob);

    [z_new, rem_idx] = PinchOuts.interpolateZSide(top_ptsC, top_nodesC, closest_mask, 'spline');
    top_nodesC(rem_idx,2) = z_new;

    polyC = PinchOuts.coordCorrectionSubgridBC(polyC, top_nodesC, bottom_nodesC, east_nodesAC);    
    polyC = interpolateHorizontal(polyC);

    polyC = logicalIndicesWestOverlap(polyC, polyA);           

    p_idxC = strcat(p_idx, 'C');
    polyC.p_idx = p_idxC;
    poly_obj.(p_idxC) = polyC;
    nodes_overlap.(p_idxC) = top_nodesC; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxC) = top_ptsC;
end

