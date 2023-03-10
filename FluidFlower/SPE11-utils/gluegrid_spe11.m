%% Load geometry
filename = 'draft_spe11_with_facies_markers.geo';
[pts, loops, facies] = parse_spe11_geo(filename);
% pts (287 x 3): 3D coords of each Point
% loops (32 x variable): ordered array of points forming the surface of loop i
% facies (7 x variable): cell array of surfaces (set of loops) connected to each facie

%% Organize polygons
polys = polyPoints(pts, loops, facies);
% Upper three surfaces: 28, 30, 32
% Sealing surface (red, facies 7): 1, 31 (NB: IS THIS WRONG? ARTICLE STATES
% THAT BOTTOM RED IS HIGH-PERM !?)
% Dark blue (facies 1): 7, 8, 9, 32
% Green (facies 5): 2, 3, 4, 5, 6
p1 = polys{1}{2}(:,1:2); % 13
p2 = polys{30}{2}(:,1:2); % 18
figure()
for i=1:size(p1,1)
    plot(p1(i,1), p1(i,2), 'ko')
    text(p1(i,1), p1(i,2), string(i))
    hold on
    %pause(0.5)
end

p1_orig = p1;
p2_orig = p2;
p1 = p1(1:end-1,:); % last row same as first - remove this
p2 = p2(1:end-1,:);

%% Sample grid
%G.nodes.coords = p1;
x_min = min(pts(:,1)); z_min = min(pts(:,2));
x_max = max(pts(:,1)); z_max = max(pts(:,2));
Lx = x_max - x_min;
Lz = z_max - z_min;
Nx = 200; % 200 intervals in x-dir
Nz = 150;
dx = Lx/Nx;
dz = Lz/Nz;
x_global = x_min + cumsum(repmat(dx, Nx+1, 1)) - dx;
z_global = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;

G_glob = cartGrid([Nx, Nz], [Lx, Lz]); % global background grid -> basis for gluing
G_glob = computeGeometry(G_glob);

if 0
    Ls_a = min(p2(:,1)); % s for subgrid
    Ls_b = max(p2(:,2));
    xs_top = p1_bottom(p1_bottom >= Ls_a & p1_bottom <= Ls_b); % Extract x-coords from top subgrid inside the matching surface below
    % Change number of nodes at top face of subgrid to match that of overlying
    % surface:
    num_xs_top = numel(xs_top);    
end

G1 = cartesianSubgrid(p1, Lx, Lz, Nx, Nz); % num_xs_top)

xG1 = G1.nodes.coords(:,1);
Nx1 = G1.cartDims(1) + 1; % +1 to include both endpoints
Nz1 = G1.cartDims(2) + 1;

%% Find top and bottom surfaces
[p1_new, bnd_idx, x_mode] = reorderPts(p1); % p1_new

sgn_diff = sign(diff(p1_new(:,1)));

% sgn_zero = find(sgn_diff == 0);
% if sgn_zero(1) == 1 && sgn_zero(end) == numel(sgn_diff)    
%     p1_end = p1_new(end,:);
%     p1_new(2:end,:) = p1_new(1:end-1,:);
%     p1_new(1,:) = p1_end;
% end

split_points = find(diff(sgn_diff)); % index of sign-change -> this separates the sub-surfaces
[~, max_ydiff_idx] = max(abs(diff(p1_new(:,2))));
split_points(sgn_diff(split_points) == 0) = [];

split_points = split_points + (max_ydiff_idx ~= split_points); % add 1 to split points if skewed comes BEFORE the split point

i = 1;
while i <= numel(split_points)
    if p1_new(split_points(i), 1) == p1_new(split_points(i)+1, 1) && ... % at other boundary of global domain
            numel(x_mode) == 1 % subgrid is in interior of global grid
        split_points(i) = [];
    else
        i = i + 1;
    end
end
num_splits = numel(split_points); % update
% -- IS +1 COREECT HERE ?? --
split_points = [0; split_points; size(p1_new,1)];
% ---------------------------

p1_xsub = cell(num_splits+1, 1);
p1_ysub = cell(num_splits+1, 1);
p1_sub = cell(num_splits+1, 1);
for i=1:num_splits+1
    p1_xsub{i} = p1_new(split_points(i)+1:split_points(i+1), 1); % +2: ...+1
    p1_ysub{i} = p1_new(split_points(i)+1:split_points(i+1), 2);
    p1_sub{i} = p1_new(split_points(i)+1:split_points(i+1), :);
end

mean_depth = cellfun(@mean, p1_ysub);
[~, top_idx] = max(mean_depth);
[~, bottom_idx] = min(mean_depth);
p1_top = p1_sub{top_idx};
p1_bottom = p1_sub{bottom_idx};

%% dx-correction
% For each surface point, find closest x-coord in subgrid and change it to
% equal the surface coordinate
bottom_mask = G1.nodes.coords(:,2) == min(G1.nodes.coords(:,2)); % only select nodes at bottom of grid
top_mask = G1.nodes.coords(:,2) == max(G1.nodes.coords(:,2));
% Bottom:
[G1, closest_bottom] = correct_dx_poly(G1, G_glob, p1_bottom, bottom_mask);
% Top:
[G1, closest_top] = correct_dx_poly(G1, G_glob, p1_top, top_mask);

%% Interpolate z-values
G1 = interpolateZ_remaining(G1, p1_bottom, closest_bottom, ...
                                bottom_mask, 'spline');
G1 = interpolateZ_remaining(G1, p1_top, closest_top, ...
                                top_mask, 'spline');

%% Interpolate internal points
xG1_bottom = G1.nodes.coords(bottom_mask, 1);
zG1_bottom = G1.nodes.coords(bottom_mask, 2);
xG1_top = G1.nodes.coords(top_mask, 1);
zG1_top = G1.nodes.coords(top_mask, 2);
% Interpolate in x+z direction for each top-bottom pair
for i=1:Nx1
    xb = xG1_bottom(i); xt = xG1_top(i);
    zb = zG1_bottom(i); zt = zG1_top(i);
    dx = (xt - xb)/(Nz1-1);
    dz = (zt - zb)/(Nz1-1);
    x = xb + cumsum(repmat(dx, Nz1, 1)) - dx;
    z = zb + cumsum(repmat(dz, Nz1, 1)) - dz;
    G1.nodes.coords(i:Nx1:Nx1*Nz1, 1) = x;
    G1.nodes.coords(i:Nx1:Nx1*Nz1, 2) = z;
end

%% Sample grid - G2
%G.nodes.coords = p1;
G2 = cartesianSubgrid(p2, Lx, Lz, Nx, Nz);

xG2 = G2.nodes.coords(:,1);
Nx2 = G2.cartDims(1) + 1; % +1 to include both endpoints
Nz2 = G2.cartDims(2) + 1;

%% Find top and bottom surfaces
p2_new = reorderPts(p2);
sgn_diff = sign(diff(p2_new(:,1)));
split_points = find(diff(sgn_diff)); % index of sign-change -> this separates the sub-surfaces
num_splits = numel(split_points);
split_points = [-1; split_points; size(p2_new,1)-1];

p2_xsub = cell(num_splits+1, 1);
p2_ysub = cell(num_splits+1, 1);
p2_sub = cell(num_splits+1, 1);
for i=1:num_splits+1
    p2_xsub{i} = p2_new(split_points(i)+2:split_points(i+1)+1, 1);
    p2_ysub{i} = p2_new(split_points(i)+2:split_points(i+1)+1, 2);
    p2_sub{i} = p2_new(split_points(i)+2:split_points(i+1)+1, :);
end

mean_depth = cellfun(@sum, p2_ysub);
[~, top_idx] = max(mean_depth);
[~, bottom_idx] = min(mean_depth);
p2_top = p2_sub{top_idx};
p2_bottom = p2_sub{bottom_idx};

%% dx-correction
% For each surface point, find closest x-coord in subgrid and change it to
% equal the surface coordinate
num_p_bottom = size(p2_bottom, 1);
num_p_top = size(p2_top, 1);

bottom_mask = G2.nodes.coords(:,2) == min(G2.nodes.coords(:,2)); % only select nodes at bottom of grid
top_mask = G2.nodes.coords(:,2) == max(G2.nodes.coords(:,2));
% Bottom:
[G2, closest_bottom] = correct_dx_poly(G2, G_glob, p2_bottom, bottom_mask);
% Top:
[G2, closest_top] = correct_dx_poly(G2, G_glob, p2_top, top_mask);

%% Interpolate z-values
G2 = interpolateZ_remaining(G2, p2_bottom, closest_bottom, ...
                                bottom_mask, 'spline'); % top face is the "remaining" parts of the polygon we want bottom surface
G2 = interpolateZ_remaining(G2, p2_top, closest_top, ...
                                top_mask, 'spline');

%% Interpolate internal points
xG2_bottom = G2.nodes.coords(bottom_mask, 1);
zG2_bottom = G2.nodes.coords(bottom_mask, 2);
xG2_top = G2.nodes.coords(top_mask, 1);
zG2_top = G2.nodes.coords(top_mask, 2);
% Interpolate in x+z direction for each top-bottom pair
for i=1:Nx2
    xb = xG2_bottom(i); xt = xG2_top(i);
    zb = zG2_bottom(i); zt = zG2_top(i);
    dx = (xt - xb)/(Nz2-1);
    dz = (zt - zb)/(Nz2-1);
    x = xb + cumsum(repmat(dx, Nz2, 1)) - dx;
    z = zb + cumsum(repmat(dz, Nz2, 1)) - dz;
    G2.nodes.coords(i:Nx2:Nx2*Nz2, 1) = x;
    G2.nodes.coords(i:Nx2:Nx2*Nz2, 2) = z;
end

%% Swap variables
G2 = G1;
p2_bottom = p1_bottom;
p2_top = p1_top;

%% Finally, plotting!
figure()
plotGrid(G1)
hold on
plot(p1_bottom(:,1), p1_bottom(:,2), 'r.', 'markersize', 20)
hold on
plot(p1_top(:,1), p1_top(:,2), 'b.', 'markersize', 20)
hold on
if 1
    plotGrid(G2)
    hold on
    plot(p2_bottom(:,1), p2_bottom(:,2), 'r.', 'markersize', 30)
    hold on
    plot(p2_top(:,1), p2_top(:,2), 'g.', 'markersize', 30)
end


%% Functions
function polys = polyPoints(pts, loops, facies)
    polys = cell(numel(loops),1);
    for f_ix = 1:numel(facies)
        facie = facies{f_ix};
        
        for poly_ix = facie' % all polys making up the facie           
            poly = pts(loops{poly_ix}, :);
            poly = [poly; poly(1,:)];
            %patch(poly(:,1), poly(:,2), colors{f_ix}); hold on
            polys{poly_ix} = {f_ix, poly}; % also store the associated facies of the polygon
        end
    end
end

function G = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, varargin)
    Lx_min = min(poly(:,1));
    Lx_max = max(poly(:,1));
    Lx = Lx_max - Lx_min;
    fac = Lx/Lx_glob;
    if nargin > 5
        Nx = varargin{1};
    else
        Nx = ceil(fac * Nx_glob);
    end
    dx_glob = Lx_glob/Nx_glob;
    % --- New ---
    Lx_new = Nx * dx_glob;
    % -----------
    %dx = Lx/Nx;
    %x = Lx_min + cumsum(repmat(dx, Nx+1, 1)) - dx;

    Lz_min = min(poly(:,2)); % 2D -> second element is depth-coord
    Lz_max = max(poly(:,2));
    Lz = Lz_max - Lz_min;
    fac = Lz/Lz_glob;
    %Nz = ceil(fac * Nz_glob);
    Nz = 5;
    %dz = Lz/Nz;
    %z = Lz_min + cumsum(repmat(dz, Nz+1, 1)) - dz;

    G = cartGrid([Nx, Nz], [Lx_new, 1]); % unit depth
    % Make nodes comply with global grid
    g_nodes_xpos = reshape(G.nodes.coords(:,1), Nx+1, []);
    g_nodes_xpos(end,:) = Lx;
    G.nodes.coords(:,1) = g_nodes_xpos(:);
    % Shift nodes to lateral location of facies
    G.nodes.coords(:,1) = G.nodes.coords(:,1) + Lx_min;   
end

function [pts_new, bnd_idx, x_mode] = reorderPts(pts)
    % Reorder points pts, starting from point with lowest x-value.
    % Returns separate arrays for monotonically increasing/decreasing
    % connected cells.
    x_pts = pts(:,1);
    y_pts = pts(:,2);
    % Special case if left and right boundaries are at global boundaries:
    U = unique(x_pts);
    count = histc(x_pts, U);
    x_mode = U(count == max(count));    

    bnd_idx = find(x_pts == min(x_mode));
   
    if numel(x_mode) > 1        
        bnd_idx = bnd_idx(2);    
    elseif numel(bnd_idx) > 1
        if (y_pts(bnd_idx(1)) < y_pts(bnd_idx(2)) && diff(bnd_idx) == 1) || ...
               (y_pts(bnd_idx(1)) > y_pts(bnd_idx(2)) && diff(bnd_idx) > 1) % revolving from bottom to top -> select top as start point
            bnd_idx = bnd_idx(2);
        elseif (y_pts(bnd_idx(1)) < y_pts(bnd_idx(2)) && diff(bnd_idx) > 1) || ...
                (y_pts(bnd_idx(1)) > y_pts(bnd_idx(2)) && diff(bnd_idx) == 1)% revolving from top to bottom
            bnd_idx = bnd_idx(1);
        end
    else
        bnd_idx = bnd_idx(1);
    end

     pts_new = pts;
     pts_start = pts_new(bnd_idx:end,:);
     pts_stop = pts_new(1:bnd_idx-1,:);
     pts_new(1:end-bnd_idx+1,:) = pts_start;
     pts_new(end-bnd_idx+2:end,:) = pts_stop;
    % Separate into monotonically increasing/decreasing parts
    sgn_diff = sign(diff(pts_new(:,1)));
    split_points = find(diff(sgn_diff));
    pts_split = split_points;
end

function [G, closest_mask] = correct_dx_poly(G, G_glob, poly_sub, boundary_mask)
    % Correct x-coords of subgrid closest to polygonal coordinates of
    % surface.
    % Input:
    %   G: cartesian subgrid
    %   G_glob: global cartesian grid
    %   poly_sub: coordinates of desired subpart (top/bottom) of polygon
    %   boundary_mask: logical mask selecting the set of nodes in subgrid that are
    %                   on desired boundary and are candidates to change.
    % Output:
    %   G: subgrid with modified nodes
    %   closest_idx: logical array of x-coords closest to each poly-point
    num_poly_pts = size(poly_sub, 1);
    closest_mask = zeros(num_poly_pts, 1);
    xs = zeros(num_poly_pts, 1);
    zs = zeros(num_poly_pts, 1);

    x_sub = G.nodes.coords(:,1);

    for i=1:num_poly_pts
        xs(i) = poly_sub(i,1);
        zs(i) = poly_sub(i,2);
        x_diff = abs(xs(i) - x_sub);
        x_diff(~boundary_mask) = inf;
        closest_mask(i) = find(x_diff == min(x_diff)); % & boundary mask)
    end

    G.nodes.coords(closest_mask, 1) = xs;
    G.nodes.coords(closest_mask, 2) = zs;
    x_sub = G.nodes.coords(:,1); % updated

    % Modification of nodes needed if poly boundary pts are not on boundary
    % of cartesian subgrid
    if min(xs) > min(x_sub) % left poly-boundary point is INSIDE subgrid -> change boundary of subgrid to conform with this
        [xs_sort, sort_idx] = sort(xs);
        zs_sort = zs(sort_idx);

        x_stop_mod = xs_sort(2); % x-coordinate to STOP modifications
        mod_idx = x_sub <= x_stop_mod & boundary_mask;
        mod_num = nnz(mod_idx);
        mod_dx = (xs_sort(2)-xs_sort(1))/(mod_num-1); % -1 to divide by correct number of intervals
        x_mod = xs_sort(1) + cumsum(repmat(mod_dx, mod_num, 1)) - mod_dx;

        % Change coordinates and update closest mask for LEFTMOST point:
        G.nodes.coords(mod_idx, 1) = x_mod;
        x_sub = G.nodes.coords(:,1); % updated
        x_diff = abs(xs(1) - x_sub);
        closest_mask(1) = find(x_diff == min(x_diff) & boundary_mask);
        G.nodes.coords(closest_mask(1), 2) = zs_sort(1); % conform z-value of end-node with depth of surface point
    end

    if max(xs) < max(x_sub) % right poly-boundary point is INSIDE subgrid -> change boundary of subgrid to conform with this
        [xs_sort, sort_idx] = sort(xs);
        zs_sort = zs(sort_idx);

        x_start_mod = xs_sort(end-1); % x-coordinate to START modifications

        mod_idx = x_sub >= x_start_mod & boundary_mask;
        mod_num = nnz(mod_idx);
        mod_dx = (xs_sort(end)-xs_sort(end-1))/(mod_num-1); 
        x_mod = xs_sort(end-1) + cumsum(repmat(mod_dx, mod_num, 1)) - mod_dx;

        % Change coordinates and update closest mask for RIGHTMOST point:
        G.nodes.coords(mod_idx, 1) = x_mod;
        x_sub = G.nodes.coords(:,1); % updated
        x_diff = abs(xs_sort(end) - x_sub); % x_diff = abs(xs(end) - x_sub);
        closest_mask(sort_idx(end)) = find(x_diff == min(x_diff) & boundary_mask);
        G.nodes.coords(closest_mask(sort_idx(end)), 2) = zs_sort(end);
    end

    % If lateral size of subgrid equals lateral size of global grid,
    % the remaining nodes are already conform with associated nodes of the
    % global grid -> no need to modify nodes!
    % But if xsize(subgrid) < xsize(global_grid) nodes must be modified.
    if G.cartDims(1) == G_glob.cartDims(1)
        return; % skip modifications
    else
        return;
    end
end

function G = interpolateZ_remaining(G, poly_face, closest_idx, face_idx, interp_method)
    % Interpolate z-values at remaining nodes (non-poly pts) of desired face
    % Input:
    %   G: cartesian subgrid
    %   poly_face: face of polygon/surface to interpolate on
    %   closest_idx: logical mask of nodes in G closest to poly pts
    %   face_idx: logical mask of face in polygon to be interpolated
    %   interp_method: interpolation method (e.g. linear, spline, ...)
    % Output:
    %   G: grid with interpolated nodes
    x_sub = G.nodes.coords(:,1);
    poly_pts_idx = zeros(G.nodes.num, 1);
    poly_pts_idx(closest_idx) = 1; % logical index for nodes closest to poly points
    x = poly_face(:,1);
    z = poly_face(:,2);
    remaining_idx = ~poly_pts_idx & face_idx; % all points on face not equal to poly points
    x_new = x_sub(remaining_idx); % remaining bottom points -> needs to be interpolated

    z_new = interp1(x, z, x_new, interp_method);

    G.nodes.coords(remaining_idx, 2) = z_new;
end