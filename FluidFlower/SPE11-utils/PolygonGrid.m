classdef PolygonGrid
    %POLYGONGRID: Generates a global composite PEBI grid conforming to the
    %geometry given in the SPE 11 benchmark dataset.
    % The SPE 11 reservoir constitutes 7 facies together forming 32
    % surface patches / polygons. Individual polygons are created as
    % instances of this class and glued together to a neighboring instance.
    
    properties
        facies
        p_idx % index of polygon
        p % data points of polygon
        p_bt % points on top and bottom faces
        p_we % points on west and east faces
        p_added % added polygonal points for gluing
        we % boolean if p_we is on west (true) or east (false) face
        p_orig
        G
        top_side % nodes at top curve of polygon
        bottom_side
        top_mask % logical mask of top nodes
        bottom_mask
        west_mask
        east_mask
        faces_top % faces connected to top_side nodes
        faces_bottom
    end
    
    methods
        function obj = PolygonGrid(polys, poly_num) % obj = poly_num
            %POLYGONGRID Discretized geometry of polygon with id poly_num
            % Input:
            %   polys: polygon data from spe11 case
            %   poly_num: id of polygon to discretize                       

            %obj.poly_num = poly_num;
            obj.facies = polys{poly_num}{1};
            obj.p = polys{poly_num}{2}(:,1:2);
            obj.p_orig = obj.p;
            obj.p = obj.p(1:end-1,:);
            obj.p_added = [];
        end
        
        function obj = cartesianSubgrid(obj, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_nodes_overlap, varargin)
            % Create cartesian subgrid of provided obj polygon.
            % Input:
            %   obj: polygon instance
            %   Lx_glob: horizontal grid size of global grid
            %   Lz_glob: vertical --- | ---
            %   Nx_glob: horizontal grid dimension of global grid
            %   Nz_glob vertical --- | ---
            %   num_nodes_overlap: number of overlapping nodes of intersecting
            %                       neighbor. If empty, assume topmost grid.            
            poly = obj.p;
            Lx_min = min(poly(:,1));
            Lx_max = max(poly(:,1));
            Lx = Lx_max - Lx_min;
            fac = Lx/Lx_glob;
            if ~isempty(num_nodes_overlap)
                Nx = num_nodes_overlap - 1; % NB: number of cells is one less than number of nodes !!
            else
                Nx = ceil(fac * Nx_glob);
            end
            dx_glob = Lx_glob/Nx_glob;
            % --- New ---
            Lx_new = Nx * dx_glob;
            % -----------            
        
            poly_b = obj.bottom_side;
            poly_t = obj.top_side;
            Lz_min = mean(poly_b(:,2)); % 2D -> second element is depth-coord
            Lz_max = mean(poly_t(:,2));
            Lz = Lz_max - Lz_min;

            if nargin > 6
                Nz = varargin{1};
            else
                fac = Lz/Lz_glob;
                Nz = ceil(fac * Nz_glob);
            end
        
            obj.G = cartGrid([Nx, 1, Nz], [Lx_new, 0.01, 1]); % unit depth

            % Make nodes comply with global grid
            g_nodes_xpos = reshape(obj.G.nodes.coords(:,1), Nx+1, []);
            g_nodes_xpos(end,:) = Lx;
            obj.G.nodes.coords(:,1) = g_nodes_xpos(:);
            % Shift nodes to lateral location of facies
            obj.G.nodes.coords(:,1) = obj.G.nodes.coords(:,1) + Lx_min;   

            if obj.G.nodes.coords(end,1) < obj.G.nodes.coords(end-1,1)
                % End-node is not eastmost node, make it rightmost node.
                % Necessary check for polygons containing pinch-out tip.
                Nxx = Nx+1; % num x nodes            
                Nzz = Nz+1; % num z nodes
                eastmost_x = obj.G.nodes.coords(Nxx:Nxx:Nxx*Nzz, 1);
                third_eastmost_x = obj.G.nodes.coords(Nxx-2:Nxx:Nxx*Nzz, 1);
                % let second eastmost ("overshot") node take average value
                % of neighboring nodes
                obj.G.nodes.coords(Nxx-1:Nxx:Nxx*Nzz, 1) = 0.5*(eastmost_x + third_eastmost_x);
            end

            obj.bottom_mask = obj.G.nodes.coords(:,2) == min(obj.G.nodes.coords(:,2));
            obj.top_mask = obj.G.nodes.coords(:,2) == max(obj.G.nodes.coords(:,2));
            obj.west_mask = obj.G.nodes.coords(:,1) == min(obj.G.nodes.coords(:,1));
            obj.east_mask = obj.G.nodes.coords(:,1) == max(obj.G.nodes.coords(:,1));
        end

        function [obj, split_points] = reorderPts(obj, edge_case)
            % Reorder points pts, starting from point with lowest x-value.
            % Returns separate arrays for monotonically increasing/decreasing
            % connected cells.            
            pts = obj.p;
            mean_x = mean(pts, 1);
            mean_x = mean_x(1);

            switch edge_case
                case 20
                    rem_idx = 2:3; % remove points along fault                                        
                case 6
                    rem_idx = 3:5; % 4:5
                case 5
                    rem_idx = 25;
                case 14
                    rem_idx = 5;
                case 4
                    rem_idx = [2:3, 7:9];                                     
                case 27
                    rem_idx = 12;
                case 17
                    rem_idx = 1;
                case 11
                    rem_idx = 1:2; 
                case 8
                    rem_idx = 1;
                case 3
                    rem_idx = 1;
                case 7
                    rem_idx = 23:25;
                case 2
                    rem_idx = 16:20;
                otherwise % no removals
                    rem_idx = [];
            end
            rem_side = pts(rem_idx, 1) < mean_x; % if false, points are on east side, if true they are on west side
            pts_we = pts(rem_idx, :);
            pts(rem_idx, :) = []; 

            obj.p_we = pts_we; % to store points along west and east sides
            obj.we = rem_side;
            
            x_pts = pts(:,1);
            y_pts = pts(:,2);
            % Special case if left and right boundaries are at global boundaries:
            U = unique(x_pts);
            count = histc(x_pts, U);
            x_mode = U(count == max(count));    
        
            bnd_idx = find(x_pts == min(x_mode));
           
            if numel(x_mode) > 1 && max(count) > 1    
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

             %stop = size(pts_new, 1);
             %pts_map = [[1:stop-bnd_idx+1, stop-bnd_idx+2:stop]; [bnd_idx:stop, 1:bnd_idx-1]]';             

            if ismember(edge_case, [2,3,4,10,14]) % [4,6,14]
                obj.p_bt = pts; % don't reorder for these edge-cases
            else
                obj.p_bt = pts_new; % ordered points from top and bottom surface
            end

            obj.p = [pts_new; pts_we];

            % After ordering, split the monotonically increasing/decreasing parts
            sgn_diff = sign(diff(pts_new(:,1)));
            split_points = find(diff(sgn_diff)); % index of sign-change -> this separates the sub-surfaces
            [~, max_ydiff_idx] = max(abs(diff(pts_new(:,2)))); % index of maximum y-difference between points in polygon
            split_points(sgn_diff(split_points) == 0) = []; % signed difference is zero => two points on same boundary => remove these from split-points
            
            split_points = split_points + (max_ydiff_idx ~= split_points); % add 1 to split points if skewed comes BEFORE the split point

            i = 1;
            while i <= numel(split_points)
                if pts_new(split_points(i), 1) == pts_new(split_points(i)+1, 1) && ... % at other boundary of global domain
                        numel(x_mode) == 1 % subgrid is in interior of global grid
                    split_points(i) = [];
                else
                    i = i + 1;
                end
            end
            
            split_points = [0; split_points; size(pts_new,1)];                        
        end

        function [p_top, p_bottom] = topAndBottomSurfaces(obj, split_pts, edge_case, pts_overlap)           
             
                num_splits = numel(split_pts) - 1; % -1 to get number of split intervals
                %pts = obj.p;
                pts = obj.p_bt;                
                p_idx = strcat('p', string(edge_case));
    
                switch edge_case
                    case 14
                        p_top = pts(1:4,:);
                        p_bottom = pts(5:end,:);
                        % NB: original indexing changes after this!
                    case 4
                        p_top = [pts(end, :); pts(1, :)];
                        p_bottom = pts(2:4, :);
%                     case 6
%                         p_top = pts(17:end, :);
%                         p_bottom = pts(4:16, :);
                    case 10
                        p_top = pts(1:5, :);
                        p_bottom = pts(6:end, :);
                    case 3
                        p_top = pts(1:4, :);
                        p_bottom = pts(5:9, :);
                    case 2
                        p_top = pts(21:end, :);
                        p_bottom = pts(1:15, :);                        
                    otherwise
                        p_xsub = cell(num_splits, 1);
                        p_zsub = cell(num_splits, 1);
                        p_sub = cell(num_splits, 1);
                        for i=1:num_splits
                            p_xsub{i} = pts(split_pts(i)+1:split_pts(i+1), 1); % +2: ...+1
                            p_zsub{i} = pts(split_pts(i)+1:split_pts(i+1), 2);
                            p_sub{i} = pts(split_pts(i)+1:split_pts(i+1), :);
                        end
                        
                        mean_depth = cellfun(@mean, p_zsub);
                        [~, top_idx] = max(mean_depth);
                        [~, bottom_idx] = min(mean_depth);
                        p_top = p_sub{top_idx};
                        p_bottom = p_sub{bottom_idx};  
                end

                % Only select poly-points overlapping with subgrid to glue
                % to, not points overlapping with other subgrids.
                if ~isempty(structfun(@numel, pts_overlap))
                    top_poly_pts = ismember(p_top, pts_overlap.(p_idx), 'rows');
                    p_top = p_top(top_poly_pts, :);                
                end
        end

        function [obj, closest_mask] = correct_dx_poly(obj, G_glob, boundary_mask, face)
            % Correct x-coords of subgrid closest to polygonal coordinates of
            % surface.
            % Input:
            %   G: cartesian subgrid
            %   G_glob: global cartesian grid            
            %   boundary_mask: logical mask selecting the set of nodes in subgrid that are
            %                   on desired boundary and are candidates to change.
            %   face: face of polygon to correct nodes for ('top' or 'bottom')
            % Output:
            %   G: subgrid with modified nodes
            %   closest_idx: logical array of x-coords closest to each poly-point
            G = obj.G;

            if strcmp(face, 'bottom')
                poly_side = obj.bottom_side;
            elseif strcmp(face, 'top')
                poly_side = obj.top_side;
            else
                error(strcat('Method not implemented for face: ', face));
            end

            num_poly_pts = size(poly_side, 1);
            closest_mask = zeros(num_poly_pts, 1);            

            xs = poly_side(:,1);
            zs = poly_side(:,2);

            % If lateral size of subgrid equals lateral size of global grid,
            % the remaining nodes are already conform with associated nodes of the
            % global grid -> no need to modify nodes!
            % But if xsize(subgrid) < xsize(global_grid) nodes must be modified.
            Nx = G.cartDims(1);
            if Nx < G_glob.cartDims(1)
                x_min = min(xs);
                x_max = max(xs);
                dx = (x_max - x_min)/Nx;
                % Uniform distribution of nodes:
                G.nodes.coords(boundary_mask, 1) = x_min + cumsum(repmat(dx, Nx+1, 1)) - dx;
            end
        
            x_sub = G.nodes.coords(:,1);
        
            for i=1:num_poly_pts
                %xs(i) = poly_side(i,1);
                %zs(i) = poly_side(i,2);
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
                %x_diff = abs(xs(1) - x_sub);
                x_diff = abs(xs_sort(1) - x_sub);
                closest_mask(sort_idx(1)) = find(x_diff == min(x_diff) & boundary_mask);
                G.nodes.coords(closest_mask(sort_idx(1)), 2) = zs_sort(1); % conform z-value of end-node with depth of surface point
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
        
            obj.G = G; % update object            
        end

        function obj = interpolateZ_remaining(obj, closest_idx, face_idx, face, interp_method)
            % Interpolate z-values at remaining nodes (non-poly pts) of desired face
            % Input:
            %   G: cartesian subgrid
            %   closest_idx: logical mask of nodes in G closest to poly pts
            %   face_idx: logical mask of face in polygon to be interpolated
            %   face: face of polygon to correct nodes for ('top' or 'bottom')
            %   interp_method: interpolation method (e.g. linear, spline, ...)
            % Output:
            %   G: grid with interpolated nodes
            if strcmp(face, 'bottom')
                poly_side = obj.bottom_side;
            elseif strcmp(face, 'top')
                poly_side = obj.top_side;
            else
                error(strcat('Method not implemented for face: ', face));
            end

            x_sub = obj.G.nodes.coords(:,1);
            poly_pts_idx = zeros(obj.G.nodes.num, 1);
            poly_pts_idx(closest_idx) = 1; % logical index for nodes closest to poly points
            x = poly_side(:,1);
            z = poly_side(:,2);
            remaining_idx = ~poly_pts_idx & face_idx; % all points on face not equal to poly points
            x_new = x_sub(remaining_idx); % remaining bottom points -> needs to be interpolated
        
            z_new = interp1(x, z, x_new, interp_method);
        
            obj.G.nodes.coords(remaining_idx, 2) = z_new;           
        end

        function obj = interpolateInternal(obj, top_mask, bottom_mask, pinch, varargin)                        
            Nx = obj.G.cartDims(1) + 1; % +1 to include both endpoints
            Nz = obj.G.cartDims(2) + 1;           

            x_bottom = obj.G.nodes.coords(bottom_mask, 1);
            z_bottom = obj.G.nodes.coords(bottom_mask, 2);
            x_top = obj.G.nodes.coords(top_mask, 1);
            z_top = obj.G.nodes.coords(top_mask, 2);
            
            if ~isempty(pinch) % Only interpolate in vertical for x-index located at pinch
                Nzb = varargin{1};
                [Nzb, sort_idx] = sort(Nzb);
                Nzb_all = [1; Nzb; Nz];
                pinch = pinch(sort_idx, :);
                ix = varargin{2};
                xb_all = [x_bottom(ix); pinch(:,1)];
                xt_all = [pinch(:,1); x_top(ix)];
                zb_all = [z_bottom(ix); pinch(:,2)];
                zt_all = [pinch(:,2); z_top(ix)];
                for i=1:size(pinch, 1)              
                    Nzb_m = Nzb_all(i); % previous Nzb
                    Nzb = Nzb_all(i+1) - Nzb_all(i) + 1;
                    Nzb_p = Nzb_all(i+2); % next Nzb
                    Nzt = Nzb_all(i+2) - Nzb_all(i+1) + 1; % +1 to include pinch-point in upper interval as well
                    xb = xb_all(i); zb = zb_all(i);
                    xt = xt_all(i+1); zt = zt_all(i+1);
                    xp = pinch(i,1); zp = pinch(i,2);
%                     xb = x_bottom(ix); xt = x_top(ix);
%                     xp = pinch(1);
%                     zb = z_bottom(ix); zt = z_top(ix);
%                     zp = pinch(2);

                    dx_b = (xp - xb)/(Nzb-1); dz_b = (zp - zb)/(Nzb-1);
                    dx_t = (xt - xp)/(Nzt-1); dz_t = (zt - zp)/(Nzt-1);
                    x_interp_b = xb + cumsum(repmat(dx_b, Nzb, 1)) - dx_b;
                    z_interp_b = zb + cumsum(repmat(dz_b, Nzb, 1)) - dz_b;
                    x_interp_t = xp + cumsum(repmat(dx_t, Nzt, 1)) - dx_t;
                    z_interp_t = zp + cumsum(repmat(dz_t, Nzt, 1)) - dz_t;
                   
%                     obj.G.nodes.coords(ix:Nx:Nx*Nzb, 1) = x_interp_b;
%                     obj.G.nodes.coords(ix:Nx:Nx*Nzb, 2) = z_interp_b;
%                     obj.G.nodes.coords(ix+Nx*(Nzb-1):Nx:Nx*Nz, 1) = x_interp_t;
%                     obj.G.nodes.coords(ix+Nx*(Nzb-1):Nx:Nx*Nz, 2) = z_interp_t;
                    lower_slice = ix+Nx*(Nzb_m-1) : Nx : Nx*Nzb_all(i+1);
                    upper_slice = ix+Nx*(Nzb_all(i+1)-1) : Nx : Nx*Nzb_p;
                    obj.G.nodes.coords(lower_slice, 1) = x_interp_b;
                    obj.G.nodes.coords(lower_slice, 2) = z_interp_b;
                    obj.G.nodes.coords(upper_slice, 1) = x_interp_t;
                    obj.G.nodes.coords(upper_slice, 2) = z_interp_t;
                end

            else % Interpolate in x+z direction for each top-bottom pair
                for i=1:Nx
                    xb = x_bottom(i); xt = x_top(i);
                    zb = z_bottom(i); zt = z_top(i);
                    dx = (xt - xb)/(Nz-1);
                    dz = (zt - zb)/(Nz-1);
                    x_interp = xb + cumsum(repmat(dx, Nz, 1)) - dx;
                    z_interp = zb + cumsum(repmat(dz, Nz, 1)) - dz;
    
                    obj.G.nodes.coords(i:Nx:Nx*Nz, 1) = x_interp; % extract correct nodes for interpolation
                    obj.G.nodes.coords(i:Nx:Nx*Nz, 2) = z_interp;
                end
            end
        end

        function obj = interpolateHorizontal(obj)
            Nx = obj.G.cartDims(1)+1;
            Nz = obj.G.cartDims(2)+1;

            for i=2:Nz-1
                %side_point = obj.G.nodes.coords(ix+Nx*(i-1), :); % (fix_idx - 1) -> (i-1)
                west_point = obj.G.nodes.coords(1+Nx*(i-1), :);
                x_west = west_point(1); z_west = west_point(2);
                %east_point = obj.G.nodes.coords((Nx-ix+1)+Nx*(i-1), :); % point on other side with same z-index                 
                east_point = obj.G.nodes.coords(Nx+Nx*(i-1), :);
                x_east = east_point(1); z_east = east_point(2);
    
                dx = (x_east - x_west)/(Nx-1);
                dz = (z_east - z_west)/(Nx-1);
                x_interp = x_west + cumsum(repmat(dx, Nx, 1)) - dx;
                z_interp = z_west + cumsum(repmat(dz, Nx, 1)) - dz;                    
                
                obj.G.nodes.coords(Nx*(i-1)+1:Nx*i, 1) = x_interp;
                %obj.G.nodes.coords(Nx*(fix_idx-1)+1:Nx*fix_idx, 2) = z_interp;
            end
        end

        function obj = interpolateSide(obj)
            Nx = obj.G.cartDims(1)+1;
            Nz = obj.G.cartDims(2)+1;

            % Interpolate points lying on x-plane of fix-point                        
            west_points = obj.p_we.*obj.we;
            west_points = west_points(any(west_points, 2),:);
            fix_idxs = zeros(size(west_points,1), 1);
            ix = 1;
            for i=1:size(west_points,1)
                fix_point = west_points(i,:);                
                                                                  
                side = obj.G.nodes.coords(ix:Nx:Nx*Nz, :);
                [~, fix_idx] = min(abs(side(:,2) - fix_point(2)));                
                fix_idxs(i) = min(max(2, fix_idx), Nz-1);                
            end
            if ~isempty(west_points)
                obj = interpolateInternal(obj, obj.top_mask, obj.bottom_mask, west_points, fix_idxs, ix);
            end

            east_points = obj.p_we.*(~obj.we);
            east_points = east_points(any(east_points, 2),:);
            fix_idxs = zeros(size(east_points,1), 1);
            ix = Nx;
            for i=1:size(east_points,1)
                fix_point = east_points(i,:);               
                                                                  
                side = obj.G.nodes.coords(ix:Nx:Nx*Nz, :);
                [~, fix_idx] = min(abs(side(:,2) - fix_point(2)));
                fix_idxs(i) = fix_idx;                
            end       
            if ~isempty(east_points)
                obj = interpolateInternal(obj, obj.top_mask, obj.bottom_mask, east_points, fix_idxs, ix);
            end                  
        end

        function [nodes, pts_overlap] = findOverlappingNodes(poly, poly_other, face)
            pp = poly.p;           
            G_other = poly_other.G;          
            pp_other = [poly_other.p; poly_other.p_we; poly_other.p_added]; % include ALL poly points of neighbor

            if strcmp(face, 'top')
                pmask_other = poly_other.bottom_mask; % Yes, bottom_mask here, because we are choosing the face of the OTHER polygon
                %pp_other = poly_other.bottom_side;
                x_or_z = 1;
            elseif strcmp(face, 'bottom')
                pmask_other = poly_other.top_mask;
                %pp_other = poly_other.top_side;
                x_or_z = 1;
            elseif strcmp(face, 'west')
                pmask_other = poly_other.east_mask;
                x_or_z = 2;
            elseif strcmp(face, 'east')
                pmask_other = poly_other.west_mask;
                x_or_z = 2;
            end

            pts_overlap = pp(ismembertol(pp, pp_other, 'ByRows', true), :);
            min_overlap = min(pts_overlap(:,x_or_z));
            max_overlap = max(pts_overlap(:,x_or_z));
            G_other_mask = G_other.nodes.coords(pmask_other, :);
            
            nodes = G_other_mask(G_other_mask(:,x_or_z) >= min_overlap-eps & ...
                                 G_other_mask(:,x_or_z) <= max_overlap+eps, :);
            
        end

        function poly = logicalIndicesUpperOverlap(poly, poly_upper)            
            % CONSIDER REMOVING THIS !
            poly.G = computeGeometry(poly.G);
            fc = poly.G.faces.centroids;              

            if numel(poly_upper) > 1 
                G = struct; fcu = struct; N = struct;
                G.A = poly_upper{1}.G;
                fcu.A = G.A.faces.centroids;
                N.A = G.A.faces.neighbors; 
                G.C = poly_upper{2}.G;
                fcu.C = G.C.faces.centroids;
                N.C = G.C.faces.neighbors; 

                cells_A = N.A(ismember(fcu.A, fc, 'rows'), :);
                cells_C = N.C(ismember(fcu.C, fc, 'rows'), :);

                cells_A_overlap = unique(cells_A(cells_A ~= 0));
                cells_C_overlap = unique(cells_C(cells_C ~= 0));

                ii_A_overlap = G.A.i(cells_A_overlap);
                ii_C_overlap = G.C.i(cells_C_overlap);
                jj_A_overlap = G.A.j(cells_A_overlap);
                jj_C_overlap = G.C.j(cells_C_overlap);

                ii_upper_overlap = [ii_A_overlap; ii_C_overlap];
                jj_upper_overlap = [jj_A_overlap; jj_C_overlap];
            else
                G_upper = poly_upper.G;
                G_upper = computeGeometry(G_upper);            
                fc_other = G_upper.faces.centroids;
                N_other = G_upper.faces.neighbors;

                cells_upper = N_other(ismembertol(fc_other, fc, 'ByRows', true), :);            
                cells_upper_overlap = unique(cells_upper(cells_upper ~= 0));        

                % Logical indices of bounding cells from upper subgrid
                ii_upper_overlap = poly_upper.G.i(cells_upper_overlap);
                jj_upper_overlap = poly_upper.G.j(cells_upper_overlap);
            end                       
        
            [ii, jj] = gridLogicalIndices(poly.G);
            jj = max(jj) - jj + 1;
            poly.G.i = ii + min(ii_upper_overlap) - 1;
            poly.G.j = jj + max(jj_upper_overlap); % NB: for glueToLowerPolygon this will not work! Needs to be modified then!
        end

        function poly = logicalIndicesUpperOverlapNew(poly, poly_upper)            
            poly.G = computeGeometry(poly.G);
            fc = poly.G.faces.centroids;   

            ii_upper_overlap = [];
            jj_upper_overlap = [];

            for i=1:numel(poly_upper)
                Gi = poly_upper{i}.G;
                fcui = Gi.faces.centroids;
                Ni = Gi.faces.neighbors;

                cells_i = Ni(ismember(fcui, fc, 'rows'), :);

                cells_i_overlap = unique(cells_i(cells_i ~= 0));

                ii_i_overlap = Gi.i(cells_i_overlap);
                jj_i_overlap = Gi.j(cells_i_overlap);

                ii_upper_overlap = [ii_upper_overlap; ii_i_overlap];
                jj_upper_overlap = [jj_upper_overlap; jj_i_overlap];            
            end                       
        
            [ii, jj] = gridLogicalIndices(poly.G);
            jj = max(jj) - jj + 1;
            poly.G.i = ii + min(ii_upper_overlap) - 1;
            poly.G.j = jj + max(jj_upper_overlap); % NB: for glueToLowerPolygon this will not work! Needs to be modified then!
        end

        function poly = logicalIndicesWestOverlap(poly, poly_west)
            poly.G = computeGeometry(poly.G);
            G_west = poly_west.G;
            G_west = computeGeometry(G_west);
            fc = poly.G.faces.centroids;
            fc_other = G_west.faces.centroids;
            N_other = G_west.faces.neighbors;    
            cells_west = N_other(ismembertol(fc_other, fc, 'ByRows', true), :);            
            cells_west_overlap = unique(cells_west(cells_west ~= 0));
        
            % Logical indices of bounding cells from western subgrid
            ii_west_overlap = poly_west.G.i(cells_west_overlap);
            jj_west_overlap = poly_west.G.j(cells_west_overlap);
        
            [ii, jj] = gridLogicalIndices(poly.G);
            jj = max(jj) - jj + 1;
            poly.G.i = ii + max(ii_west_overlap);
            poly.G.j = jj + min(jj_west_overlap) - 1; % NB: for glueToLowerPolygon this will not work! Needs to be modified then!
        end
    end

      

    methods (Static)
        function polys = polyPoints(pts, loops, facies)
            % Return a cell array of all 32 polyons with associated points
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

        function [G, x_global, z_global] = globalCartGrid(pts, nx, nz)
            x_min = min(pts(:,1)); z_min = min(pts(:,2));
            x_max = max(pts(:,1)); z_max = max(pts(:,2));
            Lx = x_max - x_min;
            Lz = z_max - z_min;
            Nx = nx;
            Nz = nz;
            dx = Lx/Nx;
            dz = Lz/Nz;

            x_global = x_min + cumsum(repmat(dx, Nx+1, 1)) - dx;
            z_global = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;           

            G = cartGrid([Nx, 1, Nz], [Lx, 0.01, Lz]); % global background grid -> basis for gluing
            G = computeGeometry(G);
        end
    end
end

