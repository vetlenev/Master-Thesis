classdef PinchOuts < PolygonGrid

    methods
        function obj = PinchOuts(polys, poly_num, side_points) % obj = poly_num
            %POLYGONGRID Discretized geometry of polygon with id poly_num
            % Input:
            %   polys: polygon data from spe11 case
            %   poly_num: id of polygon to discretize                       
            obj = obj@PolygonGrid(polys, poly_num);            
            obj.p = side_points;
%             obj.p = polys{poly_num}{2}(:,1:2);
%             obj.p_orig = obj.p;
%             obj.p = obj.p(1:end-1,:);
        end      
       
    end


   methods (Static)
       
       function [poly_side, closest_mask] = separationPoint_pinchout(poly, poly_pts, G_glob)
            Nx = poly.G.cartDims(1);
            xs = poly_pts(:,1);
            zs = poly_pts(:,2);
            poly_side = zeros(Nx+1, 2);
            if Nx < G_glob.cartDims(1) % always the case for polygons with pinch-outs
                x_min = min(xs);
                x_max = max(xs);
                dx = (x_max - x_min)/Nx;
                % Uniform distribution of nodes:
                poly_side(:,1) = x_min + cumsum(repmat(dx, Nx+1, 1)) - dx;        
            end
        
            num_poly_pts = numel(xs);
            closest_mask = zeros(num_poly_pts, 1);
            for i=1:num_poly_pts
                x_diff = abs(xs(i) - poly_side(:,1));
                closest_mask(i) = find(x_diff == min(x_diff)); % & boundary mask)
            end
          
            poly_side(closest_mask, 1) = xs;
            poly_side(closest_mask, 2) = zs;
            x_side = poly_side(:,1);
        
            if min(xs) > min(x_side) % left poly-boundary point is INSIDE subgrid -> change boundary of subgrid to conform with this
                [xs_sort, sort_idx] = sort(xs);
                zs_sort = zs(sort_idx);
        
                x_stop_mod = xs_sort(2); % x-coordinate to STOP modifications
                mod_idx = x_side <= x_stop_mod;
                mod_num = nnz(mod_idx);
                mod_dx = (xs_sort(2)-xs_sort(1))/(mod_num-1); % -1 to divide by correct number of intervals
                x_mod = xs_sort(1) + cumsum(repmat(mod_dx, mod_num, 1)) - mod_dx;
        
                % Change coordinates and update closest mask for LEFTMOST point:
                poly_side(mod_idx,1) = x_mod;
                x_side = poly_side(:,1); % updated
                x_diff = abs(xs_sort(1) - x_side);
                closest_mask(sort_idx(1)) = find(x_diff == min(x_diff));
                poly_side(closest_mask(sort_idx(1)), 2) = zs_sort(1); % conform z-value of end-node with depth of surface point
            end
        
            if max(xs) < max(x_side) % right poly-boundary point is INSIDE subgrid -> change boundary of subgrid to conform with this
                [xs_sort, sort_idx] = sort(xs);
                zs_sort = zs(sort_idx);
        
                x_start_mod = xs_sort(end-1); % x-coordinate to START modifications
        
                mod_idx = x_side >= x_start_mod;
                mod_num = nnz(mod_idx);
                mod_dx = (xs_sort(end)-xs_sort(end-1))/(mod_num-1); 
                x_mod = xs_sort(end-1) + cumsum(repmat(mod_dx, mod_num, 1)) - mod_dx;
        
                % Change coordinates and update closest mask for RIGHTMOST point:
                poly_side(mod_idx, 1) = x_mod;
                x_side = poly_side(:,1); % updated
                x_diff = abs(xs_sort(end) - x_side); % x_diff = abs(xs(end) - x_sub);
                closest_mask(sort_idx(end)) = find(x_diff == min(x_diff));
                poly_side(closest_mask(sort_idx(end)), 2) = zs_sort(end);
            end
        
            %x_pinch = pinch(:,1);
            %[~, pin_idx] = min(abs(poly_side(:,1) - x_pinch));
            %xb = nnz(poly_side(poly_side(:,1) <= poly_side(pin_idx,1), 1));
            %xt = nnz(top_side(top_side(:,1) <= sep_point_top(:,1), 1));
            %sep_point = poly_side(pin_idx,:);
       end

       function [z_new, remaining_idx] = interpolateZSide(poly_side, nodes_side, closest_idx, interp_method)           
            poly_pts_idx = zeros(size(nodes_side,1), 1);
            poly_pts_idx(closest_idx) = 1; % logical index for nodes closest to poly points
            %poly_pts_idx(pin_idx) = 1; % node closest to tip of pinch-out
            x = poly_side(:,1);
            z = poly_side(:,2);
            remaining_idx = ~poly_pts_idx; % all points on face not equal to poly points
            x_new = nodes_side(remaining_idx,1); % remaining bottom points -> needs to be interpolated
        
            z_new = interp1(x, z, x_new, interp_method);     
       end

       function [polyA, z_sep_idx] = coordCorrectionSubgridA(polyA, top_nodesA, bottom_nodesA, pinch)
           % Correct and interpolate coordinates of subgrid to the left of
           % the pinch-out.
           % Input:
           %    polyA: instance of PolygonGrid representing subgrid to left
           %            of pinch-out,
           %    top_nodesA: nodes at top surface of polyA
           %    bottom_nodesA: nodes at bottom surface of polyA
           %    pinch: coordinate of pinch-out
           % Output:
           %    polyA: polyA with updated coordinates
           %    z_sep_idx: local vertical index of separation point
            % 1. Change top and bottom sides            
            polyA.G.nodes.coords(polyA.bottom_mask, :) = bottom_nodesA;
            polyA.G.nodes.coords(polyA.top_mask, :) = top_nodesA;
        
            % 2. Change interior points
            polyA = interpolateInternal(polyA, polyA.top_mask, polyA.bottom_mask, []);
        
            % 3. Change node closest to pinch to have coordinate of pinch
            z_pinch = pinch(2);
            polyA_east = polyA.G.nodes.coords(polyA.east_mask, :);
            [~, z_sep_idx] = min(abs(polyA_east(:,2) - z_pinch));
            polyA_pinch = polyA_east(z_sep_idx, :);
            polyA_pinch_idx = ismember(polyA.G.nodes.coords, polyA_pinch, 'rows');
            polyA.G.nodes.coords(polyA_pinch_idx,:) = pinch;
        
            % 4. Interpolate nodes on east side to conform with pinch-point
            %start_i = polyA.G.cartDims(1)+1; % only interpolate east side
            polyA = interpolateInternal(polyA, polyA.top_mask, polyA.bottom_mask, pinch, z_sep_idx, polyA.G.cartDims(1)+1);
       end

       function polyBC = coordCorrectionSubgridBC(polyBC, top_nodesB, bottom_nodesB, east_nodesA)
            % 1. Change top and bottom sides
            polyBC.G.nodes.coords(polyBC.bottom_mask, :) = bottom_nodesB;
            polyBC.G.nodes.coords(polyBC.top_mask, :) = top_nodesB;
        
            % 2. Change interior points
            polyBC = interpolateInternal(polyBC, polyBC.top_mask, polyBC.bottom_mask, []);

            % 3. Let nodes of west side equal the intersecting nodes of pXA
            polyBC.G.nodes.coords(polyBC.west_mask, :) = east_nodesA;

       end

       function [nodes, pts_overlap] = findOverlappingNodesPinch(poly, poly_A, poly_C, face)
            pp = poly.p;
            pp_A = poly_A.p;
            pp_C = poly_C.p;
            G_A = poly_A.G;
            G_C = poly_C.G;

            if strcmp(face, 'top')
                pmask_A = poly_A.bottom_mask; % Yes, bottom_mask here, because we are choosing the face of the OTHER polygon
                pmask_C = poly_C.bottom_mask;
                %faces_other = poly_other.faces_bottom;
            elseif strcmp(face, 'bottom')
                pmask_A = poly_A.top_mask;
                pmask_C = poly_C.top_mask;
                %faces_other = poly_other.faces_top;
            end

            pts_overlap_A = pp(ismembertol(pp, pp_A, 'ByRows', true), :);
            pts_overlap_C = pp(ismembertol(pp, pp_C, 'ByRows', true), :);
            pts_overlap = unique([pts_overlap_A; pts_overlap_C], 'rows');
            % Here we assume subgid 
            min_overlap = min(pts_overlap_A(:,1));
            max_overlap = max(pts_overlap_C(:,1));
            G_A_mask = G_A.nodes.coords(pmask_A, :);
            G_C_mask = G_C.nodes.coords(pmask_C, :);
            
            nodesA = G_A_mask(G_A_mask(:,1) >= min_overlap & ...
                                 G_A_mask(:,1) <= max_overlap, :);
            nodesC = G_C_mask(G_C_mask(:,1) >= min_overlap & ...
                                 G_C_mask(:,1) <= max_overlap, :);

            nodes = unique([nodesA; nodesC], 'rows'); % to avoid selecting the overlapping node between A and C twice
            
       end

       function [nodes_all, pts_overlap_all] = findOverlappingNodesMultiple(poly, poly_upper, face)
            pp = poly.p;
            %pmask = cell(numel(poly_upper), 1);
            pts_overlap = cell(numel(poly_upper), 1);
            nodes = cell(numel(poly_upper), 1);

            min_overlap = inf;
            max_overlap = 0;

            for i=1:numel(poly_upper)  
                if strcmp(face, 'top')
                    pp_upper = poly_upper{i}.bottom_side;
                elseif strcmp(face, 'bottom')
                    pp_upper = poly_upper{i}.top_side;
                end                
                pts_overlap{i} = pp(ismembertol(pp, pp_upper, 'ByRows', true), :);                
                % Here we assume subgid 
                min_overlap = min(min_overlap, min(pts_overlap{i}(:,1)));
                max_overlap = max(max_overlap, max(pts_overlap{i}(:,1)));               
            end
            
            pts_overlap_all = unique(vertcat(pts_overlap{:}), 'rows');

            for i=1:numel(poly_upper)
                G_upper = poly_upper{i}.G;
                if strcmp(face, 'top')
                    pmask_upper = poly_upper{i}.bottom_mask; % Yes, bottom_mask here, because we are choosing the face of the OTHER polygon                    
                elseif strcmp(face, 'bottom')
                    pmask_upper = poly_upper{i}.top_mask;
                end
                G_upper_mask = G_upper.nodes.coords(pmask_upper, :);                
                
                nodes{i} = G_upper_mask(G_upper_mask(:,1) >= min_overlap & ...
                                     G_upper_mask(:,1) <= max_overlap, :);
            end
            
            nodes_all = unique(vertcat(nodes{:}), 'rows'); % to avoid selecting the overlapping node between A and C twice            
        end

       function [top_pts, bottom_pts] = getSides_pinchout(poly_pinch, poly_num_pinch, poly_num)
            pp = poly_pinch.p;
            switch poly_num_pinch
                case 12
                    top_pts = pp(5:8,:);
                    bottom_pts = pp(1:5,:); % tip of pinch (idx=5) included in both top and bottom sides
                case 25
                    if poly_num == 29
                        top_pts = pp(16:21,:);
                        bottom_pts = pp(11:16,:);
                    elseif poly_num == 19
                        top_pts = pp(6:9,:);
                        bottom_pts = pp(3:6,:);
                    end
            end
       end
   end
end