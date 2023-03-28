classdef Faults < PolygonGrid

    properties
        fac_scale
        bfacies
        multiple

        internal_top
        internal_bottom
        internal_west
        internal_east

        external_top
        external_bottom
        external_west
        external_east

        G_fault
        gc
        gf
        gn       
    end

    methods
        function obj = Faults(polys, poly_num, fac_scale, multiple_polys) % obj = poly_num
            %POLYGONGRID Discretized geometry of polygon with id poly_num
            % Input:
            %   polys: polygon data from spe11 case
            %   poly_num: id of polygon to discretize  
            %   fac_scale: scale fraction for resolution size (should have
            %   value > 1 to give finer cells inside fault
            obj = obj@PolygonGrid(polys, poly_num);
            obj.fac_scale = fac_scale;
            obj.multiple = multiple_polys;
            if multiple_polys % reset points and facies - must be recomputed
                obj.p = [];
                obj.p_orig = [];
                obj.p_idx = 0;
                obj.facies = 0;                
            end
        end      
       
    end


   methods (Static)              
    
        function [poly, distributed_pts] = QuasiRandomPointDistribution(poly, fault_polys, G_glob, node_density, do_plot, varargin)
            % Make quasi-random distribution of points inside triangle
            % grid, distributed based on fractional size of individual
            % triangles.
            % PARAMETERS:
            %   poly: instance of PolygonGrid
            %   fault_polys: number index for fault to distribute point in
            %   G_glob: global, virtual background grid
            %   node_density: scaling factor to adjust number of additional generating points inside faults
            %   do_plot: whether to plot resulting triangulated grid
            %           (true/false)
            % RETURNS:
            %   poly: updated polygon instance
            %   distributed_pts: coordinated of all distributed points
            A_tot = sum(G_glob.cells.volumes); % total area
            nodes_tot = G_glob.nodes.num;          

            distributed_pts = [];            
            if do_plot
                figure()
                if nargin > 5
                    colors = varargin{1};
                else
                    colors = rand(7,3);
                end
            end
            for i=1:numel(fault_polys)
                p = fault_polys(i);
                if isstring(p)
                    p_idx = fault_polys;
                else
                    p_idx = strcat('p', string(p), 'f'); % 'PEBI'
                end
                poly_internal = poly.(p_idx);
            
                % Original triangulation (only boundary nodes)
                bnodes = poly_internal.bnodes;
                np = size(bnodes,1);
                constraint = [(1:np-1)' (2:np)'; np 1]; 
                dt = delaunayTriangulation(bnodes, constraint); % NB: no constraint for added points!
                inpoly = isInterior(dt);
                dt_inpoly = dt.ConnectivityList(inpoly,:);
                G_dummy = triangleGrid(bnodes, dt_inpoly);
            
                G_int = computeGeometry(G_dummy);
                A_int = sum(G_int.cells.volumes);
                A_frac = A_int/A_tot;
                dim_frac = (G_glob.cartDims(1)*G_glob.cartDims(2))/(200*100); % [200, 100] is "standard" dimensions
                nodes_int = round(node_density*dim_frac*nodes_tot*A_frac);
            
                % loop through triangles
                int_points = [];
                for c=1:G_int.cells.num                    
                    cfaces = G_int.cells.faces(G_int.cells.facePos(c):G_int.cells.facePos(c+1)-1, :);
                    cfaces_areas = G_int.faces.areas(cfaces);
                    %scale = 1 + std(cfaces_areas)/mean(cfaces_areas);
                    scale = node_density*max(cfaces_areas)/min(cfaces_areas); % distribute more points inside very elongated triangles to make them more regularly
                    num_nodes = round(scale*nodes_int*(G_int.cells.volumes(c)/A_int)); % num nodes assigned to this triangle

                    cnodes = zeros(6, 1);
                    for j=1:numel(cfaces)
                        cnodes(2*j-1:2*j) = G_int.faces.nodes(G_int.faces.nodePos(cfaces(j)):G_int.faces.nodePos(cfaces(j)+1)-1, :);
                    end
                    ccoords = unique(G_int.nodes.coords(cnodes, :), 'rows');
    
                    for k=1:num_nodes
                        rand2pts = randperm(3,2);
                        %rand_pt = setdiff(randi(3,num_nodes,1)', repmat((1:3),num_nodes,1)', 'rows');
                        rand_ncoords = ccoords(rand2pts, :);
                        rand_fpt = UtilFunctionsFF.randomPointOnLine(rand_ncoords(1,:), rand_ncoords(2,:), 'uniform');
                        % nonuniform random point between rand_pt and opposite node
                        other_node = setdiff(1:3, rand2pts);
                        other_coord = ccoords(other_node, :);
                        
                        rand_pt = UtilFunctionsFF.randomPointOnLine(rand_fpt, other_coord, 'nonuniform');
                        int_points = [int_points; rand_pt];                       
                    end
                end
                     
                % Triangulation of boundary + quasi-random point distribution
                nodes = [bnodes; int_points];
                dt = delaunayTriangulation(nodes, constraint); % New nodes, same constraint: no constraint for added points!
                inpoly = isInterior(dt);
                dt_inpoly = dt.ConnectivityList(inpoly,:);
                GBF_int = triangleGrid(nodes, dt_inpoly);

                GBF_int = computeGeometry(GBF_int);
                GBF_int.cells.indexMap = (1:GBF_int.cells.num)';
                % Set face-orientation to undefined (NaN)
                GBF_int.cells.faces(:,2) = nan(size(GBF_int.cells.faces(:,1)));
                poly.(p_idx).G = GBF_int;
                % No logical indices for unstructured grid -> set to NaN
                poly.(p_idx).G.i = nan(GBF_int.cells.num, 1);
                poly.(p_idx).G.j = nan(GBF_int.cells.num, 1);
            
                distributed_pts = [distributed_pts; int_points];
                
                if do_plot
                    plotGrid(GBF_int, 'FaceColor', colors(poly.(p_idx).facies,:))
                    hold on
                    plot(nodes(:,1), nodes(:,2), 'r.')   
                    hold on
                end
            end
        end

        function poly = makeQuadrilaterals(all_polys, poly, poly_num, poly_neighbors, ...
                                            fac_scale, Lx_glob, Lz_glob, Nx_glob, Nz_glob)
            % Is this needed !?
            poly_fault = Faults(all_polys, poly_num, fac_scale, false); 
            p_idx_pebi = strcat('p', string(poly_num), 'f');
            poly_fault.p_idx = p_idx_pebi;
            poly.(p_idx_pebi) = poly_fault;            

            % --- part A: main vertical segment ---
            polyA = Faults(all_polys, poly_num, fac_scale, false);   
            p_idx = strcat('p', string(poly_num), 'fA');
            polyA.p_idx = p_idx; 
            % What about polyA.p
            if poly_num == 12               
                p5B = poly_neighbors{1,1};
                p5C = poly_neighbors{2,1};

                polyA = Faults.pinchOutGrid(polyA, [], p5B, p5C, ...
                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx) = polyA;
                return; % no more parts for p12 -> quit

            elseif poly_num == 25
                p32 = poly_neighbors{1,1};
                p29B = poly_neighbors{2,1};
                p29C = poly_neighbors{3,1};
                p19B = poly_neighbors{4,1};
                p19C = poly_neighbors{5,1};
                p14 = poly_neighbors{6,1};
                p4 = poly_neighbors{7,1};
                p11right = poly_neighbors{8,1};
                p5A = poly_neighbors{9,1};
                p13 = poly_neighbors{10,1};
                p18 = poly_neighbors{11,1};
                p30 = poly_neighbors{12,1};               

                partA_top = findOverlappingNodes(poly_fault, p32, 'top');

                partA_west = [p4.G.nodes.coords(p4.east_mask, :); ...
                              p14.G.nodes.coords(p14.east_mask, :); ...
                              p19C.G.nodes.coords(p19C.east_mask, :); ... 
                              p19B.G.nodes.coords(p19B.east_mask, :); ... 
                              p29C.G.nodes.coords(p29C.east_mask, :); ...                                                                                      
                              p29B.G.nodes.coords(p29B.east_mask, :)];

                p11right_sub = findOverlappingNodes(poly_fault, p11right, 'east');
                partA_east = [p11right_sub; ...
                              p5A.G.nodes.coords(p5A.west_mask, :); ...
                              p13.G.nodes.coords(p13.west_mask, :); ...
                              p18.G.nodes.coords(p18.west_mask, :); ...
                              p30.G.nodes.coords(p30.west_mask, :)];  
                              
                partA_bottom = partA_west(ismembertol(partA_west, partA_east, 'ByRows',true), :);
                
            elseif poly_num == 31
                p32 = poly_neighbors{1,1};
                p28 = poly_neighbors{2,1};
                p20 = poly_neighbors{3,1};
                p15 = poly_neighbors{4,1};
                p6 = poly_neighbors{5,1};
                p17 = poly_neighbors{6,1};
                p11left = poly_neighbors{7,1};
                p4 = poly_neighbors{8,1};
                p14 = poly_neighbors{9,1};
                p19A = poly_neighbors{10,1};
                p29A = poly_neighbors{11,1};
                
                partA_top = findOverlappingNodes(poly_fault, p32, 'top');

                p11left_sub = findOverlappingNodes(poly_fault, p11left, 'west');
                partA_west = [flip(p11left_sub);
                              p17.G.nodes.coords(p17.east_mask, :); ... 
                              p6.G.nodes.coords(p6.east_mask, :); ... 
                              p15.G.nodes.coords(p15.east_mask, :); ...                                                                                      
                              p20.G.nodes.coords(p20.east_mask, :); ...
                              p28.G.nodes.coords(p28.east_mask, :)];
                                              
                partA_east = [p4.G.nodes.coords(p4.west_mask, :); ...
                              p14.G.nodes.coords(p14.west_mask, :); ...
                              p19A.G.nodes.coords(p19A.west_mask, :); ...
                              p29A.G.nodes.coords(p29A.west_mask, :)];  

                partA_bottom = partA_west(ismembertol(partA_west, partA_east, 'ByRows',true), :);                                            
            end

            partA_west = unique(partA_west, 'rows', 'stable'); % remove duplicating overlapping nodes
            partA_east = unique(partA_east, 'rows', 'stable');
            polyA.top_side = partA_top;
            polyA.bottom_side = partA_bottom;
            
            %Nx_sub = size(partA_top, 1);
            Nx_sub = max(size(partA_top,1), size(partA_bottom,1));
            Nz_sub = max(size(partA_west,1), size(partA_east,1));
            
            main_west = size(partA_west,1) > size(partA_east,1);
            main_top = size(partA_top,1) > size(partA_bottom,1);
            if main_west % west side contains more points -> east will be collapsed side (containing most collapsed points)
                cside_x = 'east'; 
            else
                cside_x = 'west';
            end   
            if main_top
                cside_z = 'top';
            else
                cside_z = 'bottom';
            end

            polyA = Faults.collapsedGrid(polyA, cside_x, cside_z, Nx_sub, Nz_sub, ...
                                         partA_west, partA_east, partA_top, partA_bottom, ...
                                         Lx_glob, Lz_glob, Nx_glob, Nz_glob);
        
            if Nx_sub > 2 % collapse internal stacks according to closest side
                Nx_int = Nx_sub - 2;
                % Internals: copy z-coord of modified east or west side
                west_side_z = polyA.G.nodes.coords(polyA.west_mask, 2);
                west_int_z = west_side_z(2:end-1); % omit top and bottom nodes -> these are all fixed
                east_side_z = polyA.G.nodes.coords(polyA.east_mask, 2);
                east_int_z = east_side_z(2:end-1);
                for i=1:ceil(Nx_int/2)                        
                    polyA.G.nodes.coords((i+1)+Nx_sub:Nx_sub:Nx_sub*Nz_sub, 2) = west_int_z; % +Nx_sub to start at second row, not at bottom                     
                    polyA.G.nodes.coords((Nx_sub-i)+Nx_sub:Nx_sub:Nx_sub*Nz_sub, 2) = east_int_z;
                end

                % Internals: uniformly distribute x-coords between west
                % and east side
                for k=2:Nz_sub-1
                    horz_layer = (k-1)*Nx_sub+1:k*Nx_sub; % entire horizontal layer
                    x_sub = polyA.G.nodes.coords(horz_layer, 1);                        
                    x_min = x_sub(1);
                    x_max = x_sub(end);
                    dx = (x_max - x_min)/(Nx_sub-1);
                    polyA.G.nodes.coords(horz_layer, 1) = x_min + cumsum(repmat(dx, Nx_sub, 1)) - dx;
                end       
            end                               

            poly.(p_idx) = polyA;                

            if poly_num == 25 % Only pinch-outs for poly25
                % --- part B: upper pinch-out ---
                polyB = Faults(all_polys, poly_num, 2, false);   
                p_idx = strcat('p', string(poly_num), 'fB');
                polyB.p_idx = p_idx;

                polyB = Faults.pinchOutGrid(polyB, polyA, p29B, p29C, ...
                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx) = polyB;

                % --- part C: lower pinch-out ---
                polyC = Faults(all_polys, poly_num, 2, false);   
                p_idx = strcat('p', string(poly_num), 'fC');
                polyC.p_idx = p_idx;

                polyC = Faults.pinchOutGrid(polyC, polyA, p19B, p19C, ...
                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx) = polyC;
            end            
        end
   

        function [linked_idx, linked_idx_main] = collapseSides(poly, nodes, nodes_main, side, use_buff)
            % Set coordinate of collapsing nodes on given side
            buff = 0;
            if strcmp(side, 'west') || strcmp(side, 'east')
                x_or_z = 2;                
            elseif strcmp(side, 'top') || strcmp(side, 'bottom')
                x_or_z = 1;
                %main_smallest = mean(nodes(:,1)) > mean(nodes_main(:,1));
                %sgn_buff = main_smallest - (~main_smallest);
                if use_buff
                    buff = mean([nodes(1,1)-nodes_main(1,1), ...
                                 nodes(end,1)-nodes_main(end,1)]);
                end
            end
            
            % First: connect all nodes_main with nodes
            linked_idx = [1]; % top and bottom nodes already linked
            linked_idx_main = (1:size(nodes_main,1))';
            %while ix >= 1 && ix <= Nx_sub
            for k=2:size(nodes_main,1)-1
                nm = nodes_main(k,:);
                [~, min_idx] = min(abs((nm(x_or_z)+buff) - nodes(:,x_or_z)));                                
                %poly.G.nodes.coords(ix + (k-1)*Nx_sub, :) = nodes(min_zidx,:);                              
                linked_idx = [linked_idx; min_idx];
            end
            linked_idx = [linked_idx; size(nodes,1)];

            % Second: connect missed nodes to closest node in nodes_main
            diff_idx = diff(linked_idx);
            missed_idx = find(diff_idx > 1);
            missed_idx_adapt = missed_idx;
            for k=1:numel(missed_idx)
                m = missed_idx_adapt(k);
                m_orig = missed_idx(k);
                for i=1:diff_idx(m_orig)-1
                    m_miss = linked_idx(m) + i; % 5, 6
                    [~, min_idx] = min(abs(nodes(m_miss,x_or_z) - nodes_main(:,x_or_z)));
                    m_idx = m + 1;
                    linked_idx(m_idx+1:end+1) = linked_idx(m_idx:end);
                    linked_idx(m_idx) = m_miss;
                    linked_idx_main(m_idx+1:end+1) = linked_idx_main(m_idx:end);
                    linked_idx_main(m_idx) = min_idx;
                    %poly.G.nodes.coords(
                end
                missed_idx_adapt(k+1:end) = missed_idx_adapt(k+1:end) + i; % adjust missed index for added elements
            end                       

        end

        function poly = collapsedGrid(poly, cx, cz, Nx_sub, Nz_sub, ...
                                      part_west, part_east, part_top, part_bottom, ...
                                      Lx_glob, Lz_glob, Nx_glob, Nz_glob, varargin)            
            bottom_tip = true;
            use_buff = false;
            if nargin > 13
                bottom_tip = varargin{1};
            end
            % Determine what nodes should collapse
            if strcmp(cx, 'west')
                nodes_cx = part_west; % collapsed horizontal side
                nodes_mx = part_east; % main horizontal side
            elseif strcmp(cx, 'east')
                nodes_cx = part_east;
                nodes_mx = part_west;
            end
            if strcmp(cz, 'top')
                nodes_cz = part_top; % collapsed vertical side
                nodes_mz = part_bottom; % main vertical side
            elseif strcmp(cz, 'bottom')
                nodes_cz = part_bottom;
                nodes_mz = part_top;
            end
            
            if strcmp(poly.p_idx, 'p24fC') || strcmp(poly.p_idx, 'p23fA')
                use_buff = true;
            end
            
            x_collapse = size(nodes_cx,1) ~= size(nodes_mx,1);
            z_collapse = size(nodes_cz,1) ~= size(nodes_mz,1);

            if x_collapse
                [links_cx, links_mx] = Faults.collapseSides(poly, nodes_cx, nodes_mx, cx, use_buff);
                Nz_new = numel(links_cx); 
            else
                Nz_new = Nz_sub;
            end
            if strcmp(poly.p_idx, 'p25fA') || strcmp(poly.p_idx, 'p25fB') || ...
                    strcmp(poly.p_idx, 'p25fC') || strcmp(poly.p_idx, 'p12fA') || ...
                    strcmp(poly.p_idx, 'p31fA') || strcmp(poly.p_idx, 'p24fA')
                Nx_new = Nx_sub;
            elseif z_collapse
                [links_cz, links_mz] = Faults.collapseSides(poly, nodes_cz, nodes_mz, cz, use_buff);
                Nx_new = numel(links_cz);
            else
                Nx_new = Nx_sub;
            end

            % Create new conform grid with collapsed points added to the dimensions                                             
%            poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, Nx_sub, Nz_new-1);
            poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, Nx_new, Nz_new-1);

            if strcmp(cx, 'west')
                cx_mask = poly.west_mask;            
                mx_mask = poly.east_mask;
            elseif strcmp(cx, 'east')
                cx_mask = poly.east_mask;
                mx_mask = poly.west_mask;
            end
            if strcmp(cz, 'top')
                cz_mask = poly.top_mask;            
                mz_mask = poly.bottom_mask;
            elseif strcmp(cz, 'bottom')
                cz_mask = poly.bottom_mask;
                mz_mask = poly.top_mask;
            end

            % Collapse bottom nodes into tip of fault
            if bottom_tip
                poly.G.nodes.coords(poly.bottom_mask, :) = repmat(part_bottom, nnz(poly.bottom_mask), 1);
                poly.G.nodes.coords(poly.top_mask, :) = part_top;
            elseif strcmp(poly.p_idx, 'p25fA') || strcmp(poly.p_idx, 'p25fB') || ...
                    strcmp(poly.p_idx, 'p25fC') || strcmp(poly.p_idx, 'p12fA') || ...
                    strcmp(poly.p_idx, 'p31fA') || strcmp(poly.p_idx, 'p24fA')
                num_collapse = nnz(poly.bottom_mask) - size(part_bottom, 1);
                dx_right = abs(poly.internal_east(end,1) - poly.internal_east(1,1));
                dx_left = abs(poly.internal_west(end,1) - poly.internal_west(1,1));
                part_bottom_orig = part_bottom;
                if dx_right > dx_left % most irregularity on east side -> collapse on this side                    
                    for i=1:num_collapse % duplicate/collapse bottom nodes
                        part_bottom = [part_bottom; part_bottom_orig(end-i+1,:)];
                    end
                else
                    for i=1:num_collapse
                        part_bottom = [part_bottom; part_bottom_orig(i,:)];
                    end
                end
                poly.G.nodes.coords(poly.bottom_mask, :) = part_bottom;
                poly.G.nodes.coords(poly.top_mask, :) = part_top;
            elseif z_collapse             
                % Fix main vertical side
                nodes_mz_new = zeros(Nx_new,2);
                part_orig_idx = logical([1; diff(links_mz) == 1]); % if zero, we have a duplicate index indicating a collapsed node
                nodes_mz_new(part_orig_idx,:) = nodes_mz;
                counts = histcounts(links_mz)';
                collapsed_idx = find(counts > 1);
                nodes_mz_new(~part_orig_idx,:) = nodes_mz(collapsed_idx, :);
                poly.G.nodes.coords(mz_mask, :) = nodes_mz_new;
                % Fix collapsed horizontal side -> extract collapsed points from links_cx
                nodes_cz_new = nodes_cz(links_cz, :);
                poly.G.nodes.coords(cz_mask, :) = nodes_cz_new;
            else
                poly.G.nodes.coords(poly.bottom_mask, :) = part_bottom;
                poly.G.nodes.coords(poly.top_mask, :) = part_top;
            end
                     
            if x_collapse
                % Fix main horizontal side (NB: remember new points added!)
                nodes_mx_new = zeros(Nz_new,2);
                part_orig_idx = logical([1; diff(links_mx) == 1]); % if zero, we have a duplicate index indicating a collapsed node
                nodes_mx_new(part_orig_idx,:) = nodes_mx;
                counts = histcounts(links_mx)';
                collapsed_idx = find(counts > 1);
                nodes_mx_new(~part_orig_idx,:) = nodes_mx(collapsed_idx, :);
                poly.G.nodes.coords(mx_mask, :) = nodes_mx_new;
                % Fix collapsed horizontal side -> extract collapsed points from links_cx
                nodes_cx_new = nodes_cx(links_cx, :);
                poly.G.nodes.coords(cx_mask, :) = nodes_cx_new;
            else
                poly.G.nodes.coords(poly.west_mask, :) = part_west;
                poly.G.nodes.coords(poly.east_mask, :) = part_east;
            end
        end

        function poly = pinchOutGrid(poly, polyA, upper_N, lower_N, ...
                                     Lx_glob, Lz_glob, Nx_glob, Nz_glob)
            part_top = upper_N.G.nodes.coords(upper_N.bottom_mask, :);      
            part_bottom = lower_N.G.nodes.coords(lower_N.top_mask, :); 
            part_west = part_top(ismembertol(part_top, part_bottom, 'ByRows',true), :);
            if isempty(polyA) % pinch-out at boundary of global grid -> set east nodes manually
                p = poly.p;
                p_east = p(p(:,1) == max(p(:,1)), :);
                z_min = min(p_east(:,2));
                z_max = max(p_east(:,2));
                Lz = z_max - z_min;
                fac = Lz/Lz_glob;
                Nz = ceil(poly.fac_scale*fac*Nz_glob);
                dz = Lz/Nz;
                z = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;
                x = repmat(max(p(:,1)), Nz+1, 1);
                part_east = [x,z];
            else
                polyA_west = polyA.G.nodes.coords(polyA.west_mask, :); % only west side overlaps with the pinch-outs
                part_east = polyA_west(polyA_west(:,2) >= part_bottom(end,2) & ...
                                        polyA_west(:,2) <= part_top(end,2), :); % in case additional nodes have been distributed for polyA at location of pinch-out
            end
            part_east = unique(part_east, 'rows'); % to avoid getting potential duplicates from collapsing nodes
            %part_east = [part_bottom(ismembertol(part_bottom, polyA_west, 'ByRows',true), :); ...
            %              part_top(ismembertol(part_top, polyA_west, 'ByRows',true), :)];

            poly.top_side = part_top;
            poly.bottom_side = part_bottom;
            %Nx_int = 0; % start with no internal points
            %Nx_sub = Nx_int + 2;
            Nx_sub = max(size(part_top, 1), size(part_bottom, 1));
            Nz_sub = max(size(part_west,1), size(part_east,1));

            poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, Nx_sub, Nz_sub-1);
            
            poly.G.nodes.coords(poly.top_mask, :) = part_top;
            poly.G.nodes.coords(poly.bottom_mask, :) = part_bottom;
            poly.G.nodes.coords(poly.west_mask, :) = repmat(part_west, nnz(poly.west_mask), 1);
            poly.G.nodes.coords(poly.east_mask, :) = part_east;
            
            poly = interpolateInternal(poly, poly.top_mask, poly.bottom_mask, []);
        end

        
        function poly = heterogeneousFault(all_polys, poly, fac_scale, Lx_glob, Lz_glob, Nx_glob, Nz_glob)
            % --- Upper piece: poly24 ---
            poly24 = poly.p24f;
            
            % - part A -
            polyA = Faults(all_polys, 24, fac_scale, false);
            polyA.p_idx = 'p24fA';

            bottom_left_A = poly24.bnodes(8,:);
            bottom_right_A = poly24.bnodes(end,:);
            x = [bottom_left_A(1); bottom_right_A(1)];
            z = [bottom_left_A(2); bottom_right_A(2)];
            Lx = x(2) - x(1);
            fac = Lx/Lx_glob;
            Nx = ceil(polyA.fac_scale*fac*Nx_glob);
            %Nx = size(poly24.external_top, 1) - 1; % num CELLS at bottom side
            dx = Lx/Nx;
            x_new = x(1) + cumsum(repmat(dx, Nx+1, 1)) - dx;

            z_new = interp1(x, z, x_new, 'linear'); 

            partA_top = findOverlappingNodes(polyA, poly.p6, 'top');
            partA_bottom = [x_new, z_new];
            p27_east = poly.p27.G.nodes.coords(poly.p27.east_mask, :);
            partA_west = p27_east(p27_east(:,2) >= bottom_left_A(2), :);
            p17_west = poly.p17.G.nodes.coords(poly.p17.west_mask, :);
            partA_east = p17_west(p17_west(:,2) >= bottom_right_A(2), :);

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p24fA = polyA;

            % - part B -
            polyB = Faults(all_polys, 24, fac_scale, false);
            polyB.p_idx = 'p24fB';

            bottom_right_B = poly24.internal_top;
            % Define remaining parts of bottom side
            %Nx_bottom_left = size(poly.p24fA.internal_bottom,1) - size(bottom_right_B,1); % number of CELLS in bottom left part of B
            [~, xmin_idx] = min(abs(poly24.external_west(:,2) - bottom_right_B(1,2)));
            left_pt = poly24.external_west(xmin_idx, :);
            right_pt = bottom_right_B(1,:);
            x_min = left_pt(1);
            x_max = right_pt(1);
            Nx_bottom_left = ceil((x_max - x_min)/Lx_glob*Nx_glob);                        
            dx = (x_max - x_min)/Nx_bottom_left;
            x_new = x_min + cumsum(repmat(dx, Nx_bottom_left+1, 1)) - dx;
            x = [x_min, x_max];
            z = [left_pt(2), right_pt(2)];
            z_new = interp1(x, z, x_new, 'linear');
            bottom_left_B = [x_new(1:end-1), z_new(1:end-1)];

            partB_bottom = [bottom_left_B; bottom_right_B];
            partB_top = partA_bottom;
            partB_west = poly24.external_west(poly24.external_west(:,2) >= left_pt(2) & ... % bottom left of part B
                                                poly24.external_west(:,2) <= partA_bottom(1,2), :); % top left of part B
            partB_west = flip(partB_west); % to get in same order as coordinate masks
            partB_east = poly24.external_east(poly24.external_east(:,2) >= min(poly24.external_east(:,2)) & ...
                                                poly24.external_east(:,2) <= partA_bottom(end,2), :);
            
            polyB.internal_top = partB_top;
            polyB.internal_bottom = partB_bottom;
            polyB.internal_west = partB_west;
            polyB.internal_east = partB_east;

            polyB = Faults.distributeAndCollapseNodes(polyB, partB_west, partB_east, ...
                                                            partB_top, partB_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p24fB = polyB; 

            % - part C -
            polyC = Faults(all_polys, 24, fac_scale, false);
            polyC.p_idx = 'p24fC';
            
            partC_top = partB_bottom(partB_bottom(:,1) <= poly24.internal_top(1,1), :);

            
            %Nx_bottom_C = size(partC_top,1) - 1; % num CELLS at bottom
            %[x_interp, z_interp] = Faults.interpFault(poly24.internal_bottom, Nx_bottom_C);
            %partC_bottom = [x_interp, z_interp];
            partC_bottom = poly24.internal_bottom;

            partC_west = poly24.external_west(poly24.external_west(:,2) >= partC_bottom(1,2) & ...
                                              poly24.external_west(:,2) <= partC_top(1,2), :);
            partC_west = flip(partC_west);
            partC_east = poly24.internal_east;

            polyC.internal_top = partC_top;
            polyC.internal_bottom = partC_bottom;
            polyC.internal_west = partC_west;
            polyC.internal_east = partC_east;

            polyC = Faults.distributeAndCollapseNodes(polyC, partC_west, partC_east, ...
                                                            partC_top, partC_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p24fC = polyC;

            % --- Small piece: poly9 ---
            poly9 = poly.p9f;
            
            % - part A -
            polyA = Faults(all_polys, 9, fac_scale, false);
            polyA.p_idx = 'p9fA';

            partA_top = poly9.internal_top;
            %Nx_bottom_A = size(partA_top, 1) - 1; % num CELLS, not nodes
            %[x_interp, z_interp] = Faults.interpFault(poly9.internal_bottom, Nx_bottom_A);
            %partA_bottom = [x_interp, z_interp];
            partA_bottom = poly9.internal_bottom;
            partA_west = poly24.internal_east(poly24.internal_east(:,2) >= partA_bottom(1,2), :);
            partA_east = poly9.external_east;

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p9fA = polyA;

            % --- Irregular piece: poly21 ---
            poly21 = poly.p21f;
            
            % - part A -
            polyA = Faults(all_polys, 21, fac_scale, false);
            polyA.p_idx = 'p21fA';

            partA_top = poly9.internal_bottom;
            
            bottom_right_A = poly.p7small.G.nodes.coords(poly.p7small.top_mask, :);
            % Define remaining parts of bottom side            
            [~, xmin_idx] = min(abs(poly21.internal_west(:,2) - bottom_right_A(1,2)));
            left_pt = poly21.internal_west(xmin_idx, :);
            right_pt = bottom_right_A(1,:);
            x_min = left_pt(1);
            x_max = right_pt(1);
            Nx_bottom_left = ceil((x_max - x_min)/Lx_glob*Nx_glob);                        
            dx = (x_max - x_min)/Nx_bottom_left;
            x_new = x_min + cumsum(repmat(dx, Nx_bottom_left+1, 1)) - dx;
            x = [x_min, x_max];
            z = [left_pt(2), right_pt(2)];
            z_new = interp1(x, z, x_new, 'linear');
            bottom_left_A = [x_new(1:end-1), z_new(1:end-1)];

            partA_bottom = [bottom_left_A; bottom_right_A];
            partA_west = poly21.internal_west(poly21.internal_west(:,2) >= bottom_left_A(1,2), :);
            partA_east = poly21.external_east(poly21.external_east(:,2) >= bottom_right_A(end,2), :);

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p21fA = polyA;

            % - part B -
            polyB = Faults(all_polys, 21, fac_scale, false);
            polyB.p_idx = 'p21fB';

            bottom_left_B = poly24.internal_bottom(end,:);
            p7small_west = poly.p7small.G.nodes.coords(poly.p7small.west_mask, :);
            [~, min_idx] = min(abs(p7small_west(:,2) - bottom_left_B(2)));
            bottom_right_B = p7small_west(min_idx, :);
            Lx_bottom_B = bottom_right_B(1) - bottom_left_B(1);
            Nx_bottom_B = ceil(Lx_bottom_B/Lx_glob * Nx_glob);
            bottom_B = [bottom_left_B; bottom_right_B];
            [x_interp, z_interp] = Faults.interpFault(bottom_B, Nx_bottom_B);           
            partB_bottom = [x_interp, z_interp];

            partB_top = polyA.G.nodes.coords(polyA.bottom_mask, :);
            partB_top = partB_top(partB_top(:,1) <= bottom_right_A(1,1), :);

            partB_west = poly24.internal_east(poly24.internal_east(:,2) <= bottom_left_A(1,2), :);
            partB_east = p7small_west(p7small_west(:,2) >= bottom_right_B(2), :);

            polyB.internal_top = partB_top;
            polyB.internal_bottom = partB_bottom;
            polyB.internal_west = partB_west;
            polyB.internal_east = partB_east;

            polyB = Faults.distributeAndCollapseNodes(polyB, partB_west, partB_east, ...
                                                            partB_top, partB_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p21fB = polyB;

            % - part C -
            polyC = Faults(all_polys, 21, fac_scale, false);
            polyC.p_idx = 'p21fC';

            partC_top = [poly24.internal_bottom; partB_bottom(2:end, :)]; % 2:end to avoid duplicating node at overlap
            %partC_west = findOverlappingNodes(polyC, poly.p8, 'west');
            partC_west = flip(poly21.external_west);

            p7sG_bottom = poly.p7small.G.nodes.coords(poly.p7small.bottom_mask, :);
            bottom_right_C = p7sG_bottom(1,:);
            dist = sqrt((bottom_right_C(1) - partC_west(:,1)).^2 + ...
                        (bottom_right_C(2) - partC_west(:,2)).^2);
            [~, min_idx] = min(dist);
            bottom_left_C = partC_west(min_idx, :);
            Lx_bottom_C = bottom_right_C(1) - bottom_left_C(1);
            Nx_bottom_C = ceil(Lx_bottom_C/Lx_glob * Nx_glob);
            bottom_C = [bottom_left_C; bottom_right_C];
            [x_interp, z_interp] = Faults.interpFault(bottom_C, Nx_bottom_C);           
            partC_bottom = [x_interp, z_interp];
            partC_west = partC_west(partC_west(:,2) >= partC_bottom(1,2), :);

            partC_east = p7small_west(p7small_west(:,2) <= partC_top(end,2), :);

            polyC.internal_top = partC_top;
            polyC.internal_bottom = partC_bottom;
            polyC.internal_west = partC_west;
            polyC.internal_east = partC_east;

            polyC = Faults.distributeAndCollapseNodes(polyC, partC_west, partC_east, ...
                                                            partC_top, partC_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p21fC = polyC;

            % - part D -
            polyD = Faults(all_polys, 21, fac_scale, false);
            polyD.p_idx = 'p21fD';

            partD_top = [partC_bottom; p7sG_bottom(2:end,:)]; % 2:end to avoid duplicate nodes
            partD_bottom = poly21.internal_bottom;
            partD_west = poly21.external_west(poly21.external_west(:,2) <= partC_bottom(1,2), :);
            partD_west = flip(partD_west);
            partD_east = poly21.external_east(poly21.external_east(:,2) <= p7sG_bottom(end,2), :);

            polyD.internal_top = partD_top;
            polyD.internal_bottom = partD_bottom;
            polyD.internal_west = partD_west;
            polyD.internal_east = partD_east;

            polyD = Faults.distributeAndCollapseNodes(polyD, partD_west, partD_east, ...
                                                            partD_top, partD_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);            
            
            poly.p21fD = polyD;
               
            % --- Bottom piece: poly23 ---
            poly23 = poly.p23f;
            
            % - part A -
            polyA = Faults(all_polys, 23, fac_scale, false);
            polyA.p_idx = 'p23fA';

            partA_top = [poly.p26f.internal_bottom; poly.p22f.internal_bottom(2:end,:)];
            p2left_west = poly.p2left.G.nodes.coords(poly.p2left.west_mask, :);
            bottom_right_A = p2left_west(1,:);
            p3_east = poly.p3.G.nodes.coords(poly.p3.east_mask, :);
            dist = sqrt((bottom_right_A(1) - p3_east(:,1)).^2 + ...
                        (bottom_right_A(2) - p3_east(:,2)).^2);
            [~, min_idx] = min(dist);
            bottom_left_A = p3_east(min_idx, :);
            Lx_bottom_A = bottom_right_A(1) - bottom_left_A(1);
            Nx_bottom_A = ceil(Lx_bottom_A/Lx_glob * Nx_glob);
            bottom_A = [bottom_left_A; bottom_right_A];
            [x_interp, z_interp] = Faults.interpFault(bottom_A, Nx_bottom_A);           
            partA_bottom = [x_interp, z_interp];

            partA_west = p3_east(p3_east(:,2) <= partA_top(1,2) & ...
                                 p3_east(:,2) >= bottom_left_A(2), :);
            partA_east = p2left_west(p2left_west(:,2) <= partA_top(end,2), :);
            
            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob, ...
                                                            false, true);

            poly.p23fA = polyA;

            % - part B -
            polyB = Faults(all_polys, 23, fac_scale, false);
            polyB.p_idx = 'p23fB';
            % NB: since part B is to be triangulated, points must be
            % ordered in counterclockwise direction
            partB_top = flip(partA_bottom);
            partB_west = p3_east(p3_east(:,2) <= bottom_left_A(2), :);
            partB_west = flip(partB_west);
            partB_east = poly.p1mid.G.nodes.coords(poly.p1mid.top_mask, :);

            polyB.internal_top = partB_top;
            polyB.internal_west = partB_west;
            polyB.internal_east = partB_east;

            polyB.bnodes = [partB_top; partB_west; partB_east];
            polyB.bnodes = unique(polyB.bnodes, 'rows', 'stable');

            poly.p23fB = polyB;
        end

        function poly = distributeAndCollapseNodes(poly, part_west, part_east, part_top, part_bottom, ...
                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob, bottom_tip, varargin)
            part_west = unique(part_west, 'rows', 'stable'); % remove duplicating overlapping nodes
            part_east = unique(part_east, 'rows', 'stable');
            
            poly.top_side = part_top;
            poly.bottom_side = part_bottom;

            %Nx_sub = size(part_top, 1);
            Nx_sub = max(size(part_top,1), size(part_bottom,1));
            Nz_sub = max(size(part_west,1), size(part_east,1));
            
            main_west = size(part_west,1) > size(part_east,1);
            main_top = size(part_top,1) > size(part_bottom,1);
            if main_west % west side contains more points -> east will be collapsed side (containing most collapsed points)
                cside_x = 'east'; 
            else
                cside_x = 'west';
            end
            if main_top
                cside_z = 'bottom'; 
            else
                cside_z = 'top';
            end

            poly = Faults.collapsedGrid(poly, cside_x, cside_z, Nx_sub, Nz_sub, ...
                                         part_west, part_east, part_top, part_bottom, ...
                                         Lx_glob, Lz_glob, Nx_glob, Nz_glob, bottom_tip);
        
            if Nx_sub > 2 % collapse internal stacks according to closest side
                Nz_sub = poly.G.cartDims(2)+1; % vertical dim may have increased if collapsed nodes added
                Nx_sub = poly.G.cartDims(1)+1; % horizontal dim may also have changed
                Nx_int = Nx_sub - 2;
                interp = false;                
                if nargin > 10
                    interp = varargin{1};
                end
                           
                   
                if interp
                    G_dum = poly.G;
                    pt = unique(poly.G.nodes.coords(poly.west_mask, :), 'rows', 'stable');
                    G_dum.nodes.coords(poly.west_mask, :) = interparc(Nz_sub, pt(:,1), pt(:,2));
                    pt = unique(poly.G.nodes.coords(poly.east_mask, :), 'rows', 'stable');
                    G_dum.nodes.coords(poly.east_mask, :) = interparc(Nz_sub, pt(:,1), pt(:,2));
                    pt = unique(poly.G.nodes.coords(poly.top_mask, :), 'rows', 'stable');
                    G_dum.nodes.coords(poly.top_mask, :) = interparc(Nx_sub, pt(:,1), pt(:,2));
                    pt = unique(poly.G.nodes.coords(poly.bottom_mask, :), 'rows', 'stable');
                    G_dum.nodes.coords(poly.bottom_mask, :) = interparc(Nx_sub, pt(:,1), pt(:,2));

                    for k=2:Nz_sub-1
                        horz_layer = (k-1)*Nx_sub+1:k*Nx_sub; % entire horizontal layer
                        subset = G_dum.nodes.coords(horz_layer, :);  
                        x_min = subset(1,1); x_max = subset(end,1);
                        z_min = subset(1,2); z_max = subset(end,2);
                        x = [x_min; x_max];
                        z = [z_min; z_max];
                        dx = (x_max - x_min)/(Nx_sub-1);
                        x_interp = x_min + cumsum(repmat(dx, Nx_sub, 1)) - dx;
                        z_interp = interp1(round(x,10),z,round(x_interp,10),'linear'); % to avoid weird round-off error
                        %poly.G.nodes.coords(horz_layer(1:end-1), :) = [x_interp(1:end-1), z_interp(1:end-1)];
                        poly.G.nodes.coords(horz_layer(2:end-1), 1) = x_interp(2:end-1); % ONLY SET X-COORD
                    end

                    for i=2:Nx_sub-1
                        vert_layer = i:Nx_sub:(Nz_sub-1)*Nx_sub+i; % entire vertical layer
                        subset = G_dum.nodes.coords(vert_layer, :);  
                        x_min = subset(1,1); x_max = subset(end,1);
                        z_min = subset(1,2); z_max = subset(end,2);
                        x = [x_min; x_max];
                        z = [z_min; z_max];
                        dz = (z_max - z_min)/(Nz_sub-1);
                        %z_interp = poly.G.nodes.coords(vert_layer, 2);
                        z_interp = z_min + cumsum(repmat(dz, Nz_sub, 1)) - dz;
                        x_interp = interp1(round(z,10),x,round(z_interp,10),'linear'); % to avoid weird round-off error
                        %poly.G.nodes.coords(vert_layer, 1) = x_interp;
                        %poly.G.nodes.coords(vert_layer(1:end-1), :) = [x_interp(1:end-1), z_interp(1:end-1)];
                        poly.G.nodes.coords(vert_layer(2:end-1), 2) = z_interp(2:end-1);
                    end    
                else
                    % Internals: copy z-coord of modified east or west side
                    west_side_z = poly.G.nodes.coords(poly.west_mask, 2);
                    if numel(unique(west_side_z)) == 2 % internal node on west side is collapsed to top or bottom node -> set z-coord to midpoint between top and bottom
                        west_int_z = (west_side_z(1)+west_side_z(end))/2;
                    else
                        west_int_z = west_side_z(2:end-1); % omit top and bottom nodes -> these are all fixed
                    end                
                    east_side_z = poly.G.nodes.coords(poly.east_mask, 2);
                    if numel(unique(east_side_z)) == 2 %
                        %east_int_z = (east_side_z(1)+east_side_z(end))/2;
                        east_int_z = west_int_z;
                    else
                        east_int_z = east_side_z(2:end-1);
                    end
                    
                    for i=1:ceil(Nx_int/2)                        
                        poly.G.nodes.coords((i+1)+Nx_sub:Nx_sub:Nx_sub*(Nz_sub-1), 2) = west_int_z; % +Nx_sub to start at second row, not at bottom                     
                        poly.G.nodes.coords((Nx_sub-i)+Nx_sub:Nx_sub:Nx_sub*(Nz_sub-1), 2) = east_int_z;
                    end
                  
    
                    % Internals: uniformly distribute x-coords between west
                    % and east side
                    for k=2:Nz_sub-1
                        horz_layer = (k-1)*Nx_sub+1:k*Nx_sub; % entire horizontal layer
                        x_sub = poly.G.nodes.coords(horz_layer, 1);                        
                        x_min = x_sub(1);
                        x_max = x_sub(end);
    %                     if strcmp(poly.p_idx, 'p24fA')
    %                         x_max = (min(poly.G.nodes.coords(poly.east_mask,1)) + max(poly.G.nodes.coords(poly.east_mask,1)))/2;
    %                     end
                        dx = (x_max - x_min)/(Nx_sub-1);
                        horz_layer_orig = poly.G.nodes.coords(horz_layer, 1);
                        poly.G.nodes.coords(horz_layer, 1) = x_min + cumsum(repmat(dx, Nx_sub, 1)) - dx;
    %                     if strcmp(poly.p_idx, 'p24fA') % OR INTERPOLATE EACH STACK (EXCEPT COLLAPSED) BETWEEN TOP AND BOTTOM?
    %                         poly.G.nodes.coords(k*Nx_sub, 1) = horz_layer_orig(end); % retain original collapsed point
    %                     end
                    end  
                end
            end
                        
        end

        function [x_interp, z_interp] = interpFault(interp_side, Nx_side)
            left_pt = interp_side(1,:);
            right_pt = interp_side(end,:);
            x = [left_pt(1), right_pt(1)];
            z = [left_pt(2), right_pt(2)];
            dx = (x(2) - x(1))/Nx_side;
            x_interp = x(1) + cumsum(repmat(dx, Nx_side+1, 1)) - dx;
            z_interp = interp1(round(x,10), z, round(x_interp,10), 'linear');
        end

   end
end