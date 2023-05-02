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
          
            poly_fault = Faults(all_polys, poly_num, fac_scale, false); 
            p_idx_pebi = strcat('p', string(poly_num), 'f');
            poly_fault.p_idx = p_idx_pebi;
            poly.(p_idx_pebi) = poly_fault;            
                        

            if poly_num == 25 % Only pinch-outs for poly25
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
                
                % --- part F: upper pinch-out ---
                % Quadrilateral part
                polyFq = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_q = strcat('p', string(poly_num), 'fFq');
                polyFq.p_idx = p_idx_q;             
                polyFq.fac_scale = 1.5; % to not get too dense points

                % Triangular part
                polyFt = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_t = strcat('p', string(poly_num), 'fFt');
                polyFt.p_idx = p_idx_t;

                % determine transition to triangle grid
                x_tri_F = 1.31;

                [polyFq, polyFt] = Faults.pinchOutGrid(polyFq, polyFt, p29B, p29C, x_tri_F, ...
                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx_q) = polyFq;
                poly.(p_idx_t) = polyFt;

                % --- part G: lower pinch-out ---
                % Quadrilateral part
                polyGq = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_q = strcat('p', string(poly_num), 'fGq');
                polyGq.p_idx = p_idx_q;
                polyGq.fac_scale = 1.5; % - !!! -

                % Triangular part
                polyGt = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_t = strcat('p', string(poly_num), 'fGt');
                polyGt.p_idx = p_idx_t;

                x_tri_G = 1.235;
                [polyGq, polyGt] = Faults.pinchOutGrid(polyGq, polyGt, p19B, p19C, x_tri_G, ...
                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx_q) = polyGq;
                poly.(p_idx_t) = polyGt;

                % After pinch-outs are defined, continue with vertical
                % segment.
                polyA = Faults(all_polys, poly_num, fac_scale, false);   
                p_idxA = strcat('p', string(poly_num), 'fA');
                polyA.p_idx = p_idxA; 

                partA_top = findOverlappingNodes(poly_fault, p32, 'top');

                partA_west = [polyFq.internal_east; ...
                                p29B.G.nodes.coords(p29B.east_mask, :)];
                %partA_west = unique(partA_west, 'rows', 'stable');

                bottom_left_A = polyFq.internal_east(1,:);
                p30_bottom = p30.G.nodes.coords(p30.bottom_mask, :);
                bottom_right_A = p30_bottom(1,:);
                Lx_bottom_A = bottom_right_A(1) - bottom_left_A(1);
                Nx_bottom_A = ceil(Lx_bottom_A/Lx_glob * Nx_glob);
                bottom_A = [bottom_left_A; bottom_right_A];
                [x_interp, z_interp] = Faults.interpFault(bottom_A, Nx_bottom_A);           
                partA_bottom = [x_interp, z_interp];

                partA_east = p30.G.nodes.coords(p30.west_mask, :);

                polyA.internal_top = partA_top;
                polyA.internal_bottom = partA_bottom;
                polyA.internal_west = partA_west;
                polyA.internal_east = partA_east;
    
                polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                                partA_top, partA_bottom, ...
                                                                Lx_glob, Lz_glob, Nx_glob, Nz_glob);
                                                                
                               
                poly.(p_idxA) = polyA;

                % - part B - 
                polyB = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx = strcat('p', string(poly_num), 'fB');
                polyB.p_idx = p_idx; 

                partB_top = partA_bottom;
                partB_west = [polyGq.internal_east; ...
                              p19B.G.nodes.coords(p19B.east_mask, :); ...
                              p29C.G.nodes.coords(p29C.east_mask, :)];
                partB_west = unique(partB_west, 'rows', 'stable');
                
                bottom_left_B = polyGq.internal_east(1,:);
                p18_west = p18.G.nodes.coords(p18.west_mask, :);
                dist = Faults.euclideanDist(bottom_left_B, p18_west);
                [~, min_idx] = min(dist);
                bottom_right_B = p18_west(min_idx,:);
                Lx_bottom_B = bottom_right_B(1) - bottom_left_B(1);
                Nx_bottom_B = ceil(Lx_bottom_B/Lx_glob * Nx_glob);
                bottom_B = [bottom_left_B; bottom_right_B];
                [x_interp, z_interp] = Faults.interpFault(bottom_B, Nx_bottom_B);           
                partB_bottom = [x_interp, z_interp];

                partB_east = p18_west(p18_west(:,2) >= bottom_right_B(2)-eps, :);

                polyB.internal_top = partB_top;
                polyB.internal_bottom = partB_bottom;
                polyB.internal_west = partB_west;
                polyB.internal_east = partB_east;
    
                polyB = Faults.distributeAndCollapseNodes(polyB, partB_west, partB_east, ...
                                                                partB_top, partB_bottom, ...
                                                                Lx_glob, Lz_glob, Nx_glob, Nz_glob);
                poly.(p_idx) = polyB; 

                % - part C -
                polyC = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx = strcat('p', string(poly_num), 'fC');
                polyC.p_idx = p_idx; 

                partC_top = partB_bottom;
                p14_east = p14.G.nodes.coords(p14.east_mask, :);
                partC_west = [p14_east; ...
                              p19C.G.nodes.coords(p19C.east_mask, :)];

                bottom_left_C = p14_east(1,:);
                p5_west = p5A.G.nodes.coords(p5A.west_mask, :);
                dist = Faults.euclideanDist(bottom_left_C, p5_west);
                [~, min_idx] = min(dist);
                bottom_right_C = p5_west(min_idx,:);
                Lx_bottom_C = bottom_right_C(1) - bottom_left_C(1);
                Nx_bottom_C = ceil(Lx_bottom_C/Lx_glob * Nx_glob);
                bottom_C = [bottom_left_C; bottom_right_C];
                [x_interp, z_interp] = Faults.interpFault(bottom_C, Nx_bottom_C);           
                partC_bottom = [x_interp, z_interp];

                partC_east = [p5_west(p5_west(:,2) >= bottom_right_C(2)-eps, :); ...
                              p13.G.nodes.coords(p13.west_mask, :); ...
                              p18_west(p18_west(:,2) <= bottom_right_B(2)+eps, :)];
                
                polyC.internal_top = partC_top;
                polyC.internal_bottom = partC_bottom;
                polyC.internal_west = partC_west;
                polyC.internal_east = partC_east;
    
                polyC = Faults.distributeAndCollapseNodes(polyC, partC_west, partC_east, ...
                                                                partC_top, partC_bottom, ...
                                                                Lx_glob, Lz_glob, Nx_glob, Nz_glob);
                poly.(p_idx) = polyC;

                % - part D -
                polyD = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx = strcat('p', string(poly_num), 'fD');
                polyD.p_idx = p_idx; 

                partD_top = partC_bottom;

                bottom_right_D = p5_west(1,:);
                p4_east = p4.G.nodes.coords(p4.east_mask, :);
                dist = Faults.euclideanDist(bottom_right_D, p4_east);
                [~, min_idx] = min(dist);
                bottom_left_D = p4_east(min_idx,:);
                Lx_bottom_D = bottom_right_D(1) - bottom_left_D(1);
                Nx_bottom_D = ceil(Lx_bottom_D/Lx_glob * Nx_glob);
                bottom_D = [bottom_left_D; bottom_right_D];
                [x_interp, z_interp] = Faults.interpFault(bottom_D, Nx_bottom_D);           
                partD_bottom = [x_interp, z_interp];

                partD_west = p4_east(p4_east(:,2) >= bottom_left_D(2)-eps, :);
                partD_east = p5_west(p5_west(:,2) <= bottom_right_C(2)+eps, :);

                polyD.internal_top = partD_top;
                polyD.internal_bottom = partD_bottom;
                polyD.internal_west = partD_west;
                polyD.internal_east = partD_east;
    
                polyD = Faults.distributeAndCollapseNodes(polyD, partD_west, partD_east, ...
                                                                partD_top, partD_bottom, ...
                                                                Lx_glob, Lz_glob, Nx_glob, Nz_glob);
                
                poly.(p_idx) = polyD;
                
                % - part E -
                polyE = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx = strcat('p', string(poly_num), 'fE');
                polyE.p_idx = p_idx; 

                % Order counterclockwise
                polyE.internal_top = flip(partD_bottom);
                polyE.internal_west = p4_east(p4_east(:,2) <= bottom_left_D(2), :);
                polyE.internal_west = flip(polyE.internal_west);
                p11right_west = p11right.G.nodes.coords(p11right.west_mask, :);
                polyE.internal_east = p11right_west(p11right_west(:,2) >= polyE.internal_west(end,2), :);

                polyE.bnodes = [polyE.internal_top; polyE.internal_west; polyE.internal_east];
                [~, unique_idx] = uniquetol(polyE.bnodes, 'ByRows',true);
                polyE.bnodes = polyE.bnodes(sort(unique_idx), :);

                poly.(p_idx) = polyE;
            
            elseif poly_num == 12     
                % Quadrilateral part
                polyAq = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_q = strcat('p', string(poly_num), 'fAq');
                polyAq.p_idx = p_idx_q;
                % Triangular part
                polyAt = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_t = strcat('p', string(poly_num), 'fAt');
                polyAt.p_idx = p_idx_t;

                p5B = poly_neighbors{1,1};
                p5C = poly_neighbors{2,1};

                x_tri = 2.30;

                [polyAq, polyAt] = Faults.pinchOutGrid(polyAq, polyAt, p5B, p5C, x_tri, ...
                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob);

                poly.(p_idx_q) = polyAq;
                poly.(p_idx_t) = polyAt;

            elseif poly_num == 31
                polyAq = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_q = strcat('p', string(poly_num), 'fAq');
                polyAq.p_idx = p_idx_q; 

                polyAt = Faults(all_polys, poly_num, fac_scale, false);   
                p_idx_t = strcat('p', string(poly_num), 'fAt');
                polyAt.p_idx = p_idx_t; 

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

                z_tri = 0.69; % hard-coded location !

                p11left_sub = findOverlappingNodes(poly_fault, p11left, 'west');
                p11left_q = p11left_sub(p11left_sub(:,2) >= z_tri, :);                
                p11left_q = flip(p11left_q);
                partA_west = [p11left_q;
                              p17.G.nodes.coords(p17.east_mask, :); ... 
                              p6.G.nodes.coords(p6.east_mask, :); ... 
                              p15.G.nodes.coords(p15.east_mask, :); ...                                                                                      
                              p20.G.nodes.coords(p20.east_mask, :); ...
                              p28.G.nodes.coords(p28.east_mask, :)];
                                              
                p4_west = p4.G.nodes.coords(p4.west_mask, :);
                p4_q = p4_west(p4_west(:,2) >= z_tri, :);
                partA_east = [p4_q; ...
                              p14.G.nodes.coords(p14.west_mask, :); ...
                              p19A.G.nodes.coords(p19A.west_mask, :); ...
                              p29A.G.nodes.coords(p29A.west_mask, :)];
                
                bottom_left_A = partA_west(1,:);
                bottom_right_A = partA_east(1,:);
                Lx_bottom_A = bottom_right_A(1) - bottom_left_A(1);
                Nx_bottom_A = ceil(Lx_bottom_A/Lx_glob * Nx_glob);
                bottom_A = [bottom_left_A; bottom_right_A];
                [x_interp, z_interp] = Faults.interpFault(bottom_A, Nx_bottom_A);           
                partA_bottom = [x_interp, z_interp];


                polyAq.internal_top = partA_top;
                polyAq.internal_bottom = partA_bottom;
                polyAq.internal_west = partA_west;
                polyAq.internal_east = partA_east;
    
                polyAq = Faults.distributeAndCollapseNodes(polyAq, partA_west, partA_east, ...
                                                                partA_top, partA_bottom, ...
                                                                Lx_glob, Lz_glob, Nx_glob, Nz_glob);
                               
                poly.(p_idx_q) = polyAq;

                % Triangular bottom part
                partA_top = flip(polyAq.internal_bottom);
                partA_west = p11left_sub(p11left_sub(:,2) < z_tri, :); 
                %partA_west = flip(partA_west);
                partA_east = p4_west(p4_west(:,2) < z_tri, :);

                polyAt.bnodes = [partA_top; partA_west; partA_east];
                [~, unique_idx] = uniquetol(polyAt.bnodes, 'ByRows', true); % stable to give in same order as list of neighbors
                polyAt.bnodes = polyAt.bnodes(sort(unique_idx), :);

                poly.(p_idx_t) = polyAt;

            end
                          
                     
        end
   

        function [linked_idx, linked_idx_main] = collapseSides(poly, nodes, nodes_main, side, use_buff)
            % Set coordinate of collapsing nodes on given side
            buff = 0;
            if strcmp(side, 'west') || strcmp(side, 'east')
                x_or_z = 2;                
            elseif strcmp(side, 'top') || strcmp(side, 'bottom')
                x_or_z = 1;
                if use_buff % used for very skewed subgrids
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
            use_buff = false;            
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
            
            if strcmp(poly.p_idx, 'p24fC') || strcmp(poly.p_idx, 'p23fA') ...
                    || strcmp(poly.p_idx, 'p25fA') || strcmp(poly.p_idx, 'p25fB') ...
                    || strcmp(poly.p_idx, 'p25fC') || strcmp(poly.p_idx, 'p25fD') ...
                    || strcmp(poly.p_idx, 'p31fAq')
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
           
            if z_collapse
                [links_cz, links_mz] = Faults.collapseSides(poly, nodes_cz, nodes_mz, cz, use_buff);
                Nx_new = numel(links_cz);
            else
                Nx_new = Nx_sub;
            end

            % Create new conform grid with collapsed points added to the dimensions 
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
                       
            if z_collapse             
                % Fix main vertical side
                nodes_mz_new = zeros(Nx_new,2);
                part_orig_idx = logical([1; diff(links_mz) == 1]); % if zero, we have a duplicate index indicating a collapsed node
                nodes_mz_new(part_orig_idx,:) = nodes_mz;
                ulinks_mz = unique(links_mz);                
                counts = histc(links_mz, ulinks_mz)'; % use histc, I don't understand last element of histcounts
                %counts = histcounts(links_mz)';
                collapsed_idx = find(counts > 1);
                nodes_mz_new(~part_orig_idx,:) = nodes_mz(collapsed_idx, :);
                poly.G.nodes.coords(mz_mask, :) = nodes_mz_new;
                % Fix collapsed horizontal side -> extract collapsed points
                % from links_cz
                nodes_cz_new = nodes_cz(links_cz, :);
                poly.G.nodes.coords(cz_mask, :) = nodes_cz_new;
            else
                poly.G.nodes.coords(poly.bottom_mask, :) = part_bottom;
                poly.G.nodes.coords(poly.top_mask, :) = part_top;
            end
                     
            if x_collapse
                % Fix main horizontal side (NB: new points may have been added!)
                nodes_mx_new = zeros(Nz_new,2);
                part_orig_idx = logical([1; diff(links_mx) == 1]); % if zero, we have a duplicate index indicating a collapsed node
                nodes_mx_new(part_orig_idx,:) = nodes_mx;
                ulinks_mx = unique(links_mx);                
                counts = histc(links_mx, ulinks_mx)';                
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

        function [poly_quad, poly_tri] = pinchOutGrid(poly_quad, poly_tri, upper_N, lower_N, x_tri, ...
                                     Lx_glob, Lz_glob, Nx_glob, Nz_glob)
            % Make quadrilateral grid
            upper_bottom = upper_N.G.nodes.coords(upper_N.bottom_mask, :);                
            lower_top = lower_N.G.nodes.coords(lower_N.top_mask, :);
            % find line segment with x-coords closest to x_tri
            [~, min_idx] = min(abs(upper_bottom(:,1)-x_tri) + abs(lower_top(:,1)-x_tri));
            part_bottom = lower_top(min_idx:end, :);
            part_top = upper_bottom(min_idx:end, :);

            upper_east = upper_N.G.nodes.coords(upper_N.east_mask, :);
            lower_east = lower_N.G.nodes.coords(lower_N.east_mask, :);
            point_top_east = part_top(ismembertol(part_top, upper_east, 'ByRows',true), :);
            point_bottom_east = part_bottom(ismembertol(part_bottom, lower_east, 'ByRows',true), :);
           
            point_top_west = part_top(1,:);
            point_bottom_west = part_bottom(1,:);

            p_east = [point_bottom_east; point_top_east];
            x = p_east(:,1);
            z = p_east(:,2);
            z_min = min(z); z_max = max(z);
            Lz = z_max - z_min;
            fac = Lz/Lz_glob;
            Nz = ceil(poly_quad.fac_scale*fac*Nz_glob);
            dz = Lz/Nz;
            z_interp = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;     
            x_interp = interp1(round(z,10),x,round(z_interp,10),'linear');
            part_east = [x_interp, z_interp];

            p_west = [point_bottom_west; point_top_west];
            x = p_west(:,1);
            z = p_west(:,2);
            z_min = min(z); z_max = max(z);
            Lz = z_max - z_min;
            dz = Lz/Nz; % use dimension of east part
            z_interp = z_min + cumsum(repmat(dz, Nz+1, 1)) - dz;        
            x_interp = interp1(round(z,10),x,round(z_interp,10),'linear');
            part_west = [x_interp, z_interp];
              
            poly_quad.external_top = part_top;
            poly_quad.external_west = part_west;
            poly_quad.external_bottom = part_bottom;
            poly_quad.external_east = part_east;
            poly_quad.internal_east = poly_quad.external_east;

            poly_quad.top_side = part_top;
            poly_quad.bottom_side = part_bottom; 

            Nx_sub = max(size(part_top, 1), size(part_bottom, 1));
            Nz_sub = max(size(part_west,1), size(part_east,1));

            poly_quad = cartesianSubgrid(poly_quad, Lx_glob, Lz_glob, Nx_glob, Nz_glob, Nx_sub, Nz_sub-1);
            
            poly_quad.G.nodes.coords(poly_quad.top_mask, :) = part_top;
            poly_quad.G.nodes.coords(poly_quad.bottom_mask, :) = part_bottom;
            poly_quad.G.nodes.coords(poly_quad.west_mask, :) = part_west;
            poly_quad.G.nodes.coords(poly_quad.east_mask, :) = part_east;
            
            poly_quad = interpolateInternal(poly_quad, poly_quad.top_mask, poly_quad.bottom_mask, []);
            poly_quad.G.i = nan(poly_quad.G.cells.num, 1); % unstructured -> set to nan
            poly_quad.G.j = nan(poly_quad.G.cells.num, 1);

            % Make triangular grid
            part_top = upper_bottom(upper_bottom(:,1) < x_tri, :);
            part_top = flip(part_top);
            part_bottom = lower_top(lower_top(:,1) < x_tri, :);
            part_east = poly_quad.external_west;

            poly_tri.external_top = part_top;
            poly_tri.external_bottom = part_bottom;
            poly_tri.internal_east = part_east;

            poly_tri.bnodes = [part_top; part_bottom; part_east];
            poly_tri.bnodes = unique(poly_tri.bnodes, 'rows', 'stable');
        end

        
        function poly = heterogeneousFault(all_polys, poly, fac_scale, Lx_glob, Lz_glob, Nx_glob, Nz_glob)
            % Make composite quadrilateral-triangulated grid of bottom
            % heteroegenous fault. 
            % Yep, it gets pretty nasty...
            
            % --- Upper piece: poly24 ---
            poly24 = poly.p24f;
            
            % - part A -
            polyA = Faults(all_polys, 24, fac_scale, false);
            polyA.p_idx = 'p24fA';

            left_target =  [0.487, 0.720];
            right_target = [0.508, 0.727];

            dist_left = Faults.euclideanDist(left_target, poly24.bnodes);
            dist_right = Faults.euclideanDist(right_target, poly24.bnodes);
            [~, min_left] = min(dist_left);
            [~, min_right] = min(dist_right);
            bottom_left_A = poly24.bnodes(min_left, :);
            bottom_right_A = poly24.bnodes(min_right, :);
            x = [bottom_left_A(1); bottom_right_A(1)];
            z = [bottom_left_A(2); bottom_right_A(2)];
            Lx = x(2) - x(1);
            fac = Lx/Lx_glob;
            Nx = ceil(polyA.fac_scale*fac*Nx_glob);
            dx = Lx/Nx;
            x_new = x(1) + cumsum(repmat(dx, Nx+1, 1)) - dx;

            z_new = interp1(x, z, x_new, 'linear'); 

            partA_top = findOverlappingNodes(polyA, poly.p6, 'top');
            partA_bottom = [x_new, z_new];
            p27_east = poly.p27.G.nodes.coords(poly.p27.east_mask, :);
            partA_west = p27_east(p27_east(:,2) >= bottom_left_A(2)-eps, :);
            p17_west = poly.p17.G.nodes.coords(poly.p17.west_mask, :);
            partA_east = p17_west(p17_west(:,2) >= bottom_right_A(2)-eps, :);

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.p24fA = polyA;

            % - part B -
            polyB = Faults(all_polys, 24, fac_scale, false);
            polyB.p_idx = 'p24fB';

            bottom_right_B = poly24.internal_top;
            % Define remaining parts of bottom side
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
            partB_west = poly24.external_west(poly24.external_west(:,2) >= left_pt(2)-eps & ... % bottom left of part B
                                                poly24.external_west(:,2) <= partA_bottom(1,2)+eps, :); % top left of part B
            partB_west = flip(partB_west); % to get in same order as coordinate masks
            partB_east = poly24.external_east(poly24.external_east(:,2) >= min(poly24.external_east(:,2))-eps & ...
                                                poly24.external_east(:,2) <= partA_bottom(end,2)+eps, :);
            
            polyB.internal_top = partB_top;
            polyB.internal_bottom = partB_bottom;
            polyB.internal_west = partB_west;
            polyB.internal_east = partB_east;

            polyB = Faults.distributeAndCollapseNodes(polyB, partB_west, partB_east, ...
                                                            partB_top, partB_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.p24fB = polyB; 

            % - part C -
            polyC = Faults(all_polys, 24, fac_scale, false);
            polyC.p_idx = 'p24fC';
            
            partC_top = partB_bottom(partB_bottom(:,1) <= poly24.internal_top(1,1)+eps, :);
          
            partC_bottom = poly24.internal_bottom;

            partC_west = poly24.external_west(poly24.external_west(:,2) >= partC_bottom(1,2)-eps & ...
                                              poly24.external_west(:,2) <= partC_top(1,2)+eps, :);
            partC_west = flip(partC_west);
            partC_east = poly24.internal_east;

            polyC.internal_top = partC_top;
            polyC.internal_bottom = partC_bottom;
            polyC.internal_west = partC_west;
            polyC.internal_east = partC_east;

            polyC = Faults.distributeAndCollapseNodes(polyC, partC_west, partC_east, ...
                                                            partC_top, partC_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.p24fC = polyC;

            % --- Small piece: poly9 ---
            poly9 = poly.p9f;
            
            % - part A -
            polyA = Faults(all_polys, 9, fac_scale, false);
            polyA.p_idx = 'p9fA';

            partA_top = poly9.internal_top;
            partA_bottom = poly9.internal_bottom;
            partA_west = poly24.internal_east(poly24.internal_east(:,2) >= partA_bottom(1,2)-eps, :);
            partA_east = poly9.external_east;

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

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
            partA_west = poly21.internal_west(poly21.internal_west(:,2) >= bottom_left_A(1,2)-eps, :);
            partA_east = poly21.external_east(poly21.external_east(:,2) >= bottom_right_A(end,2)-eps, :);

            polyA.internal_top = partA_top;
            polyA.internal_bottom = partA_bottom;
            polyA.internal_west = partA_west;
            polyA.internal_east = partA_east;

            polyA = Faults.distributeAndCollapseNodes(polyA, partA_west, partA_east, ...
                                                            partA_top, partA_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

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
            partB_top = partB_top(partB_top(:,1) <= bottom_right_A(1,1)+eps, :);

            partB_west = poly24.internal_east(poly24.internal_east(:,2) <= bottom_left_A(1,2)+eps, :);
            partB_east = p7small_west(p7small_west(:,2) >= bottom_right_B(2)-eps, :);

            polyB.internal_top = partB_top;
            polyB.internal_bottom = partB_bottom;
            polyB.internal_west = partB_west;
            polyB.internal_east = partB_east;

            polyB = Faults.distributeAndCollapseNodes(polyB, partB_west, partB_east, ...
                                                            partB_top, partB_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.p21fB = polyB;

            % - part C -
            polyC = Faults(all_polys, 21, fac_scale, false);
            polyC.p_idx = 'p21fC';

            partC_top = [poly24.internal_bottom; partB_bottom(2:end, :)]; % 2:end to avoid duplicating node at overlap
            %partC_west = findOverlappingNodes(polyC, poly.p8, 'west');
            partC_west = flip(poly21.external_west);

            p7sG_bottom = poly.p7small.G.nodes.coords(poly.p7small.bottom_mask, :);
            bottom_right_C = p7sG_bottom(1,:);
            dist = Faults.euclideanDist(bottom_right_C, partC_west);
            [~, min_idx] = min(dist);
            bottom_left_C = partC_west(min_idx, :);
            Lx_bottom_C = bottom_right_C(1) - bottom_left_C(1);
            Nx_bottom_C = ceil(Lx_bottom_C/Lx_glob * Nx_glob);
            bottom_C = [bottom_left_C; bottom_right_C];
            [x_interp, z_interp] = Faults.interpFault(bottom_C, Nx_bottom_C);           
            partC_bottom = [x_interp, z_interp];
            partC_west = partC_west(partC_west(:,2) >= partC_bottom(1,2)-eps, :);

            partC_east = p7small_west(p7small_west(:,2) <= partC_top(end,2)+eps, :);

            polyC.internal_top = partC_top;
            polyC.internal_bottom = partC_bottom;
            polyC.internal_west = partC_west;
            polyC.internal_east = partC_east;

            polyC = Faults.distributeAndCollapseNodes(polyC, partC_west, partC_east, ...
                                                            partC_top, partC_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.p21fC = polyC;

            % - part D -
            polyD = Faults(all_polys, 21, fac_scale, false);
            polyD.p_idx = 'p21fD';

            partD_top = [partC_bottom; p7sG_bottom(2:end,:)]; % 2:end to avoid duplicate nodes
            partD_bottom = poly21.internal_bottom;
            partD_west = poly21.external_west(poly21.external_west(:,2) <= partC_bottom(1,2)+eps, :);
            partD_west = flip(partD_west);
            partD_east = poly21.external_east(poly21.external_east(:,2) <= p7sG_bottom(end,2)+eps, :);

            polyD.internal_top = partD_top;
            polyD.internal_bottom = partD_bottom;
            polyD.internal_west = partD_west;
            polyD.internal_east = partD_east;

            polyD = Faults.distributeAndCollapseNodes(polyD, partD_west, partD_east, ...
                                                            partD_top, partD_bottom, ...
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);
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
            dist = Faults.euclideanDist(bottom_right_A, p3_east);
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
                                                            Lx_glob, Lz_glob, Nx_glob, Nz_glob);

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
            [~, unique_idx] = uniquetol(polyB.bnodes, 'ByRows',true);
            polyB.bnodes = polyB.bnodes(sort(unique_idx), :);

            poly.p23fB = polyB;
        end

        function poly = distributeAndCollapseNodes(poly, part_west, part_east, part_top, part_bottom, ...
                                                    Lx_glob, Lz_glob, Nx_glob, Nz_glob, varargin)
            [~, unique_idx] = uniquetol(part_west, 'ByRows', true); % remove duplicating overlapping nodes
            part_west = part_west(sort(unique_idx), :);
            [~, unique_idx] = uniquetol(part_east, 'ByRows', true);
            part_east = part_east(sort(unique_idx), :);          
            
            poly.top_side = part_top;
            poly.bottom_side = part_bottom;

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
                                         Lx_glob, Lz_glob, Nx_glob, Nz_glob);

            poly.G.i = nan(poly.G.cells.num, 1);
            poly.G.j = nan(poly.G.cells.num, 1);
        
            if Nx_sub > 2 % collapse internal stacks according to closest side
                Nz_sub = poly.G.cartDims(2)+1; % vertical dim may have increased if collapsed nodes added
                Nx_sub = poly.G.cartDims(1)+1; % horizontal dim may also have changed
                Nx_int = Nx_sub - 2;
                
                G_dum = poly.G;
                west_side = poly.G.nodes.coords(poly.west_mask, :);
                pt = unique(west_side, 'rows', 'stable');
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
                    poly.G.nodes.coords(vert_layer(2:end-1), 2) = z_interp(2:end-1); % ONLY SET Z-COORD
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

        function dist = euclideanDist(pt, line)
            dist = sqrt((pt(1) - line(:,1)).^2 + ...
                        (pt(2) - line(:,2)).^2);
        end

   end
end