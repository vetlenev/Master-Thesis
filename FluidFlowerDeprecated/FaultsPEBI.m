classdef Faults < PolygonGrid

    properties
        fac_scale
        bfacies
        multiple

        internal_top
        internal_bottom
        internal_west
        internal_east

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
            if multiple_polys
                obj.p = [];
                obj.p_orig = [];
                obj.p_idx = 0;
                obj.facies = 0;                
            end
        end      
       
    end


   methods (Static)
       
       function obj = globalSubgrid(obj, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_nodes_overlap, varargin)
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
            Lx_new = Nx * dx_glob;                    
                    
            Lz_min = min(poly(:,2)); % 2D -> second element is depth-coord
            Lz_max = max(poly(:,2));
            Lz = Lz_max - Lz_min;
            
            fac = Lz/Lz_glob;
            Nz = ceil(fac * Nz_glob);
                          
            loop = [obj.bnodes; obj.bnodes(1,:)]; % add first node to create loop                        

            loop(:,1) = loop(:,1) - Lx_min;
            loop(:,2) = loop(:,2) - Lz_min;

            if strcmp(obj.p_idx, 'p25PEBI') % Remove faces at tip of fault -> too small to give PEBI cells
                loop = {loop(1:20,:), loop([20,24],:), loop(24:65,:), loop([65,69],:), loop(69:end-1,:), loop(end-1:end,:)};                
                % Update bnodes
                obj.bnodes = vertcat(loop{:});
                obj.bnodes = obj.bnodes(1:end-1,:); % last element same as first -> remove this
                obj.bnodes(:,1) = obj.bnodes(:,1) + Lx_min;
                obj.bnodes(:,2) = obj.bnodes(:,2) + Lz_min;
                res_size = Lz/80;
            elseif strcmp(obj.p_idx, 'pBFPEBI')                
                loop = {loop(1:62,:), loop(65:end,:)};
                % Update bnodes
                obj.bnodes = vertcat(loop{:});
                obj.bnodes = obj.bnodes(1:end-1,:); % last element same as first -> remove this
                obj.bnodes(:,1) = obj.bnodes(:,1) + Lx_min;
                obj.bnodes(:,2) = obj.bnodes(:,2) + Lz_min;
                res_size = Lz/90;
            else
                loop = {loop(1:end-1,:), loop(end-1:end,:)};
                res_size = Lz/80;
            end

            %obj.G = compositePebiGrid2D([Lx/20, Lz/20], [Lx, Lz], 'faceConstraints', loop);
            obj.G = pebiGrid2D(res_size, [1.2*Lx, 1.2*Lz], 'faceConstraints', loop, ...
                                                'FCFactor', 1, ...
                                                'FCRefinement', true, ...
                                                'circleFactor', 0.55);%, 'interpolateFC', true);
            
            % Shift nodes to lateral location of facies
            obj.G.nodes.coords(:,1) = obj.G.nodes.coords(:,1) + Lx_min;
            % Shift nodes to vertical location of facies
            obj.G.nodes.coords(:,2) = obj.G.nodes.coords(:,2) + Lz_min;

            obj.bottom_mask = obj.G.nodes.coords(:,2) == min(obj.G.nodes.coords(:,2));
            obj.top_mask = obj.G.nodes.coords(:,2) == max(obj.G.nodes.coords(:,2));
            obj.west_mask = obj.G.nodes.coords(:,1) == min(obj.G.nodes.coords(:,1));
            obj.east_mask = obj.G.nodes.coords(:,1) == max(obj.G.nodes.coords(:,1));
        end
       
     
        function obj = cellsInsideFault(obj)            
            obj.G = computeGeometry(obj.G);
            %boundary_coords = G.faces.centroids(G.faces.tag, :);            
            bnodes = obj.bnodes;

            xv = bnodes(:,1);
            zv = bnodes(:,2);
            xq = obj.G.cells.centroids(:,1);
            zq = obj.G.cells.centroids(:,2);

            inside_fault = inpolygon(xq,zq,xv,zv);
            
            [obj.G_fault, obj.gc, obj.gf, obj.gn] = extractSubgrid(obj.G, inside_fault);
        end

        function [poly, bneighbors] = assignFacies(poly, poly_neighbors)
            pG = poly.G;
            %pGF = poly.p31PEBI.G_fault;
            
            bneighbors = pG.faces.neighbors(pG.faces.tag, :);
            bfacies = zeros(size(bneighbors,1), 1);

            for i=1:size(bneighbors,1)
                d = 1;
                while d < 3 % go through first and second index
                    bn = bneighbors(i,d);
                    bnc = pG.cells.centroids(bn, :);
                    for j=1:size(poly_neighbors,1)
                        pn = poly_neighbors{j,1};
                        bnodes = [pn.G.nodes.coords(pn.top_mask,:); ...
                                     pn.G.nodes.coords(pn.west_mask,:); ...
                                     pn.G.nodes.coords(pn.bottom_mask,:); ...                         
                                     pn.G.nodes.coords(pn.east_mask,:)];
                        in_poly = inpolygon(bnc(:,1), bnc(:,2), bnodes(:,1), bnodes(:,2));
                        if in_poly
                            bfacies(i) = pn.facies;
                            d = 3;
                            break;
                        end
                    end
                    d = d + 1;
                end
            end

            poly.bfacies = bfacies;
        end

        function [poly, bneighbors] = assignFacies25(poly, poly_neighbors)
            pG = poly.G;
            %pGF = poly.p31PEBI.G_fault;
            
            bneighbors = pG.faces.neighbors(pG.faces.tag, :);
            bfacies = zeros(size(bneighbors,1), 1);

            pinch_map = struct('A', 0.1, 'B', 0.2, 'C', 0.3);

            for i=1:size(bneighbors,1)
                d = 1;
                while d < 3 % go through first and second index
                    bn = bneighbors(i,d);
                    bnc = pG.cells.centroids(bn, :);
                    for j=1:size(poly_neighbors,1)
                        pn = poly_neighbors{j,1};
                        bnodes = [pn.G.nodes.coords(pn.top_mask,:); ...
                                     pn.G.nodes.coords(pn.west_mask,:); ...
                                     pn.G.nodes.coords(pn.bottom_mask,:); ...                         
                                     pn.G.nodes.coords(pn.east_mask,:)];
                        in_poly = inpolygon(bnc(:,1), bnc(:,2), bnodes(:,1), bnodes(:,2));
                        if in_poly       
                            pinch_type = regexp(pn.p_idx, 'p\d*([ABC])', 'tokens');
                            if ~isempty(pinch_type)
                                sub_facie = pinch_map.(pinch_type{1});
                            else
                                sub_facie = 0;
                            end
                            bfacies(i) = pn.facies + sub_facie;
                            d = 3;
                            break;
                        end
                    end
                    d = d + 1;
                end
            end

            poly.bfacies = bfacies;
        end

        function poly = removeRedundantFaces(poly, bneighbors, tip_faces, varargin) 
            pG = poly.G;
            pGF = poly.G_fault;
            faces2keep = [];
            if nargin > 3
                faces2keep = varargin{1};
            end

            bcells_idx = ismember(poly.gc, bneighbors);
            bcells = find(bcells_idx); % local cells (for G_fault)
            bcells_glob = poly.gc(bcells_idx); % global cells (for G)
            bfaces_glob = cell(numel(bcells), 1);
            for i=1:numel(bcells)
                c = bcells_glob(i);
                bfaces_glob{i} = pG.cells.faces(pG.cells.facePos(c):pG.cells.facePos(c+1)-1, :);
            end
            bfaces_glob = vertcat(bfaces_glob{:});
            ubfaces = unique(bfaces_glob);
            edge_faces_glob = ubfaces(histc(bfaces_glob, ubfaces) > 1);
            edge_faces = find(ismember(poly.gf, edge_faces_glob)); % global -> local mapping
            
            % If face inside fault does NOT share node with neighboring polygon, 
            % AND neighboring cells are from different facies,
            % add this face to list of faces that should be removed.
            remove_faces = [];
            bfacies = poly.bfacies;
            for i=1:numel(edge_faces)                
                f = edge_faces(i);
                f_glob = edge_faces_glob(i);
            
                bnodes = pGF.faces.nodes(pGF.faces.nodePos(f):pGF.faces.nodePos(f+1)-1, :);
                bnodes = pGF.nodes.coords(bnodes, :);
                shared_node = ismembertol(bnodes, poly.bnodes, 1e-5, 'ByRows', true);
            
                bneigh = pG.faces.neighbors(f_glob, :);
                bn_idx1 = any(ismember(bneighbors, bneigh(1)), 2);
                bn_idx2 = any(ismember(bneighbors, bneigh(2)), 2);
            
                % Tip cell connected to multiple facies - must be removed manually
                if nnz(bn_idx1) > 1        
                    bn_idx1 = bn_idx2;
                elseif nnz(bn_idx2) > 1
                    bn_idx2 = bn_idx1;
                end
                
                if ~any(shared_node) && bfacies(bn_idx1) == bfacies(bn_idx2) ...
                        && ~ismember(f, faces2keep)
                    remove_faces = [remove_faces; f];
                end
            end
            
            % Remove provided faces at tip of fault
            remove_faces = unique([remove_faces; tip_faces]);
            
            % Remove discontinuous face transitions
            for i=1:numel(remove_faces)
                [pGF, remove_faces] = mergeCellsPEBI(pGF, remove_faces);
            end
            poly.G_fault = pGF;
        end

        function poly = fixUnmatchedFaces(poly)
            pGF = poly.G_fault;

            top_pts = poly.bnodes(1:2,:);
            pt = pGF.faces.centroids;
            d = Faults.dist_to_line(pt, top_pts(1,:), top_pts(2,:));
            faces_fix = find(d < 1e-3);
            
            nodes_fix = zeros(numel(faces_fix)*2, 2); % *2 since each of the faces have two nodes
            nodes_idx = zeros(numel(faces_fix)*2, 1);
            for i=1:numel(faces_fix)
                f = faces_fix(i);
                node_idx = pGF.faces.nodes(pGF.faces.nodePos(f):pGF.faces.nodePos(f+1)-1, :);
                nodes_fix(i*2-1:i*2, :) = pGF.nodes.coords(node_idx, :);
                nodes_idx(i*2-1:i*2) = node_idx;
            end
            
            tol = 1e-4; % tolerance needed in case of very small faces
            nodes_q = nodes_fix(nodes_fix(:,1) > min(nodes_fix(:,1)) & ...
                                nodes_fix(:,1) < max(nodes_fix(:,1)), :);
            nodes_qidx = nodes_idx(nodes_fix(:,1) > min(nodes_fix(:,1)) & ...
                                nodes_fix(:,1) < max(nodes_fix(:,1)));
            
            [nodes_q, qidx] = unique(nodes_q, 'rows');
            nodes_qidx = nodes_qidx(qidx);
            
            z_new = interp1(top_pts(:,1), top_pts(:,2), nodes_q(:,1));
            nodes_q = [nodes_q(:,1), z_new];
            
            pGF.nodes.coords(nodes_qidx, :) = nodes_q;

            pGF = computeGeometry(pGF);
            poly.G_fault = pGF;
        end

        function poly = fixUnmatchedFaces25(poly, point_set)
            pGF = poly.G_fault;
            pt = pGF.faces.centroids;

            for p=1:numel(point_set)
                pts = point_set{p};
                line_pts = poly.bnodes(pts,:);                
                
                if p == 1
                    d = Faults.dist_to_line(pt, line_pts(1,:), line_pts(2,:));
                    faces_fix = find(d < 1e-3);
                elseif p==2
                    faces_fix = find(pGF.faces.centroids(:,1) > 1.3694 & ...
                                      pGF.faces.centroids(:,1) < 1.3696 & ...
                                      pGF.faces.centroids(:,2) > 0.675 & ...
                                      pGF.faces.centroids(:,2) < 0.677);
                elseif p == 3
                    faces_fix = find(pGF.faces.centroids(:,1) > 1.3705 & ...
                                      pGF.faces.centroids(:,1) < 1.371 & ...
                                      pGF.faces.centroids(:,2) > 0.671 & ...
                                      pGF.faces.centroids(:,2) < 0.673);
                end
                
                nodes_fix = zeros(numel(faces_fix)*2, 2); % *2 since each of the faces have two nodes
                nodes_idx = zeros(numel(faces_fix)*2, 1);
                for i=1:numel(faces_fix)
                    f = faces_fix(i);
                    node_idx = pGF.faces.nodes(pGF.faces.nodePos(f):pGF.faces.nodePos(f+1)-1, :);
                    nodes_fix(i*2-1:i*2, :) = pGF.nodes.coords(node_idx, :);
                    nodes_idx(i*2-1:i*2) = node_idx;
                end
                
                tol = 1e-4; % tolerance needed in case of very small faces
                % CHANGE THIS (nodes_q and nodes_qidx) !! Can't use
                % min(nodes_fix(:,1)) -> only applicable for top side
                nodes_q = nodes_fix(nodes_fix(:,1) > min(line_pts(:,1))+tol & ...
                                    nodes_fix(:,1) < max(line_pts(:,1))-tol, :);
                nodes_qidx = nodes_idx(nodes_fix(:,1) > min(line_pts(:,1))+tol & ...
                                    nodes_fix(:,1) < max(line_pts(:,1))-tol);
                
                [nodes_q, qidx] = unique(nodes_q, 'rows');
                nodes_qidx = nodes_qidx(qidx);
                
                z_new = interp1(line_pts(:,1), line_pts(:,2), nodes_q(:,1));
                nodes_q = [nodes_q(:,1), z_new];
                
                pGF.nodes.coords(nodes_qidx, :) = nodes_q;
            end

            pGF = computeGeometry(pGF);
            poly.G_fault = pGF;
        end

        function poly = changeNodesForFaces(poly, face, true_coord_left, true_coord_right)
            pGF = poly.G_fault;
            node_idx = pGF.faces.nodes(pGF.faces.nodePos(face):pGF.faces.nodePos(face+1)-1, :);
            old_coords = pGF.nodes.coords(node_idx,:);
            [~, row_idx] = sortrows(old_coords);
            
            new_node_left = find(sqrt(sum((pGF.nodes.coords - true_coord_left).^2, 2)) < 1e-4);
            new_node_right = find(sqrt(sum((pGF.nodes.coords - true_coord_right).^2, 2)) < 1e-4);
            
            new_nodes = [new_node_left; new_node_right];
            %pGF.faces.nodes(pGF.faces.nodePos(old_face):pGF.faces.nodePos(old_face+1)-1, :) = new_nodes(row_idx);
            pGF.nodes.coords(node_idx(row_idx),:) = pGF.nodes.coords(new_nodes,:);
            poly.G_fault = pGF;
        end
        
        function d = dist_to_line(pt, p1, p2)   
            % Add dummy y-coord
            pt = [pt(:,1), ones(size(pt,1),1), pt(:,2)];
            p1 = [p1(1) 1 p1(2)];
            p2 = [p2(1) 1 p2(2)]; 
        
            p1 = repmat(p1, size(pt,1), 1);
            p2 = repmat(p2, size(pt,1), 1);
        
            line = p2 - p1;
            pt2 = pt - p1;
            d = sqrt(sum(cross(line,pt2,2).^2,2)) ./ sqrt(sum(line.^2,2));
        end

        function pts_inside = uniformPointDistribution(poly, num_pts)
            % Uniform distribution of points inside polygon bounded by
            % bnodes of grid G.            
            G = poly.G;
            p = poly.bnodes;

            G = computeGeometry(G);
            x_min = min(G.cells.centroids(:,1));
            x_max = max(G.cells.centroids(:,1));
            z_min = min(G.cells.centroids(:,2));
            z_max = max(G.cells.centroids(:,2));
            xx = linspace(x_min, x_max, num_pts);
            zz = linspace(z_min, z_max, num_pts);
            
            [X, Z] = meshgrid(xx,zz);
            points = [X(:), Z(:)];

            inpoly = inpolygon(points(:,1), points(:,2), p(:,1), p(:,2));
            pts_inside = [points(inpoly,1), points(inpoly,2)];
        end
    
        
   
   end
end