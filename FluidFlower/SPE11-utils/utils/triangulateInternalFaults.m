function [poly_obj, nodes_overlap] = triangulateInternalFaults(poly_obj, poly_num, ...
                                                                nodes_overlap, pts_overlap, G_glob, varargin)  
   
% Generate nodes for internal polygon of bottom fault.
    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,1));
    N_glob = G_glob.cartDims;
    Nx_glob = N_glob(1); Nz_glob = N_glob(2);

    p_idx_pebi = strcat('p', string(poly_num), 'PEBI'); 
    poly = poly_obj.(p_idx_pebi);
    external_nodes = poly.bnodes;

    internal_pts = poly.p(~ismembertol(poly.p, poly_obj.pBFPEBI.p, 'ByRows', true), :);
    internal_nodes = [];
    % Interpolate nodes for internal sides
    if poly_num == 24
        % bottom side
        ext_swap = find(external_nodes(:,1) == min(external_nodes(:,1))); % swap from external to internal nodes oocurs at min x-coord
        p_bottom = [external_nodes(ext_swap,:); internal_pts(1,:)];
        [x,z] = interpolateInternalSide(poly, p_bottom, Lx_glob, Nx_glob, 1); % z-interpolation horizontally along bottom side        
        poly.internal_bottom = [x,z];
        internal_nodes = [internal_nodes; poly.internal_bottom];       
        % internal east side
        p_east = internal_pts;
        [x,z] = interpolateInternalSide(poly, p_east, Lz_glob, Nz_glob, 2); % x-interpolation vertically along east side
        poly.internal_east = [x,z];
        internal_nodes = [internal_nodes; poly.internal_east];
        % internal top side
        p_top = [internal_pts(end,:); external_nodes(ext_swap+1,:)];
        [x,z] = interpolateInternalSide(poly, p_top, Lx_glob, Nx_glob, 1); % z-interpolation horizontally
        poly.internal_top = [x,z];
        internal_nodes = [internal_nodes; poly.internal_top];

        % update boundary nodes for internal polygon to include internal points
        poly.bnodes = [external_nodes(1:ext_swap,:); internal_nodes; external_nodes(ext_swap+1:end,:)];

        nodes_overlap.(p_idx_pebi) = double.empty(0,2); % no overlap with other internal polygons, since nodes haven't been assigned to the others yet
    elseif poly_num == 9
        % top side -> copy nodes from p24
        poly.internal_top = findInternalOverlap(poly, poly_obj.p24PEBI, 'top');
        internal_nodes = [internal_nodes; flip(poly.internal_top)];
        % west side -> copy associated nodes from p24 (east side of p24)
        poly.internal_west = findInternalOverlap(poly, poly_obj.p24PEBI, 'east');
        internal_nodes = [internal_nodes; flip(poly.internal_west)]; 
        % bottom side
        p_bottom = [internal_pts(1,:); external_nodes(1,:)];
        [x,z] = interpolateInternalSide(poly, p_bottom, Lx_glob, Nx_glob, 1); % z-interpolation horizontally along bottom side        
        poly.internal_bottom = [x,z];
        internal_nodes = [internal_nodes; poly.internal_bottom];              

        poly.bnodes = [internal_nodes; external_nodes];
        nodes_overlap.(p_idx_pebi) = [poly.internal_top; poly.internal_west];
    elseif poly_num == 21
        % upper top side -> copy from bottom of p9
        internal_top_upper = findInternalOverlap(poly, poly_obj.p9PEBI, 'bottom');
        internal_nodes = [internal_nodes; flip(internal_top_upper)];
        % west side -> overlap with p24
        poly.internal_west = findInternalOverlap(poly, poly_obj.p24PEBI, 'east');
        internal_nodes = [internal_nodes; flip(poly.internal_west)];
        % lower top side -> copy from bottom of p24
        internal_top_lower = findInternalOverlap(poly, poly_obj.p24PEBI, 'bottom');
        poly.internal_top = [internal_top_upper; internal_top_lower];
        internal_nodes = [internal_nodes; flip(internal_top_lower)];
        % bottom side: no nodes -> interpolate
        p_added = poly_obj.p26PEBI.p(ismember(poly_obj.p26PEBI.p, poly_obj.p22PEBI.p, 'rows'), :);
        p_added = p_added(p_added(:,2) == max(p_added(:,2)), :); % this node is for some reason not included in the upper neighbor ... add manually
        ext_swap = find(external_nodes(:,1) == min(external_nodes(:,1))); % for this fault we know node with smallest x-coord gives correct index for external swap
        p_bottom = [external_nodes(ext_swap,:); p_added; external_nodes(ext_swap+1,:)];
        [x,z] = interpolateInternalSide(poly, p_bottom, Lx_glob, Nx_glob, 1); % z-interpolation horizontally along bottom side        
        poly.internal_bottom = [x,z];
        internal_nodes = [internal_nodes; poly.internal_bottom]; 
          
        int_swap = find(ismembertol(internal_nodes, poly.bnodes(1,:), 'ByRows',true));        
        poly.bnodes = [internal_nodes(1:int_swap,:); external_nodes(1:ext_swap,:); internal_nodes(int_swap+1:end,:); external_nodes(ext_swap+1:end,:)];
        nodes_overlap.(p_idx_pebi) = [poly.internal_top; poly.internal_west];
    elseif poly_num == 26
        % top side -> copy from associated part of 21
        poly.internal_top = findInternalOverlap(poly, poly_obj.p21PEBI, 'bottom');
        internal_nodes = [internal_nodes; flip(poly.internal_top)];
        % west side -> only external nodes
        % bottom side -> interpolate
        p_bottom = [external_nodes(end,:); internal_pts(end,:)];
        [x,z] = interpolateInternalSide(poly, p_bottom, Lx_glob, Nx_glob, 1); % z-interpolation horizontally along bottom side        
        poly.internal_bottom = [x,z];
        internal_nodes = [internal_nodes; poly.internal_bottom]; 
        % east side -> interpolate
        p_east = flip(internal_pts(1:end,:));
        [x,z] = interpolateInternalSide(poly, p_east, Lz_glob, Nz_glob, 2);
        poly.internal_east = [x,z];
        internal_nodes = [internal_nodes; poly.internal_east]; 

        int_swap = find(ismembertol(internal_nodes(:,1), max(external_nodes(:,1)), 'ByRows',true));
        poly.bnodes = [internal_nodes(1:int_swap,:); external_nodes; internal_nodes(int_swap+1:end,:)];
        nodes_overlap.(p_idx_pebi) = poly.internal_top;
    elseif poly_num == 22
        % top side -> copy from p21
        poly.internal_top = findInternalOverlap(poly, poly_obj.p21PEBI, 'bottom');
        internal_nodes = [internal_nodes; flip(poly.internal_top)];
        % west side -> copy from p26
        poly.internal_west = findInternalOverlap(poly, poly_obj.p26PEBI, 'east');
        internal_nodes = [internal_nodes; flip(poly.internal_west)];
        % bottom side -> interpolate
        p_bottom = [internal_pts(end,:); external_nodes(1,:)];
        [x,z] = interpolateInternalSide(poly, p_bottom, Lx_glob, Nx_glob, 1);
        poly.internal_bottom = [x,z];
        internal_nodes = [internal_nodes; poly.internal_bottom]; 
        % east side -> only externals

        poly.bnodes = [internal_nodes; external_nodes];
        nodes_overlap.(p_idx_pebi) = [poly.internal_top; poly.internal_west];
    elseif poly_num == 23
        % top side -> copy from p26 AND p22
        internal_top_left = findInternalOverlap(poly, poly_obj.p26PEBI, 'bottom');
        internal_top_right = findInternalOverlap(poly, poly_obj.p22PEBI, 'bottom');
        internal_top = [flip(internal_top_right); flip(internal_top_left)];
        [~, unique_idx] = uniquetol(internal_top, 'ByRows',true);
        poly.internal_top = internal_top(sort(unique_idx), :);
        internal_nodes = [internal_nodes; poly.internal_top];
        % west/bottom/east side -> only externals

        poly.bnodes = [internal_nodes; external_nodes];
        nodes_overlap.(p_idx_pebi) = poly.internal_top;
    end
 
    [~, unique_idx] = uniquetol(poly.bnodes, 'ByRows',true);
    poly.bnodes = poly.bnodes(sort(unique_idx), :);

    poly_obj.(p_idx_pebi) = poly;    
end

function [x_all,z_all] = interpolateInternalSide(poly, p_side_all, L_glob, N_glob, x_or_z)
    % Linearly interpolate interior points on internal side.
    % NB: Assumes p_side is order from lower to higher value (x or z)
    x_all = [];
    z_all = [];
    for i=1:size(p_side_all,1)-1 % loop through each line segment
        p_side = p_side_all(i:i+1,:);
        p_min = p_side(p_side(:,x_or_z) == min(p_side(:,x_or_z)), :);
        p_max = p_side(p_side(:,x_or_z) == max(p_side(:,x_or_z)), :);
        L = p_max(:,x_or_z) - p_min(:,x_or_z);
        fac = L/L_glob;
        N = ceil(poly.fac_scale*fac*N_glob);
        if all(N == 1) || all(all(p_min == p_max)) % two all's to check over rows and cols
            x = p_side(:,1);
            z = p_side(:,2);
        else        
            dX = L/N;
            xv = p_side(:,1);
            zv = p_side(:,2);
            if x_or_z == 1 % interpolate z-values horizontally
                x_interp = xv(1) + cumsum(repmat(dX, N+1, 1)) - dX;
                z_interp = interp1(xv, zv, x_interp);
            elseif x_or_z == 2 % interpolate x-values vertically
                z_interp = zv(1) + cumsum(repmat(dX, N+1, 1)) - dX;       
                x_interp = interp1(zv, xv, z_interp);
            else
                error('Only two dimensions allowed.')
            end
        
            x = [xv(1); x_interp(2:end-1); xv(end)];
            z = [zv(1); z_interp(2:end-1); zv(end)];  
        end
        x_all = [x_all; x];
        z_all = [z_all; z];       
    end

    [~, unique_idx] = uniquetol([x_all,z_all], 'ByRows',true);
    x_all = x_all(sort(unique_idx), :);
    z_all = z_all(sort(unique_idx), :);

end

function nodes = findInternalOverlap(poly, poly_other, side)
    % Find nodes from neighboring polygon at given side.
    % NB: Parameter 'side'  is the side from neighboring polygon.
    p = poly.p;

    if strcmp(side, 'east')
        nodes_other = poly_other.internal_east;
        x_or_z = 2;
    elseif strcmp(side, 'top')
        nodes_other = poly_other.internal_top;
        x_or_z = 1;
    elseif strcmp(side, 'west')
        nodes_other = poly_other.internal_west;
        x_or_z = 2;
    elseif strcmp(side, 'bottom')
        nodes_other = poly_other.internal_bottom;
        x_or_z = 1;
    end

    p_overlap = p(ismembertol(p, nodes_other, 'ByRows',true), :);
    min_p = p_overlap(p_overlap(:,x_or_z) == min(p_overlap(:,x_or_z)), :);
    max_p = p_overlap(p_overlap(:,x_or_z) == max(p_overlap(:,x_or_z)), :);
    nodes = nodes_other(nodes_other(:,x_or_z) >= min_p(:,x_or_z)-eps & ...
                        nodes_other(:,x_or_z) <= max_p(:,x_or_z)+eps, :);

end