function [topz, bottomz, height, topcells] = getCellHeights(G, opt)
    cNo = rldecode(1 : G.cells.num, ...
                    diff(G.cells.facePos), 2) .';
    
    hf_pt = G.faces.centroids(G.cells.faces(:, 1), :);
    
    hf_z = hf_pt(:, 2);
    hf_top = hf_z;
    hf_bottom = hf_z;
    if size(G.cells.faces, 2) > 1 && (any(strcmpi(G.type, 'cartGrid')) || any(strcmpi(G.type, 'processGRDECL')))
        ftype = G.cells.faces(:, 2);
        % Manually estimate ftype for triangulated cells
        tri_faces = isnan(ftype);
        fno_tri = G.cells.faces(tri_faces, 1);       
        ftype_tri = ones(size(fno_tri));    
        cNo_tri = cNo(tri_faces);       
        normals = abs(G.faces.normals(fno_tri, :));
        nz = normals(:, 2);
        nx = normals(:, 1);
        zc = G.cells.centroids(cNo_tri, 2); % same centroid compared for all faces of a cell
        z = G.faces.centroids(fno_tri, 2);
        
        ftype_tri(nz > nx & z < zc) = 3;
        ftype_tri(nz > nx & z > zc) = 4;
        ftype(tri_faces) = ftype_tri;      
    else
        % Manually estimate ftype for all faces
        fno = G.cells.faces(:,1);
        ftype = ones(size(fno));
        normals = abs(G.faces.normals(fno, :));
        nz = normals(:, 2);
        nx = normals(:, 1);
        zc = G.cells.centroids(cNo, 2); % same centroid compared for all faces of a cell
        z = G.faces.centroids(fno, 2);
        
        ftype(nz > nx & z < zc) = 3;
        ftype(nz > nx & z > zc) = 4;        
        warning('Grid is not corner point or cartGrid derived. Unable to guess column structure. Results may be inaccurate!');
    end
    % Avoid picking "wrong" type of face for columns
    maskTop    = ftype == 4; % 3
    maskBottom = ftype == 3; % 4

    hf_top(~maskTop) = -inf;
    hf_bottom(~maskBottom) = inf; 

    topIx = accumarray(cNo, hf_top, [G.cells.num, 1],  @maxIndex); % accumulate over all faces belonging to given cell
    bottomIx = accumarray(cNo, hf_bottom, [G.cells.num, 1],  @minIndex); % NB: minIndex, not maxIndex, because z-coord is HEIGHT, not DEPTH
    
    % We have local indexes, increment with the number of faces
    offsets = [0; cumsum(diff(G.cells.facePos(1:end-1)))];
    top = hf_pt(topIx + offsets, :); % all face centroids for half-faces corresponding to top cells
    topSubs = topIx + offsets;
    topcells = G.faces.neighbors(G.cells.faces(topSubs, 1), :);
    cells = (1:G.cells.num)';
    for i = 1:2
        topcells(topcells(:, i) == cells, i) = 0;
    end
    topcells = sum(topcells, 2);
    topcells(topcells == 0) = cells(topcells == 0);
    
    bottom = hf_pt(bottomIx + offsets, :);
    topz = top(:, 2);
    bottomz = bottom(:, 2);

    if opt.useDepth
        %height = bottomz - topz;
        height = topz - bottomz;
    else
        height = sqrt(sum((top - bottom).^2, 2));
    end
    
    triangle_cells = diff(G.cells.facePos) == 3;
    imperm_cells = G.facies == 7;
    nan_cells = isnan(G.i);
    structured_flow = ~triangle_cells & ~imperm_cells & ~nan_cells; % all structured (cartesian) cells with nonzero fluid flow
    assert(all(height(structured_flow) >= 0))
    height(height < 0) = 0;
    %height = abs(height);
end

function ix = minIndex(x)
    [~, ix] = min(x);
end

function ix = maxIndex(x)
    [~, ix] = max(x);
end