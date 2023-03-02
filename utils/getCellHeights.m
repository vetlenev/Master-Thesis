function [topz, bottomz, height, topcells] = getCellHeights(G, useDepth)
    faceNo = rldecode(1 : G.cells.num, ...
                    diff(G.cells.facePos), 2) .';
    
    hf_pt = G.faces.centroids(G.cells.faces(:, 1), :);
    
    hf_z = hf_pt(:, 3);
    hf_top = hf_z;
    hf_bottom = hf_z;
    if size(G.cells.faces, 2) > 1 && (any(strcmpi(G.type, 'cartGrid')) || any(strcmpi(G.type, 'processGRDECL')))
        ftype = G.cells.faces(:, 2);
    else
        fno = G.cells.faces(:, 1);
        ftype = ones(size(fno));
        cno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
        normals = abs(G.faces.normals(fno, :));
        nz = normals(:, 3);
        nxy = max(normals(:, 1), normals(:, 2));
        zc = G.cells.centroids(cno, 3);
        z = G.faces.centroids(fno, 3);
        
        ftype(nz > nxy & z < zc) = 5;
        ftype(nz > nxy & z > zc) = 6;
        warning('Grid is not corner point or cartGrid derived. Unable to guess column structure. Results may be inaccurate!');
    end
    % Avoid picking "wrong" type of face for columns
    maskTop    = ftype == 5;
    maskBottom = ftype == 6;

    hf_top(~maskTop) = inf;
    hf_bottom(~maskBottom) = -inf;

    topIx = accumarray(faceNo, hf_top, [G.cells.num, 1],  @minIndex);
    bottomIx = accumarray(faceNo, hf_bottom, [G.cells.num, 1],  @maxIndex);
    
    % We have local indexes, increment with the number of faces
    offsets = [0; cumsum(diff(G.cells.facePos(1:end-1)))];
    top = hf_pt(topIx + offsets, :);
    topSubs = topIx + offsets;
    topcells = G.faces.neighbors(G.cells.faces(topSubs, 1), :);
    cells = (1:G.cells.num)';
    for i = 1:2
        topcells(topcells(:, i) == cells, i) = 0;
    end
    topcells = sum(topcells, 2);
    topcells(topcells == 0) = cells(topcells == 0);
    
    bottom = hf_pt(bottomIx + offsets, :);
    topz = top(:, 3);
    bottomz = bottom(:, 3);
    if useDepth
        height = bottomz - topz;
    else
        height = sqrt(sum((top - bottom).^2, 2));
    end
    assert(all(height >= 0))
end

function ix = minIndex(x)
    [~, ix] = min(x);
end

function ix = maxIndex(x)
    [~, ix] = max(x);
end
