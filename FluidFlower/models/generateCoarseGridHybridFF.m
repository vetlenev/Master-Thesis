function CG = generateCoarseGridHybridFF(G, VEGroups, varargin)
%Create a coarse hybrid VE grid from specified VE regions
%
% SYNOPSIS:
%   CG = generateCoarseGridMultiVE(G, VEGroups
%
% REQUIRED PARAMETERS:
%   G        - Aggregated grid from glue2DGrid
%
%   VEGroups - Indicator array with one entry per fine cell, containing the
%              index of the VE region where that cell belongs.
%
% RETURNS:
%   CG       - A modified MRST coarse grid with additional fields required
%              by the hybrid VE solvers.
%
% SEE ALSO:
%   convertToMultiVEModel

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('useDepth', true);
    opt = merge_options(opt, varargin{:});
    [G.cells.topDepth, G.cells.bottomDepth, G.cells.height, G.cells.topCells] = getCellHeights(G, opt);
    
    require coarsegrid
    r = zeros(G.cells.num, 1);
    p = zeros(G.cells.num, 1);
    
    isFine = VEGroups == 0;
    r(isFine) = 1;
        
    ii = G.i; % assuming G is global glued grid with aggregated i-indices from subgrids
    % Logical indexing follows cell ordering given by G.cells.indexMap
    for i = 1:max(VEGroups)
        g = VEGroups == i;
        r(g) = i + 1;
        p(g) = max(p) + ii(g) + 1;
    end
    p(isFine) = max(p) + (1:nnz(isFine))';
    p = processPartition(G, p);
    p = compressPartition(p);
    
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
    
    coarse = accumarray(p, r, [], @max);
    
    % Coarse cells, discretization type
    CG.cells.discretization = coarse;
    % Fine cells, discretization type
    CG.parent.cells.discretization = r;
    
    % Lateral connections between VE types should be set to zero    
    % Vertical connections between types should be flagged
    
    [CG.cells.topDepth, CG.cells.bottomDepth, CG.cells.height] = deal(zeros(CG.cells.num, 1));
    if opt.useDepth
        CG.cells.height = accumarray(CG.partition, CG.parent.cells.height, [CG.cells.num, 1], @sum);
        CG.cells.topDepth = accumarray(CG.partition, CG.parent.cells.topDepth, [CG.cells.num, 1], @max);
        CG.cells.bottomDepth = accumarray(CG.partition, CG.parent.cells.bottomDepth, [CG.cells.num, 1], @min);
    else
        for i = 1:CG.cells.num
            cells = CG.partition == i;
            CG.cells.height(i) = sum(CG.parent.cells.height(cells));
            % Alternate mode: We define depths from the top as metric
            % distance along the piecewise linear curve through the faces
            cells = find(cells);
            % Top and bottom is parametric distance from top
            td = max(CG.parent.cells.topDepth(cells));
            CG.cells.topDepth(i) = td;
            CG.cells.bottomDepth(i) = td - CG.cells.height(i);
            % Then fix top and bottom of cells as the cumulative sum along
            % the parametrization of the column
            d = CG.parent.cells.topDepth(cells);
            [~, sortIx] = sort(d);
            cells = cells(sortIx);
            t = cumsum(CG.parent.cells.height(cells));
            b = [0; t(1:end-1)];

            CG.parent.cells.topDepth(cells) = t + td;
            CG.parent.cells.bottomDepth(cells) = b + td;
            % Sort by depth
            % Assign new top and bottoms from this...
        end
    end
end

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