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
