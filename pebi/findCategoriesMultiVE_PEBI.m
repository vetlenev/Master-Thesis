function [categories, c_h, colNo] = findCategoriesMultiVE_PEBI(G, trans, ve_cols, fine_cols)
%Classify cells according to VE region
%
% SYNOPSIS:
%   [categories, c_h, colNo] = findCategoriesMultiVE(G, trans)
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   trans - Transmissibility. Only zero/non-zero information is used.
%
% RETURNS:
%   categories - Classification of each cell in the fine grid into VE
%                regions, where the same region indicator for a pair of
%                cells mean that they belong to the same horizontal flow
%                region where (potentially) VE can be valid
%
%   c_h        - Indicator of which column subset each entry belongs to. A
%                column may be split into multiple parts if the
%                transmissibility is zero, or geometry connectivity is
%                missing.
%
%   colNo      - Column number of each cell (uses IJK indexing, no account
%                for disconnected layers etc)
%

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
    
    mrstModule add matlab_bgl coarsegrid
    require matlab_bgl coarsegrid
    N = getNeighbourship(G);
    n1 = N(:, 1);
    n2 = N(:, 2);
    
    act = all(G.faces.neighbors > 0, 2); 
    %[ii, jj, kk] = gridLogicalIndices(G);
    
    % Make graph with horizontal connections removed
    T_h = trans(act);
    %T_h(ii(n1) ~= ii(n2)) = 0; % only for structured grids
    %T_h(jj(n1) ~= jj(n2)) = 0;
    xx = G.cells.centroids(:,1);
    %yy = G.cells.centroids(:,2); % uncomment for 3D grid
    dist_tol = @(x,y,tol) abs(x - y) < tol;
    horz_conn = ~dist_tol(xx(n1),xx(n2),1e-7);
    T_h(horz_conn) = 0; % NB: will set trans zero for pebi cells since these are not structurally aligned vertically
    %T_h(dist_tol(yy(n1),yy(n2),1e-7)) = 0; % uncomment for 3D grid
    %T_h(G.faces.tag > 0) = 0; % sealing faces
    
    % Compute graph partitioning. Gives us every column, divided up by
    % sealing transmissibilities or missing connections
    n = G.cells.num;
    A = sparse(n1, n2, T_h, n, n);
    c_h = components(A + A');
    
    % Generate tenative coarse grid
    CG = generateCoarseGrid(G, c_h);
    %colNo = ii + G.cartDims(1).*jj;
    ve_count = sum(cellfun(@numel, ve_cols));
    colNo = zeros(G.cells.num, 1); % zero for fine cells
    for i=1:numel(ve_cols)
        colNo(ve_cols{i}) = i;
    end        
    fine_cells = true(G.cells.num, 1);
    all_ve_cols = vertcat(ve_cols{:});
    fine_cells(all_ve_cols) = false;
    colNo(fine_cells) = mcolon(max(colNo)+1, max(colNo)+nnz(fine_cells));
    for i=1:size(fine_cols, 1)
        colNo(fine_cols(i,1)) = fine_cols(i,2);
    end
    
    Nc = getNeighbourship(CG);
    faceno = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2)';
    % Map into columns
    coarseColNo = zeros(CG.cells.num, 1);
    coarseColNo(CG.partition) = colNo;

    % For each coarse block, set connections to zero if zero on fine scale
    % or if that block is connected to multiple blocks in the same column
    %trans_c = trans_test(act);
    %trans_c(horz_conn) = 0;

    T_c_all = accumarray(faceno, trans(CG.faces.fconn));
    T_c_int = T_c_all(all(CG.faces.neighbors > 0, 2));  
    
    for i = 1:CG.cells.num
        act = Nc(:, 1) == i | Nc(:, 2) == i;

        nn = Nc(act, :);
        cn = coarseColNo(nn);
        if size(nn, 1) == 1
            cn = cn';
        end
        cn = sort(cn, 2);
        [C, JC, IC] = unique(cn, 'rows');
        
        if i == 173
            test = 0;
        end
        rem = find(accumarray(IC, 1) > 1);
        if any(rem)
            fa = find(act);
            for j = 1:numel(rem)
                T_c_int(fa(IC == rem(j))) = 0;
            end
        end
    end
    A = sparse(Nc(:, 1), Nc(:, 2), T_c_int, CG.cells.num, CG.cells.num);
    
    cC = components(A + A');
    categories = cC(CG.partition);
    categories = compressPartition(categories);
end