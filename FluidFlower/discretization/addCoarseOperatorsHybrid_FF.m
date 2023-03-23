function model = addCoarseOperatorsHybrid_FF(model)
% Add various extra coarse-scale information for a multi-VE model
%
% SYNOPSIS:
%   model = addCoarseOperatorsMultiVE(model)
%
% REQUIRED PARAMETERS:
%   model - Output from generateCoarseGridMultiVE
%
% OPTIONAL PARAMETERS:
%   None.
%
% RETURNS:
%   model - A modified model.
%
% SEE ALSO:
%   addCoarseOperatorsMultiVE

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

    G = model.G;
    if ~isfield(G, 'partition') || ~isfield(G.cells, 'discretization')
        error('Function only applicable to coarse hybrid VE grids from generateCoarseGridMultiVE');
    end
    N = model.operators.N;
    n1 = N(:, 1);
    n2 = N(:, 2);

    d = G.cells.discretization;
    % Change in discretization type
    transition = d(n1) ~= d(n2);

    % Internal connections inside same VE zone
    model.operators.connections.veInternalConn = ...
        d(n1) > 1 & ~transition;
    % Between two fine cells
    model.operators.connections.fineInternalConn =...
        d(n1) == 1 & ~transition;
    % Between two VE zones horizontally
    ve_trans = transition & d(n1) > 1 & d(n2) > 1; % > 0 & > 0
    diff_column = G.cells.columns(n1) ~= G.cells.columns(n2);
    model.operators.connections.veTransitionHorizontalConn = ...
        ve_trans & diff_column;
    % Between two VE zones vertically
    model.operators.connections.veTransitionVerticalConn = ...
        ve_trans & ~diff_column;
    % From one VE zone to fine
    model.operators.connections.veToFineConn = ...
        transition & (d(n1) == 1 | d(n2) == 1);
    % From VE zone to fine horizontally
    model.operators.connections.veToFineHorizontal = ...
        transition & (d(n1) == 1 | d(n2) == 1) & diff_column;
    % From VE zone to fine vertically
    model.operators.connections.veToFineVertical = ...
        transition & (d(n1) == 1 | d(n2) == 1) & ~diff_column;

    globalCoarseFace = find(all(G.faces.neighbors > 0, 2));

    nt = numel(globalCoarseFace);

    transFaceTop = zeros(nt, 2);
    transFaceBottom = zeros(nt, 2);
    faceHeights = zeros(nt, 2);
    for i = 1:nt
        f = globalCoarseFace(i); % connection number i
        fine = G.faces.fconn(G.faces.connPos(f):G.faces.connPos(f+1)-1); % all fine faces associated with coarse connection f
        C = G.faces.neighbors(f, :); % blocks connected to connection f

        c0 = G.parent.faces.neighbors(fine, :); % fine cells connected to associated fine faces
        blockNo = G.partition(c0); % block indices of given fine cells
        if numel(fine) == 1
            blockNo = blockNo';
        end
        c = c0;

        sorted = blockNo(:, 1) == C(1); % C(1) -> only one coarse face considered per loop iter?
        c(:, 1) = sorted.*c0(:, 1) + ~sorted.*c0(:, 2);
        c(:, 2) = sorted.*c0(:, 2) + ~sorted.*c0(:, 1);
        if numel(fine) == 1
            t = G.parent.cells.topDepth(c)';
            b = G.parent.cells.bottomDepth(c)';            
        else
            t = max(G.parent.cells.topDepth(c)); % NB: changed min -> max when z is height, not depth
            b = min(G.parent.cells.bottomDepth(c)); % max -> min
        end
        H = G.cells.height(C)';
        B = G.cells.bottomDepth(C)';
        T = G.cells.topDepth(C)';
    
        % Think this is ok - scale the height by the fraction.
        % Sort of assumes that columns are vector-like with no
        % deviations.
        faceHeights(i, :) = H.*(t - b)./(T - B);

        transFaceTop(i, :) = t;
        transFaceBottom(i, :) = b;
    end
    % Top of the columns connected to the face
    model.operators.connections.faceTopDepth = transFaceTop;
    % Bottom of the columns connected to the face
    model.operators.connections.faceBottomDepth = transFaceBottom;
    % Height of the cell-column connected to the specific face
    model.operators.connections.faceHeight = faceHeights;
    
    % Partition connections for parent grid
end