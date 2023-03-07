function [poroG, permG, k] = poroPermToGrid(G, M, f, sand, clay)
%
%
%

% 1) Map diagonals to Grid indexing: Grid indexing starts at bottom left,
%    columns (x) move faster. Matrix starts counting at top left and rows
%    (y) move faster. This is for 2D grid. For 3D grid, the x-z plane
%    starts at top left, so we just need to transpose.

% figure(9); spy(sparse(M.vals))       % check that smears are correctly placed in M.
% 2D grid (xy) mapping
isSmear = reshape(transpose(flipud(M.vals)), G.cells.num, 1);
whichUnit = reshape(transpose(flipud(M.units)), G.cells.num, 1);

% 2) Assign permeability and porosity. Note that we receive FZ sand/clay poro
%    and permeability, the latter given perpendicular and parallel to the
%    to the bedding. This is then directly assumed to be the permeability
%    along/across the clay/sand smears (principal directions). This diagonal
%    permeability tensor needs to be rotated to the fault local axes,
%    since clay/sand smears are not parallel to the fault.

% Poro
poroG = ones(G.cells.num, 1)*sand.poro;
poroG(isSmear) = clay.poro;

% Perm
T = [cosd(f.alpha) sind(f.alpha); -sind(f.alpha) cosd(f.alpha)];
nk = numel(clay.perm(: ,1));
if nk > 1                       % we have different perms for each smear.
    assert(nk == sum(M.isclayIn));
    % sand
    sPermsGrid = mean(sand.perm, 1);
    S = diag([sPermsGrid(1), sPermsGrid(2)]);
    kmat = transformTensor(S, T, 'posDef');
    k.sand = [kmat(1,1), kmat(1,2), kmat(2,2)];
    permG = repmat(k.sand, [G.cells.num, 1]);
    
    % clay
    cPermsGrid = clay.perm;
    if isfield(M, 'idSmearInRemoved')
        cPermsGrid(M.idSmearInRemoved, :) = [];
        nk = numel(cPermsGrid(:, 1));
    end
    assert(nk == numel(M.Psmear))
    cunit = M.unit(M.isclay);
    for n=1:nk
        C = diag([cPermsGrid(n, 1), cPermsGrid(n, 2)]);
        kmat = transformTensor(C, T, 'posDef');
        k.cSmears(n, :) = [kmat(1,1), kmat(1,2), kmat(2,2)];
        isSmearInUnit = all([whichUnit == cunit(n), isSmear], 2);
        permG(isSmearInUnit,:) = repmat(k.cSmears(n, :), [sum(isSmearInUnit), 1]);
    end
else
    C = diag([clay.perm(1), clay.perm(2)]);
    S = diag([sand.perm(1), sand.perm(2)]);
    kmats = transformTensor({C,S}, T, 'posDef');
    k.vals = [kmats{1}(1,1), kmats{1}(1,2), kmats{1}(2,2);
        kmats{2}(1,1), kmats{2}(1,2), kmats{2}(2,2)];
    permG = repmat(k.vals(2,:), [G.cells.num, 1]);               % sand
    permG(isSmear, :) = repmat(k.vals(1,:), [sum(isSmear), 1]);  % clay
end


end