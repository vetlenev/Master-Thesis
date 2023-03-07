%% Compute permeability of a fault with clay smears
%
% In the configuration used here, FW on the left and HW is on the right.
%
% The smears are consistent with (1) the number of sand and clay layers at
% each side of the fault; (2) the thickness of these layers; (3) fault 
% displacement and thickness; (4) the residual friction angle of the sands
% and clays; (5) the SSFc; and (6) the maximum length of a smear segment in
% each subdomain (to control whether smears are broken or not).
% Based on these inputs (1) to (6), the following are calculated: (I.) the
% probability of encountering a cell with smear in each smear subdomain;
% (II.) the length of each smear segment within each subdomain; and (III.)
% the thickness of each smear (and accordingly, the thickness of each
% subdomain).
%
%__________________________________________________________________________
clc, clear, close all
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils


%% 1. Initialize
% Input independent variables
% Note: fw and hw layers' thickness (T) is the true layer thickness
%       to agree with Egholm et a. (2008). Relation to the apparent thickness
%       is Tap = T/(cosBl*sinB) where Bl is layer dip angle and B is
%       fault dip angle.
f    = struct('T', 5, 'dip', 46.66);            % fThickness not from grid
sand = struct('perm', [100, 30]*milli*darcy, 'poro', 0.3, 'phi', 35);
clay = struct('perm', [0.01, 1e-4]*milli*darcy, 'poro', 0.1, 'phi', 10, 'SSFc', 5);
fw   = struct('T', repelem(25,4), 'isclay', logical([0 1 0 1]), 'Bl', 0);
hw   = struct('T', 100, 'isclay', false, 'Bl', 0);
g    = struct('resL', 1, 'resT', 0.1);          % grid resolution [m]
smear.inL  = 'Lsmear';                          % Smear segment length [m] 
smear.modeledAs = 'objects';                    % algorithm to place smears
smear.tolerance = 0.01;                         % tolerance in P matching  
displayFigures = [0 1 1];                       % [grid, smear location, perm] 

%% 2. Create and populate grid
% Cartesian grid in MRST has indexing that starts at (0, 0) coordinate
% and runs faster on the X axis. Mapping grid and values for properties 
% are stored in structure M.

% 2.1 Compute initial values
[fw, hw, f, clay, smear] = smearOrganizeIn(f, sand, clay, fw, hw, smear);

% 2.2 Create and plot grid
ydim = 5;
nels = max([round(f.T/g.resT), round(f.D/g.resL)]);
G = computeGeometry(cartGrid([nels,10,nels], [f.T ydim f.D]));
G.xzFaceDim = [f.T/nels, f.D/nels];

% 2.3 Mapping matrix structure
M = smearMap(G, fw, hw, f, clay, smear);

% 2.4 Smears (final 1s)
[M, smear] = smearContinuity(fw, f, M, clay, smear);

if any(M.Psmear(M.isclayIn) < 1) % place objects (segments) of Lsmear, and match P(smear)
    assert(strcmp(smear.modeledAs, 'objects'));
    M = smearObjects(M, smear, f, G);
else
    disp('Continuous smears: P(smear) = 1')
end

% 2.5 To compromise spdiags truncation rules with smears chopped parallel
% to f.T (what we want) and not to f.L, the last step is to transpose. 
% Note that this needs to be done since the diagonal indices (diagIds) used 
% in spdiags during M.vals and Mtest assignment are flipped. Mtest is 
% directly transposed above, so that P is computed and matched with the 
% same exact configuration as in the definitive mapping matrix (M.vals).
M.vals = transpose(M.vals);


% 2.5 Map diagonals to Grid indexing: Grid indexing starts at bottom left, 
%     columns (x) move faster. Matrix starts counting at top left and rows 
%     (y) move faster. This is for 2D grid. For 3D grid, the x-z plane 
%     starts at top left (z positive downwards), so we just need to transpose.

% figure(9); spy(sparse(M.vals))       % check placement of smears in M
idGridFrontXZ = reshape(transpose(M.vals), G.cartDims(1)*G.cartDims(3), 1);
MunitsFrontXZ = reshape(transpose(M.units), G.cartDims(1)*G.cartDims(3), 1);
xyNum = G.cartDims(1)*G.cartDims(2);
idGrid.ids = zeros(G.cells.num,1);
idGrid.Munits3D = zeros(G.cells.num,1);
for n=1:G.cartDims(3)
    ids = (n-1)*xyNum+1:n*xyNum;
    idsXZ = (n-1)*G.cartDims(1)+1:n*G.cartDims(1);
    idGrid.ids(ids) = repmat(idGridFrontXZ(idsXZ),G.cartDims(2),1);
    idGrid.Munits3D(ids) = repmat(MunitsFrontXZ(idsXZ),G.cartDims(2),1);
end
idGrid.ids = logical(idGrid.ids);

% 2.6 Assign permeability and porosity. Note that we receive FZ sand/clay poro
%     and permeability, the latter given perpendicular and parallel to the
%     to the bedding. This is then directly assumed to be the permeability
%     along/across the clay/sand smears (principal directions). This diagonal
%     permeability tensor needs to be rotated to the fault local axes,
%     since clay/sand smears are not parallel to the fault.
% Poro
rock.poro = ones(G.cells.num, 1)*sand.poro;
rock.poro(idGrid.ids) = clay.poro;

% Perm
T = [cosd(f.alpha) 0 sind(f.alpha); 0 1 0; -sind(f.alpha) 0 cosd(f.alpha)];
C = diag([clay.perm(2), clay.perm(1), clay.perm(1)]);
S = diag([sand.perm(2), sand.perm(1), sand.perm(1)]);
k.mats = transformTensor({C,S}, T, 'posDef');
k.vals = [k.mats{1}(1,1), k.mats{1}(1,2), k.mats{1}(1,3) ...
          k.mats{1}(2,2), k.mats{1}(2,3), k.mats{1}(3,3);
          k.mats{2}(1,1), k.mats{2}(1,2), k.mats{2}(1,3) ...
          k.mats{2}(2,2), k.mats{2}(2,3), k.mats{2}(3,3)];
rock.perm = repmat(k.vals(2,:), [G.cells.num, 1]);                   % sand
rock.perm(idGrid.ids,:) = repmat(k.vals(1,:), [sum(idGrid.ids), 1]); % clay


%% 3. Compute equivalent/upscaled permeability and porostiy
% Porosity (additive)
% Here, for 1 cell in upscaled grid and regular cells we could just do 
% Poro = (sum(sum(M.vals))*clay.poro +  ...
%        (G.cartDims(1)^2 - sum(sum(M.vals)))*sand.poro) / G.cartDims(1)^2;
p = partitionCartGrid(G.cartDims, [1 1 1]);
CG = generateCoarseGrid(G, p);
Poro = accumarray(p, rock.poro)./accumarray(p,1);
%plotGrid(CG,1,'FaceColor','none','EdgeColor','k','LineWidth',1.5);
%set(gcf,'PaperPositionMode','auto'); axis equal;

% Permeability (NOT additive)

% Plots
smearFigs(displayFigures, G, M, f, fw, hw, idGrid, rock)