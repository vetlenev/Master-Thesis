%% Load GRDECL geometry
grdecl = readGRDECL('deck/CSP11A.GRDECL');
units = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, units);

G = processGRDECL(grdecl);
G = computeGeometry(G);

%% Read data of SPE 11. Extract rock and fluid properties.
%rock = grdecl2Rock(grdecl, G.cells.indexMap);
deck = readEclipseDeck('deck/CSP11A.DATA');
deck = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);
%fluid = assignRelPerm(fluid);

permX = deck.GRID.PERMX(G.cells.indexMap); % only choose active cells
permY = deck.GRID.PERMY(G.cells.indexMap);
permZ = deck.GRID.PERMZ(G.cells.indexMap);
poro = deck.GRID.PORO(G.cells.indexMap);
pv = deck.GRID.PORV(G.cells.indexMap);

perm = [permX, permY, permZ];
rock = makeRock(G, perm, poro);

%% Define model
model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);

%% Initial state
z = G.cells.centroids(:,3);
state0 = initResSol(G, fluid.rhoWS*norm(gravity)*z, [1,0]);
state0.sGmax = state0.s(:,2);

%% Operational schedule
bc = addBC([], 'Top', 'pressure', 1.1*10^5*Pascal, 'sat', [1,0]);

tot_time = 5*day; % 400*year
inj_start_W1 = 0;
inj_stop_W1 = 5*hour;
inj_start_W2 = 2.5*hour;
inj_stop_W2 = 5*hour;

inj_rate_W1 = 1.6*10^(-7) * kilogram/second;
inj_rate_W2 = 1.6*10^(-7) * kilogram/second;

min_x = @(x) min(abs(G.cells.centroids(:,1) - x));
min_z = @(z) min(abs(G.cells.centroids(:,3) - z));
idx_x = @(x) find(G.cells.centroids(:,1) - x == min_x(x));
idx_z = @(z) find(G.cells.centroids(:,3) - z == min_z(z));

min_x1 = min_x(0.9); idx_x1 = idx_x(0.9);
min_x2 = min_x(1.7); idx_x2 = idx_x(1.7);
min_z1 = min_z(0.3); idx_z1 = idx_z(0.3);
min_z2 = min_z(0.7); idx_z2 = idx_z(0.7);

W1_cell = intersect(idx_x1, idx_z1);
[ii, ~, kk] = gridLogicalIndices(G);
W1_x = ii(W1_cell);
W1_z = kk(W1_cell);

W1 = verticalWell([], G, rock, W1_x, 1, W1_z, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W1, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1]);

W2_cell = intersect(idx_x2, idx_z2);
W2_x = ii(W2_cell);
W2_z = kk(W2_cell);

W2 = verticalWell([], G, rock, W2_x, 1, W2_z, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W2, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1]);


%% Plot grid
figure()
plotGrid(G)
plotCellData(G, poro)
view(0, 0)
axis equal tight
colorbar('southoutside')
