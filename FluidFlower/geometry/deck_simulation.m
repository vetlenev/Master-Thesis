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

permX = deck.GRID.PERMX(G.cells.indexMap);% .* 1e-5; % only choose active cells
permY = deck.GRID.PERMY(G.cells.indexMap);% .* 1e-5;
permZ = deck.GRID.PERMZ(G.cells.indexMap);% .* 1e-5;

poro = deck.GRID.PORO(G.cells.indexMap);
pv = deck.GRID.PORV(G.cells.indexMap);

perm = [permX, permY, permZ];
rock = makeRock(G, perm, poro);

facies = deck.REGIONS.SATNUM(G.cells.indexMap);
rock.regions.saturation = facies;

%% Rename fluid field
[fluid.krO] = fluid.krOG;
fluid = rmfield(fluid, 'krOG');
[fluid.krPts.o] = fluid.krPts.og;
fluid.krPts = rmfield(fluid.krPts, 'og');

b = 1; % unit formation volume factor
fluid.bO = @(p, varargin) 1 + 0.*p; %b*constantReciprocalFVF(p, varargin{:});
fluid.bG = @(p, varargin) 1 + 0.*p; %b*constantReciprocalFVF(p, varargin{:});
fluid.muO = @(p, varargin) 3e-5*Pascal*second + 0.*p;
fluid.muG = @(p, varargin) 8e-4*Pascal*second + 0.*p;

%% Define model
%model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);
model = GenericBlackOilModel(G, rock, fluid);
model.water = false;
%model.OutputStateFunctions{end+1} = {'SurfaceDensity'};
model = validateModel(model);

%model.FlowPropertyFunctions.RelativePermeability = RelativePermeabilityHys(model);

%% Initial state
z = G.cells.centroids(:,3);
state0 = initResSol(G, fluid.rhoWS*norm(gravity)*z, [1,0]);
state0.sGmax = state0.s(:,2);

%% Operational schedule
bf = boundaryFaces(G);
bfz = G.faces.centroids(bf, 3);
top_faces = bf(bfz == min(bfz));
%bc = pside([], G, 'Top', 1.1*10^5*Pascal, 'sat', [1,0]);
bc = addBC([], top_faces, 'pressure', 1.1*10^5*Pascal, 'sat', [1,0]);

ds = deck.SCHEDULE;
dt = ds.step.val;
t = cumsum(dt);
tot_time = t(end);

tot_time = 5*day;
t_control2 = 2.5*hour;% 2.5*hour;
t_control3 = 5*hour; %5*hour;
dt1 = rampupTimesteps(t_control2, t_control2/70, 10);
dt2 = rampupTimesteps(t_control3-t_control2, (t_control3-t_control2)/120, 12);
dt3 = rampupTimesteps(tot_time-t_control3, (tot_time-t_control3)/120, 12);
dt = [dt1; dt2; dt3];
t = cumsum(dt);

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

W = verticalWell([], G, rock, W1_x, 1, W1_z, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W1, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1], ...
                'status', true); % starts injecting

W2_cell = intersect(idx_x2, idx_z2);
W2_x = ii(W2_cell);
W2_z = kk(W2_cell);

W = verticalWell(W, G, rock, W2_x, 1, W2_z, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W2, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1], ...
                'status', false); % starts off

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

schedule.control(2:3) = schedule.control(1);
schedule.control(2).W(2).status = true;
schedule.control(3).W(1).status = false;
schedule.control(3).W(2).status = false;

schedule.step.control(t <= t_control2) = 1;
schedule.step.control(t > t_control2 & t <= t_control3) = 2;
schedule.step.control(t > t_control3) = 3;

%% Run simulation
nls = NonLinearSolver('maxIterations', 70);
%[ws_dummy, states_dummy, report_dummy] = simulateScheduleAD(state0, model, schedule);
[ws, states, report] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

%% Plot grid
N = G.faces.neighbors;
mask = all(N > 0, 2);
pnn = zeros(G.faces.num, 2);
pnn(mask,:) = permZ(N(mask,:));
perm_diff = pnn(:,1) - pnn(:,2);
perm_change = find(perm_diff);
faces_facie = false(G.faces.num, 1);
faces_facie(perm_change) = true;

figure(5)
plotOutlinedGrid(G,W,bc,faces_facie)
plotCellData(G, permZ)
view(0, 0)
axis equal tight
colorbar('southoutside')

%% Plot results
c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');
rootdir = strrep(ROOTDIR, '\', '/');

for i=1:10:numel(states)
   f2 = figure(2); clf
   sf = states{i}.s(:,2);
   plotCellData(G, sf, 'edgec', 'none')
   plotFaces(G, faces_facie, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
   view(0, 0); colormap(cc); caxis([0,1]);
   colorbar('location', 'southoutside');
   axis tight on    
   pause(0.2)
   title({'Fine saturation', sprintf('hours: %.1f',t(i)/hour)}) 

   saveas(f2, sprintf(strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/output/sat_%d'), i), 'png')
end
