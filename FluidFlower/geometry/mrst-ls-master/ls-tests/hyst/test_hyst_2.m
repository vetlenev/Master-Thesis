%% Black Oil Test: relative permeability hysteresis
%
% 3D mesh with 2 fluids (water and CO2)
%
%
%
clear, close all
mrstModule add ls-tests ls-utils ad-props deckformat mrst-gui ad-core ad-blackoil co2lab linearsolvers
%global BOOSTPATH AMGCLPATH  
%BOOSTPATH = '';
%AMGCLPATH = '';
gravity reset on;

% Mesh
meshpath = 'mrst-dev/mrst-ls/ls-tests/hyst/nodes_coordinates.dat';
vertices = dlmread(meshpath);
x = vertices(:,1); y = vertices(:, 2).*-1;
thick = [linspace(200,20,5) 15 linspace(20,200,5)];
G = makeLayeredGrid(triangleGrid([x, y]), thick);
G.nodes.coords = G.nodes.coords(:,[3 1 2]);
G = computeGeometry(G);

% Rock & regions
rock = makeRock(G, 100*milli*darcy, 0.3); % homogeneous perm
id_seal = G.cells.centroids(:,3) < 360;
rock.perm(id_seal) = 1e-3*milli*darcy;

% Input deck and fluid.
%fn   = 'mrst-2019b/mytests/CASE2_props.DATA';      % 'WG' model, tested.
%fn   = 'mrst-2019b/myprojects/imb_co2brine.DATA';   % 'OG' model
fn   = 'mrst-dev/mrst-ls/ls-proj/Pc_imb_co2brine.DATA';   % 'OG' model, pcHyst


% Regions
if strcmp(fn, 'mrst-dev/mrst-ls/ls-tests/hyst/CASE2_props.DATA')
    rock.regions.saturation = ones(G.cells.num,1);
    rock.regions.imbibition = ones(G.cells.num,1)*2;
else
    rock.regions.saturation = ones(G.cells.num,1);
    rock.regions.saturation(id_seal) = 2;
    rock.regions.imbibition = ones(G.cells.num,1)*3;
    rock.regions.imbibition(id_seal) = 4;
end

% Deck
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.SATNUM = rock.regions.saturation;
deck.REGIONS.IMBNUM = rock.regions.imbibition;
fluid = initDeckADIFluid(deck);
if ~isfield(fluid, 'bO')
    fluid.bO = [];
end
fluid = addScanKr(fluid, rock.regions.imbibition);
fluid = addScanPc(fluid, rock.regions.imbibition);
%fluid = rmfield(fluid, 'krHyst');                       % disable hysteresis
fluid.krHyst = 4;
fluid.pcHyst = 4;
%pause

% Initialize reservoir
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
if isfield(fluid, 'muO') %'OG'
    rho_wr = fluid.rhoOS*kilogram/meter^3;
elseif isfield(fluid, 'bW') %'WG'
    rho_wr = fluid.rhoWS;
end
p_r = g*rho_wr*z_0; 
if isfield(fluid, 'muO')
    equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
    rs0    = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
elseif isfield(fluid, 'bW')
    equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
    rs0    = 0;
end
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% Wells
simTime = 100*year;

numwell = 1;                                  % 4 horizontal injectors
injmass = 0.5*10^9*(simTime/year);           % 0.5 Mt/year of gas
injrate = injmass/(fluid.rhoGS(1)*simTime*numwell);
%wellLoc = [25,2300,500; 25,2360,500; 25,2420,500; 25,2480,500];
%wellLoc = [25,2360,500]; 1.3375e+03
wellLoc = [sum(thick)/2, 2360, 500]; 
dirw = 'z';
distx = pdist2(G.cells.centroids(:,1), wellLoc(:,1));
disty = pdist2(G.cells.centroids(:,2), wellLoc(:,2));
distz = pdist2(G.cells.centroids(:,3), wellLoc(:,3));
dist = sqrt(distx.^2 + disty.^2 + distz.^2);
[~, wellInx] = min(dist,[],1);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', dirw, ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1]); % order is always 'WOG' ('WG' in GenericBlackOilModel)
        
% BC
idb = G.faces.normals(:, 3) == 0;
cellsn0 = G.faces.neighbors(idb, :);
idb1 = cellsn0(:,2) == 0; idb2 = cellsn0(:,1) == 0;
cellsb = unique([cellsn0(idb1, 1); cellsn0(idb2, 2)]);

% Changes to model
if isfield(fluid, 'muO')
    model = GenericBlackOilModel(G, rock, fluid, 'water', false, ...
                                 'disgas', true);
elseif isfield(fluid, 'bW')
    model = GenericBlackOilModel(G, rock, fluid, 'oil', false);
end
model = model.validateModel();
ncomp  = model.getNumberOfComponents();
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
solver = AMGCL_CPRSolverAD('tolerance', 1e-5, 'block_size', ncomp, ...
                           'useSYMRCMOrdering', true, 'coarsening', ...
                           'aggregation', 'relaxation', 'ilu0');
nls = NonLinearSolver('LinearSolver', solver);

%nls = NonLinearSolver('tolerance', 1e-5, 'maxIterations', 10, 'maxTimestepCuts', 8);  % converges with current times setup.
%nls = [];

model = model.validateModel();
model.operators.pv(cellsb) = model.operators.pv(cellsb)*1e3;
%model.FlowPropertyFunctions = MyFlowPropertyFunctions(model);
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', BlackOilCapillaryPressure(model));
%model.OutputStateFunctions = {'ComponentTotalMass','RelativePermeability','CapillaryPressure'};
model = model.validateModel(); % Set up the state function groups
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
model.OutputStateFunctions = {'ComponentTotalMass','RelativePermeability', 'Density', ...
                              'Mobility', 'PhasePressures', 'PoreVolume', 'ShrinkageFactors', 'Viscosity'};

% Timesteps & Schedule
%times = [0.1:.1:3 3.2:.2:4 4.5:.5:10 10.1:0.1:12 12.5:.5:20 21:50 52:2:100]*year;
%times = [0.1:.1:4 4.2:.2:6 6.5:.5:10 10.1:0.1:12 12.5:.5:20 21:50 52:2:100]*year;
times = [0.001 0.004 0.01 0.02 0.04 0.07 0.1:.1:2 2.2:.25:5 ...
         5.5:.5:10 10.1:.1:11 11.25:.25:13 13.5:.5:20 21:50 52:2:100]*year;
timesteps = [times(1) diff(times)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', []);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % field 2 empty well (after injection)
schedule.control(2).W.val = 0;
schedule.step.control(cumsum(schedule.step.val) > 10*year) = 2;             % inject 10 years only

% Preliminary plots
% figure;
% plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
%              'FaceAlpha', 0.95, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
% plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
% title('Permeability (mD)')
% axis tight, view(-40, 70), c = colorbar('SouthOutside'); 
% caxis([0.001 100]); c.Ruler.Scale = 'log';

% Simulation
%nls = NonLinearSolver('useLinesearch', true);
maxNumCompThreads(1)
solver.amgcl_setup.nthreads = 1;
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'nonlinearsolver', nls, ...
                                                'Verbose', true);

% Results
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.pc = pc{1,2};
end

figure(10)
plotToolbar(G, states); view(-40, 70); set(gca,'YDir','reverse')

figure(11)
allCells = 1:G.cells.num;
plotToolbar(G, states, allCells(~id_seal)); view(-40, 70); 
set(gca,'YDir','reverse'); caxis([0 1]); colormap jet; axis equal

figure(12)
allCells = 1:G.cells.num;
plotToolbar(G, states); view(-40, 70); 
set(gca,'YDir','reverse'); caxis([0 1]); colormap jet; axis equal

figure(13)
cellsAtCenter = (1:G.layerSize) + G.layerSize*((G.numLayers-1)/2);
plotCellData(G, states{end}.s(:, 2), cellsAtCenter, ...
             'FaceAlpha', 0.9, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
set(gca,'YDir','reverse')
title('CO_2 saturation, t=500y'); colormap(jet)
axis tight, view(90, 0), c = colorbar('SouthOutside'); caxis([0 0.7])
axis equal

%% kr Hyst plots 
ns = 50;
s = linspace(0, 0.7, ns);
kr = zeros(ns, ns);
kr_d = zeros(ns, ns);
for i = 1:ns
    for j = 1:ns
        sg = s(i);
        sg = initVariablesADI(sg);
        sgmax = s(j);
        sgmax = max(sgmax, sg);
       
        kg = fluid.krGi{end}(sg, sgmax);
        kr(i, j) = value(kg);
        kr_d(i, j) = kg.jac{1}(1, 1);
    end
end
%%
figure(14)
subplot(1,2,1)
plot(s, kr(:,10), 'color', [0.5 0.5 0.5], 'linewidth', 4)
hold on
plot(s, kr(:,35), 'b', 'linewidth', 1.5)
plot(s, kr, '.-');
title('kr_g')

subplot(1,2,2)
stairs(s, kr_d(:,10), 'color', [0.5 0.5 0.5], 'linewidth', 4)
hold on
stairs(s, kr_d(:,35), 'b', 'linewidth', 1.5)
stairs(s, kr_d);
title('d kr_g / d sg')