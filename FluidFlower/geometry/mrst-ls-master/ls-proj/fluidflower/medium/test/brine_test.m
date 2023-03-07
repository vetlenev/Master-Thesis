%% Brine test to equilibrium (blackoil module)
%--------------------------------------------------------------------------
% 3 saturation regions, kr hysteresis, fault Pc scaled
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
mrstVerbose on


%% Options
hyster = 'off';
topboundary = 'straight';
bctype = 'pvm';                                                % constant pressure (cp) or pore volume multiplier (pvm)
depth = 'surface'; % 'surface' or '1km'


%% Save data
folderName = ['brine_test_' depth];
topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\tests\';

%% Plots
fig3D = @() figure('Position', [0, 0, 1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );


%% Mesh
% Load simpleExtrudedMesh
G = cartGrid([1, 100, 100], [0.01, 0.1, 0.1]);
if strcmp(depth, '1km')
    G.nodes.coords(:,3) = G.nodes.coords(:,3)+1000;
end
G = computeGeometry(G);
idleft = G.cells.centroids(:, 2) < 0.05;

% Plot
% fig3D(); plotGrid(G); plotGrid(G, idleft, 'edgecolor', 'r'); 
% view([90, 0])


%% Rock
% Layers follow uGr (bot to top, left to right for faulted, fault at the end)
nc = G.cells.num;
rock.poro = nan(G.cells.num, 1);
rock.perm = nan(G.cells.num, 1);    % isotropic
poroMean = [0.35];
r = 0.03; % half range
rock.poro = (poroMean - r) + rand(nc, 1)*(2*r);
permMean = 5000;   % [mD]
logr = 0.5;
rock.perm = 10.^(rand(nc, 1)*(2*logr) + (log10(permMean) - logr));
rock.perm = rock.perm*(milli*darcy); % SI units

% Plot permeability
% fig3D(); plotCellData(G, log10(rock.perm/(milli*darcy)), 'edgealpha', 0.2)
% setAxProps(gca), colormap(turbo), c = colorbar; %caxis([3 4])
% axis equal off; view([90 0]), %ylim([0.3 0.6]); zlim([0.42 0.47])


%% Fluids
% 3 Saturation regions
rock.regions.saturation = ones(G.cells.num, 1);                              % reservoir units

% 1 Regions for rock compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% load deck and initialize fluid
if strcmp(depth, 'surface')
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium');
    fn  = fullfile(pth, 'krHyst_fPcFault_co2wat_test.DATA');
else
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/test/brine_test_props');
    fn  = fullfile(pth, 'brinetest_1km.DATA');
end
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);

% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Hysteresis
% if isfield(fluid, 'krHyst') && fluid.krHyst == 1
%    numReg = max(rock.regions.saturation);
%    fluid.krHyst = 1 + numReg;                % porous rock is first imb region (4), fault is third (6)
%    rock.regions.imbibition = rock.regions.saturation + numReg;
%    fluid = addScanKr(fluid, rock.regions.imbibition);
% end

% Plot density, viscosity
% f = fluid;
% p = linspace(90, 110, 50)*barsa;
% rs   = f.rsSat(p);
% rhoB = f.bO(p', rs', true(numel(p), 1)).*(rs.*f.rhoGS + f.rhoOS)';
% rhoB_noCO2 = f.bO(p', zeros(numel(p), 1), false(numel(p), 1)).*f.rhoOS;
% muB  = f.muO(p, rs, true(numel(p), 1));
% plot(p, rhoB, '-b'); hold on; plot(p, rhoB_noCO2, '--b')
% 
% p2 = 100*barsa*ones(50,1);
% rs = linspace(0,f.rsSat(p2(1)),50)';
% rhoB2 = [f.bO(p2(1:end-1),rs(1:end-1),false(49,1)); ...
%          f.bO(p2(1),rs(end),true)].*(rs.*f.rhoGS + f.rhoOS);
% plot(rs, rhoB2, '-k'); xlabel('rs'), ylabel('\rho [kg/m^3]')

%% Initialize
gravity reset on
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
water_column = 30; % m
p_r = 100*barsa + g*rho_wr*water_column; % p at shallowest z
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
rs0 = zeros(G.cells.num, 1);
rs0(idleft) = fluid.rsSat(p_r);
issat = false(G.cells.num, 1);
issat(idleft) = true;
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rv0 = 0;             
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% %plot
% fig3D(); plotCellData(G, p0, 'edgealpha', 0.2)
% setAxProps(gca), colormap(jet), c = colorbar;
% axis equal off


%% No Wells        
simtime = 5*day;
reportTimes = [1 10:10:24*60 1470:30:7200]*minute;                             
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simtime, 'sum of timesteps must equal simTime')


%% Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'vapoil', false, 'water', false);
                         
% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true);
model.toleranceCNV = 5e-3;
model = model.validateModel();

% nls = getNonLinearSolver(model, 'TimestepStrategy', 'ds', ...
%                          'useCPR', true);
nls = getNonLinearSolver(model, 'TimestepStrategy', 'none', ...
                         'useCPR', true);
stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, ...
                                     'targetChangeAbs', 0.25);
nls.timeStepSelector = stepSel;
if strcmp(bctype, 'pvm')
    nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-3, 'Solver', 'bicgstab');
end
nls.LinearSolver.maxIterations = 50;
%nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'gmres', ...
%                                          'preconditioner', 'amg', ...
%                                          'coarsening', 'smoothed_aggregation');
%nls.useRelaxation = 1;
%nls.relaxationType = 'sor';
%nls.enforceResidualDecrease = 1;
nls.useLinesearch = 1;
nls.maxIterations = 12;
nls.maxTimestepCuts = 10;
nls.acceptanceFactor = 5;
%nls.alwaysUseStabilization=1;
%nls.convergenceIssues=true;
%nls = [];

% Model changes (after acceleration!)
model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
if strcmp(hyster, 'on')
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
    %model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));
end
model.minimumPressure = min(state0.pressure);

%% BCs -- closed
% f = false(G.faces.num, 1);
% f(boundaryFaces(G)) = true;
% f = all([G.faces.centroids(:, 3) < 1e-4, f], 2);
% bc = addBC([], find(f), 'pressure', p_r, 'sat', [1 0]);

%fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r')
%setAxProps(gca), colormap(jet), c = colorbar; axis equal off


%% Schedule
schedule = simpleSchedule(timesteps, 'W', [], 'bc', []);                 % Simple schedule, this would inject for the total simTime


%% Simulation
N = 2;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                           'NonLinearSolver', nls, ...
%                                           'Verbose', true);
 
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  


%% Results
% compute quantitites
model = model.validateModel();
%bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;
    %rho = model.getProp(states{n}, 'ComponentPhaseDensity');
    %states{n}.FlowProps.ComponentPhaseDensity = rho; % 2 components in brine phase.
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    %bhp(n) = wellSols{n}.bhp;
end
states{1}.reg = model.rock.regions.saturation;


% basic overview
%rescells = find(ismember(G_dat.compartID, uGr{2}));
fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); 
plotGrid(G, idleft, 'edgecolor', [0.5 0.5 0.5]);
setAxProps(gca), colormap(turbo), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]) %ylim([0.42 0.48]), zlim([0.40 0.47])