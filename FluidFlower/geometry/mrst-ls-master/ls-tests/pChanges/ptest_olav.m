%% CO2 injection through single well in a rectangular prism
%
% Lluis Salo, Sept 2020, MIT
% 
% Summary
% Check p changes at the beginning of injection and right after injection
% end.
%
%
clc, clear, close all force
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
%% Initial options
saveData     = 0;
vaporizedOil = 0;               % gas dissolution in brine always


%% Create grid
% 
Lx = 5000;
Ly = 5000;
Lz = 1500;
nx = 50;
ny = 50;
nz = 30;
G  = cartGrid([nx, ny, nz], [Lx Ly Lz]);

% Visualize grid
%figure(1)
%plotGrid(G, 'FaceColor',[.8 .8 .8]); view(25,40), axis tight

% Bring top to 1km depth
G.nodes.coords(:,3) = G.nodes.coords(:,3);

% Compute centroids, etc
G = computeGeometry(G);


%% Set up simulation
% Gravity
gravity reset on

% Rock
idSeal = G.cells.centroids(:,3) < 1000;
rock.poro = 0.25*ones(G.cells.num, 1);
rock.poro(idSeal) = 0.1;
rock.perm = 100*milli*darcy*ones(G.cells.num, 1);
rock.perm(idSeal) = 1e-3*milli*darcy;
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.rocknum    = ones(G.cells.num, 1);

% Fluid
if vaporizedOil == 1
    fn = 'ls-proj/radialCO2storage/co2water40CwetGas.DATA';    
else
    fn = 'co2water40C.DATA'
%     fn = 'ls-tests/pChanges/co2water40C.DATA'; 
    %fn = 'ls-proj/radialCO2storage/co2water50C.DATA';
end
deck  = convertDeckUnits(readEclipseDeck(fn));
cP    = deck.PROPS.ROCK(:, 2);
fluid = initDeckADIFluid(deck);

% Model
% fluid = assignPvMult(fluid, cP, rock.regions.rocknum);
if vaporizedOil == 1
    model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'vapoil', true, 'water', false);
else
    model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
end

% Set gravity and water col
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
p_r = g*rho_wr*50*meter;  % add 50m water column

% Initialize
z = G.cells.centroids(:,3);
[z_0, z_max] = deal(min(z), max(z));
if 0
    equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
    p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
    s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
    rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
    rv0 = 0;                                 % dry gas
    state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
else
    reg = getInitializationRegionsBlackOil(model, 0, 'datum_depth', z_0, 'datum_pressure', p_r, 'rs', @(varargin) 0);
    state0 = initStateBlackOilAD(model, reg);
end
% Wells and timesteps (well perforation thickness = 50m, i.e. all depth)
pdist2 = @(x, y) sqrt(sum((y - x).^2, 2));
dist = pdist2(G.cells.centroids, [2500 2500 1350]);
% dist = sqrt(sum((G.cells.centroids - [2500, 2500, 1350]).^2, 2));
[~, perforatedCell] = min(dist,[],1);

injTime = 10*year;
injRate = 5*10^5*meter^3/day;                                               % [Sm^3/s]
rhoInj  = fluid.rhoGS;
mRate   = injRate*rhoInj*year/10^9;                                         % Mt/y
W = addWell([ ], G, rock, perforatedCell, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injRate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', G.cells.centroids(perforatedCell, G.griddim), 'radius', 0.1);    
        
simTime     = 30*year;
reportTimes = [[1, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
               240, 300, 365, 456.25, 547.5, 638.75, 730]*day, ...
               [2.25:.25:5 5.5:0.5:12, 13:20 22:30]*year];
timesteps   = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')


% model.minimumPressure = min(p0);

% Acceleration
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'deferredAssembly', true);
model = model.validateModel();
nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, 'targetChangeAbs', 0.25);
nls.timeStepSelector = stepSel;
nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');
%nls.enforceResidualDecrease = 1;
%nls.useLinesearch = 1;
%nls = [];

% Changes to model
model.operators.p0          = state0.pressure;                              % Needed for MyPvMult
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));

% BCs - Laterally infinite, so we use volume multipliers at the boundary
bc = [];
east = G.cells.centroids(:,1) < Lx/nx;
west = G.cells.centroids(:,1) > Lx-(Lx/nx);
north = G.cells.centroids(:,2) < Ly/ny;
south = G.cells.centroids(:,2) > Ly-(Ly/ny);
cellsExt = any([east, west, north, south], 2);
% check: figure(1); hold on; plotCellData(G, rock.poro, cellsExt)
model.operators.pv(cellsExt) = model.operators.pv(cellsExt)*1e3;

%% Property plots
latx  = {'Interpreter','latex'};
clrs  = [65,2,0; 164,30,53; 125,52,162; 49,148,206; 220,220,220]./255;
%Sbr   = linspace(f.krPts.og(1, 2), f.krPts.og(1, 3), 50);
%Sbs   = linspace(f.krPts.og(2, 2), f.krPts.og(2, 3), 50);
 
% Density, viscosity
%  p = linspace(50, 250, 50)*barsa;
%  rhoG = fluid.bG(p).*fluid.rhoGS;
%  muG  = fluid.muG(p);
%  rs   = fluid.rsSat(p);
%  rhoB = fluid.bO(p, rs, true(numel(p), 1)).*(rs.*fluid.rhoGS + fluid.rhoOS);
%  muB  = fluid.muO(p, rs, true(numel(p), 1));
%  rhoB = fluid.bO(p, zeros(numel(p), 1), false(numel(p), 1)).*fluid.rhoOS;
%  muB  = fluid.muO(p, zeros(numel(p), 1), false(numel(p), 1));
%  
%  h = figure(35);
%  subplot(1,2,1)
%  plot(p/barsa, rhoB, '-', 'color', clrs(4,:), 'linewidth', 1.5);
%  hold on; plot(p/barsa, rhoG, '-', 'color', clrs(3,:), 'linewidth', 1.5);
%  hold off,  xlim([min(p) max(p)]./barsa), grid on; ax = gca;
%  xlabel('p [bar]','fontsize', 14, latx{:})
%  ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
%  set(ax,'Xtick',50:25:250);
%  set(ax,'Ytick',0:100:1100); ylim([0 1100]);
%  legend({'Brine', 'CO$_2$'}, latx{:}, 'fontSize', 10)
%  title('Density')
%  
%  subplot(1,2,2)
%  semilogy(p/barsa, muB/(centi*poise), 'color', clrs(4,:), 'LineStyle', '-', 'linewidth', 1.5);
%  hold on; semilogy(p/barsa, muG/(centi*poise), '-', 'color', clrs(3,:), 'linewidth', 1.5);
%  hold off;  xlim([min(p) max(p)]./barsa)
%  ylabel('$\mu$ [cP]','fontsize', 14, latx{:}); ax = gca;
%  set(ax,'Ytick',[0 0.15 0.3 0.45 0.6]); 
%  ylim([0.01 0.9])
%  set(ax,'Xtick',50:25:250);
%  set(ax,'Ytick', [0.01:0.01:0.1 0.2:.1:.9]);
%  set(h, 'Position', [600, 600, 380, 320])
%  grid on
%  title('Viscosity')

% Schedule
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules 
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % nr of wells must be the same for the entire simulation
schedule.control(2).W.val = 0;                                              % field 2 rate 0 (after injection)
schedule.control(2).W.status = false;                                              % field 2 rate 0 (after injection)

schedule.step.control(cumsum(schedule.step.val) > injTime) = 2;                % set timesteps after injTime to use Well field 2

%% Run simulation
% we run a parallel simulation with N threads
%pause
N = 4;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                   % Specify threads manually
%schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);                      % Simple schedule, this injects for the total simTime
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'NonLinearSolver', nls, ...
                                                'Verbose', true);

%% Plots
model = model.validateModel();
bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    bhp(n) = wellSols{n}.bhp;
    states{n}.FlowProps.Viscosity = model.getProp(states{n}, 'Viscosity');
    states{n}.FlowProps.phaseDensity = model.getProp(states{n}, 'Density');
end

% General
% cmap = flipud(cmap_agu());
cmap = colormap();
figure(3)
colormap(cmap)
plotToolbar(G, states, 'FaceAlpha', 0.85);
view(53, 63), camproj perspective; axis equal tight
plotWell(G, W, 'color', 'k', 'height', 5);
c = colorbar; caxis([0 1])

%Cross section
idCenter = G.cells.centroids(:,2) == G.cells.centroids(W.cells,2);
figure(4)
colormap(cmap)
plotToolbar(G, states, idCenter,'FaceAlpha', 0.85);
view(-18, 27), camproj perspective; axis equal tight
plotWell(G, W, 'color', 'k', 'height', 5);
c = colorbar; caxis([0 1])
%%
% BHP
p0 = state0.pressure;
figure(5)
subplot(2,1,1)
plot([0 reportTimes(1:end)/year], [p0(W.cells)./10^5 bhp(1:end)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
%plot([0 reportTimes/year], [105 bhp./10^5], '.-k')
xlabel('time [y]')
ylabel('p [bar]')
%ylim([100 171])
%yticks([100 110 120 130 140 150 160 170])
grid on
title('BHP during injection period (50y)')

subplot(2,1,2)
plot([0 reportTimes(1:16)/day], [p0(W.cells)./10^5 bhp(1:16)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
xlabel('time [days]')
xticks([0 7 30 60 90 120 150 180 210 240 270 300 330 365])
yticks([100 110 120 130 140 150 160 165])
%ylim([100 165])
hold on
%text(0.75,188,'6h')
ylabel('p [bar]')
grid on
xlim([0 365])
title('Detail during first year')

% Monitoring cell (right above injector cell)
dist = pdist2(G.cells.centroids, [2500 2500 G.cells.centroids(W.cells,3)-(Lz/nz)]);
[~, monitoringCell] = min(dist,[],1);

dpCell = zeros(numel(states), numel(monitoringCell));
sCell = zeros(numel(states), numel(monitoringCell));
for n=1:numel(states)
    dpCell(n, :) = states{n}.dp(monitoringCell);
    sCell(n, :) = states{n}.s(monitoringCell, 2);
end

latx = {'Interpreter', 'latex'};
f16 = figure(16);
left_color = [0 0 0];
right_color = [0 0 1];
set(f16,'defaultAxesColorOrder',[left_color; right_color]);

fill([0 0 10 10 0], [-2 15 15 -2 -2], [0.85 0.85 0.85], 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dpCell(1:end,1)./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 3)
%annotation(f16,'doublearrow', [0.130697210391006 0.157422969187675], [0.282333333333333 0.280952380952381]);
%text(5, 0, 'Inj.', latx{:}, 'fontSize', 12)
hold off
xlabel('$t$ [y]', latx{:})
xticks([0 5 10 20 30])
ylim([-2 15])
ylabel('$\Delta p$ [bar]', latx{:})
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sCell(:, 1)], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
ylabel('$S_\mathrm{g}$ [-]', latx{:})
ylim([0 0.7])
grid on
%legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
%ylim([0 100])
title('Monitoring cell', latx{:}, 'fontSize', 12)

