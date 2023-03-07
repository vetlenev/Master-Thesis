%% CO2 injection through single well in radial domain
%
% Lluis Salo, July 2020, MIT
% 
% Background:
% Request from ExxonMobil to compare with analytical
% solution by Azizi and Cinar (2013; based on the earlier work by 
% Noh, Bryant and Lake, 2007) and CMG's solution.
%
% Summary:
% Run simulation of CO2 injection through a single central injector, during
% 50 years, in a 50km-radius radial domain, initially saturated with water. 
% Dissolution of CO2 in water active.
%
%
clc, clear, close all force
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
%% Initial options
saveData = 0;
vaporizedOil = 1;               % gas dissolution in brine always


%% Create radial grid 
% (see fig. 3.25 in pg 86 and sect 3.5.2 in MRST book).
% Scale grid to 50km radius and 50m thickness, and move top depth to 1000m.
P = []; 
for r = exp([-9.25:.1:-5 -4.95:0.05:-3 -2.975:0.025:0])             
    [x,y,z] = cylinder(r,125); 
    P = [P [x(1,:); y(1,:)]];
end
P = unique([P'; 0 0],'rows'); 
P = P.*50*kilo*meter;
G = makeLayeredGrid(pebi(triangleGrid(P)), repelem(12.5,4)); 
%k  = G.nodes.coords(:, 3) > 0;
%G.nodes.coords(k, 3) = 50;
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000;

% Visualize grid
figure(1)
plotGrid(G,'FaceColor',[.8 .8 .8]); view(25,40), axis equal tight

figure(2)
subplot(1,2,1)
plotGrid(G,'FaceColor',[.8 .8 .8]); view(0,90), axis equal tight
xlim([-10000 10000]), ylim([-10000 10000])

subplot(1,2,2)
plotGrid(G,'FaceColor',[.8 .8 .8]); view(0,90), axis equal tight
xlim([-10 10]), ylim([-10 10])

% Compute centroids, etc
G = computeGeometry(G);


%% Set up simulation
% Gravity
gravity reset off

% Rock
rock.poro = 0.2*ones(G.cells.num, 1);
rock.perm = 80*milli*darcy*ones(G.cells.num, 1);
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.rocknum    = ones(G.cells.num, 1);

% Fluid
if vaporizedOil == 1
    fn = 'ls-proj/radialCO2storage/co2water40CwetGas.DATA';    
else
    fn = 'ls-proj/radialCO2storage/co2water40C.DATA'; 
end
deck  = convertDeckUnits(readEclipseDeck(fn));
cP    = deck.PROPS.ROCK(:, 2);
fluid = initDeckADIFluid(deck);

% Initialize
p0 = 10.5*10^6*ones(G.cells.num, 1);
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             
rv0 = zeros(G.cells.num, 1);                                 
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% Wells and timesteps (well perforation thickness = 50m, i.e. all depth)
dist = pdist2(G.cells.centroids, [0 0 1000]);
[~, perforatedCell] = min(dist,[],1);
if G.numLayers > 1
    perforatedCell = [perforatedCell zeros(1, G.numLayers-1)];
    for n=2:G.numLayers
       perforatedCell(n) = perforatedCell(1)+G.layerSize*(n-1); 
    end
end

injTime = 50*year;
injRate = 10^6*meter^3/day;                                                 % [Sm^3/s]
rhoInj  = fluid.rhoGS;
mRate   = injRate*rhoInj*year/10^9;                                         % Mt/y
W = addWell([ ], G, rock, perforatedCell, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injRate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', min(G.faces.centroids(:,3)), 'radius', 0.1);    
        
simTime     = 50*year;
reportTimes = [[0.25, 1, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
               240, 300, 365, 456.25, 547.5, 638.75, 730, 821.25]*day, ...
               [2.5:0.5:22.5, 23:50]*year];
timesteps   = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

% Model
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);
if vaporizedOil == 1
    model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'vapoil', true, 'water', false);
else
    model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
end

% Acceleration
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
model = model.validateModel();

ncomp   = model.getNumberOfComponents();
solver  = AMGCL_CPRSolverAD('tolerance', 1e-6, 'block_size', ncomp, ...
                            'useSYMRCMOrdering', true, ...
                            'coarsening', 'aggregation', 'relaxation', 'ilu0');
nls     = NonLinearSolver('LinearSolver', solver);
%nls = [];

% Changes to model
model.operators.p0          = state0.pressure;                              % Needed for MyPvMult
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));

% BCs - Laterally infinite, so we use volume multipliers at the boundary
bc = [];
vertExtFaceId = abs(sqrt(G.faces.centroids(:,1).^2 + G.faces.centroids(:,2).^2) - 5e+4) < 50;
cellsExt = unique([G.faces.neighbors(vertExtFaceId, 1); G.faces.neighbors(vertExtFaceId, 2)]); 
cellsExt(cellsExt==0) = [];
% check: figure(1); hold on; plotCellData(G, rock.poro, cellsExt)
model.operators.pv(cellsExt) = model.operators.pv(cellsExt)*1e4;


%% Run simulation
% we run a parallel simulation with N threads
N = 4;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                   % Specify threads manually
schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);                      % Simple schedule, this injects for the total simTime
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'NonLinearSolver', nls, ...
                                                'Verbose', true);

%% Plots
model = model.validateModel();
bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    bhp(n) = wellSols{n}.bhp;
    %states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
end

cmap = flipud(cmap_agu());
figure(3)
colormap(cmap)
plotToolbar(G, states, 'FaceAlpha', 0.85);
view(53, 63), camproj perspective; axis equal tight
plotWell(G, W, 'color', 'k', 'height', 5);
c = colorbar; caxis([0 1])

figure(4)
colormap(cmap)
plotToolbar(G, states, G.layerSize+1:2*G.layerSize);
view(53, 63), camproj perspective; axis equal tight
plotWell(G, W, 'color', 'k', 'height', 5);
c = colorbar; caxis([0 1]); xlim([-50 50]); ylim([-50 50])

figure(5)
subplot(1,2,1)
plot([0 reportTimes/year], [105 bhp./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5])
xlabel('Time [y]')
ylabel('p [bar]')
grid on
title('BHP during injection period (50y)')

subplot(1,2,2)
plot([0 reportTimes(1:6)/day], [105 bhp(1:6)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5])
xlabel('Time [days]')
xticks([0 1 7 14 21 30])
hold on
text(0.75,188,'6h')
ylabel('p [bar]')
grid on
title('Detail during first month')

if saveData == 1
    disp('Data will be saved when any key is pressed. To cancel, press ctrl+c')
    pause
    if vaporizedOil == 1
        save('data_vap') 
    else
        save('data_noVap')
    end
end

