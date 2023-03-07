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
saveData     = 1;
vaporizedOil = 0;               % gas dissolution in brine always
importTriangleGrid = 1;


%% Create radial grid
if importTriangleGrid ~= 1
    % (see fig. 3.25 in pg 86 and sect 3.5.2 in MRST book).
    % Scale grid to 50km radius and 50m thickness, and move top depth to 1000m.
    P = [];
    for r = exp([-9:.1:-5 -4.95:0.05:-3 -2.975:0.025:0])
        [x,y,z] = cylinder(r,110);
        P = [P [x(1,:); y(1,:)]];
    end
    P = unique([P'; 0 0],'rows');
    P = P.*50*kilo*meter;
    G = makeLayeredGrid(pebi(triangleGrid(P)), 1);
    k  = G.nodes.coords(:, 3) > 0;
    G.nodes.coords(k, 3) = 50;
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000;
else
    meshpath = 'mrst-dev/mrst-ls/ls-proj/radialCO2storage/mesh/nodes_coordinates.dat';
    %meshpath = 'mrst-dev/mrst-ls/ls-proj/radialCO2storage/mesh/triangle2/nodes_coordinates.dat';
    vertices = dlmread(meshpath);
    x = vertices(:,1); y = vertices(:, 2).*-1;
    G = makeLayeredGrid(triangleGrid([x, y]), [17.5 15 17.5]);
    %k  = G.nodes.coords(:, 3) > 0;
    %G.nodes.coords(k, 3) = 50;
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000;
end
% Visualize grid
% figure(1)
% plotGrid(G, 'FaceColor',[.8 .8 .8]); view(25,40), axis tight
% 
% figure(2)
% subplot(1,2,1)
% plotGrid(G,'FaceColor',[.8 .8 .8]); view(0,90), axis equal tight
% xlim([-1000 1000]), ylim([-1000 1000])
% 
% subplot(1,2,2)
% plotGrid(G,'FaceColor',[.8 .8 .8]); view(0,90), axis equal tight
% xlim([-50 50]), ylim([-50 50])

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
    %fn = 'ls-proj/radialCO2storage/co2water50C.DATA';
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
reportTimes = [[0.5 1 6]*hour [1, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
               240, 300, 365, 456.25, 547.5, 638.75, 730]*day, ...
               [2.2:.2:5 5.5:0.5:20, 21:50]*year];
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
% nls = [];

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


%% Run simulation
% we run a parallel simulation with N threads
%pause
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
    states{n}.FlowProps.Viscosity = model.getProp(states{n}, 'Viscosity');
    states{n}.FlowProps.phaseDensity = model.getProp(states{n}, 'Density');
end

% General
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
c = colorbar; caxis([0 1]); xlim([-1300 1300]); ylim([-1300 1300])

% BHP
figure(5)
subplot(2,1,1)
plot([0 reportTimes(4:end)/year], [105 bhp(4:end)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
%plot([0 reportTimes/year], [105 bhp./10^5], '.-k')
xlabel('time [y]')
ylabel('p [bar]')
ylim([100 171])
yticks([100 110 120 130 140 150 160 170])
grid on
title('BHP during injection period (50y)')

subplot(2,1,2)
plot([0 reportTimes(4:16)/day], [105 bhp(4:16)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
xlabel('time [days]')
xticks([0 7 30 60 90 120 150 180 210 240 270 300 330 365])
yticks([100 110 120 130 140 150 160 165])
ylim([100 165])
hold on
%text(0.75,188,'6h')
ylabel('p [bar]')
grid on
xlim([0 365])
title('Detail during first year')

% Dp vs time at different distances 
% (center depth and 10m, 100m, 1000m, 5000m, 25km, 50km)
centralLayerIds = transpose(1:G.cells.num);
centralLayerIds = centralLayerIds(G.cells.centroids(:,3) == 1025);

cellsClosestTo = [10, 0, 1025; 100, 0, 1025; 1000, 0, 1025; 5000, 0, 1025;...
                  25000, 0, 1025; 50000, 0, 1025];
dist = pdist2(G.cells.centroids, cellsClosestTo);
[~, monitoringCells] = min(dist,[],1);      

dpCells = zeros(numel(states), numel(monitoringCells));
for n=1:numel(states)
    dpCells(n, :) = states{n}.dp(monitoringCells);
end

figure(6)
subplot(2,1,1)
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,1)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'markerSize', 3)
hold on
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,2)./10^5], '-om', 'markerFaceColor', [255, 125, 255]/255, 'markerSize', 3)
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,3)./10^5], '-or', 'markerFaceColor', [255, 125, 125]/255, 'markerSize', 3)
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,4)./10^5], '-og', 'markerFaceColor', [125, 255, 125]/255, 'markerSize', 3)
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,5)./10^5], '-ob', 'markerFaceColor', [125, 125, 255]/255, 'markerSize', 3)
plot([0 reportTimes(4:end)/year], [0; dpCells(4:end,6)./10^5], '-oc', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 3)
hold off
xlabel('time [y]')
ylabel('\Deltap [bar]')
grid on
%legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
ylim([0 50])
title('Pressure change vs time at different distances from injection well')

% subplot(3,1,2)
% plot([0 reportTimes/year], [0; dpCells(:,1)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'markerSize', 2)
% hold on
% plot([0 reportTimes/year], [0; dpCells(:,2)./10^5], '-om', 'markerFaceColor', [255, 125, 255]/255, 'markerSize', 2)
% plot([0 reportTimes/year], [0; dpCells(:,3)./10^5], '-or', 'markerFaceColor', [255, 125, 125]/255, 'markerSize', 2)
% plot([0 reportTimes/year], [0; dpCells(:,4)./10^5], '-og', 'markerFaceColor', [125, 255, 125]/255, 'markerSize', 2)
% plot([0 reportTimes/year], [0; dpCells(:,5)./10^5], '-ob', 'markerFaceColor', [125, 125, 255]/255, 'markerSize', 2)
% plot([0 reportTimes/year], [0; dpCells(:,6)./10^5], '-oc', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 2)
% hold off
% xlabel('time [y]')
% ylabel('\Deltap [bar]')
% grid on
% %legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
% xlim([4 20])
% title('Zoom to oscillations')

subplot(2,1,2)
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,1)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'markerSize', 3)
hold on
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,2)./10^5], '-om', 'markerFaceColor', [255, 125, 255]/255, 'markerSize', 3)
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,3)./10^5], '-or', 'markerFaceColor', [255, 125, 125]/255, 'markerSize', 3)
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,4)./10^5], '-og', 'markerFaceColor', [125, 255, 125]/255, 'markerSize', 3)
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,5)./10^5], '-ob', 'markerFaceColor', [125, 125, 255]/255, 'markerSize', 3)
plot([0 reportTimes(4:16)/day], [0; dpCells(4:16,6)./10^5], '-oc', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 3)
hold off
legend({'9.6m', '95m', '984m', '4968m', '24.841km', '49.683km'}, 'box', 'on', 'fontsize', 8)
xlabel('time [day]')
ylabel('\Deltap [bar]')
grid on
xticks([0 7 30 60 90 120 150 180 210 240 270 300 330 365])
ylim([0 50])
xlim([0 365])
title('Detailed first year view')

% p vs distance at t = 15y
tstep = find(reportTimes/year == 15);
ycoords = [0.4:0.5:2000 2001:1:50000];   
cellsClosestTo = [zeros(numel(ycoords), 1), ycoords', 1025*ones(numel(ycoords), 1)];
dist = pdist2(G.cells.centroids, cellsClosestTo);
[~, monitoringCells2] = min(dist,[],1);      
monitoringCells2 = unique(monitoringCells2);
[~, idm] = sort(G.cells.centroids(monitoringCells2, 2));
monitoringCells2 = monitoringCells2(idm);
rmids = false(numel(monitoringCells2), 1);
for n=2:numel(monitoringCells2)
    distx = abs(G.cells.centroids(monitoringCells2(n), 1) - G.cells.centroids(monitoringCells2(n-1), 1));
    disty = abs(G.cells.centroids(monitoringCells2(n), 2) - G.cells.centroids(monitoringCells2(n-1), 2));
    if distx > disty
        rmids(n) = true;
    end
end
monitoringCells2(rmids) = [];
pCells = states{tstep}.pressure(monitoringCells2);
sCells = states{tstep}.s(monitoringCells2, 2);
dist = pdist2(G.cells.centroids(monitoringCells2, 1:2), G.cells.centroids(perforatedCell(1), 1:2));

figure(7)
subplot(2,1,1)
plot(dist, pCells/10^5, '-b', 'linewidth', 1)
grid on
xlabel('distance [m]')
ylabel('p [bar]')
title('Reservoir p vs distance from wellbore at t=15y')
yticks([100 110 120 130 140 150 160])

subplot(2,1,2)
semilogx(dist, pCells/10^5, '-b', 'linewidth', 1)
grid on
xlabel('distance [m]')
ylabel('p [bar]')
title('Reservoir p vs distance from wellbore at t=15y')
xlim([0.1 10^5])

% S vs distance 
figure(8)
subplot(2,1,1)
plot(dist, sCells, '-', 'linewidth', 1.25, 'color', [0.3 0.8 0.3])
grid on
xlabel('distance [m]')
ylabel('p [bar]')
title('CO_2 saturation vs distance from wellbore at t=15y')
xlim([0.1 10^5])
yticks(0:.05:.5)
ylim([0 0.5])

subplot(2,1,2)
semilogx(dist, sCells, '-g', 'linewidth', 1.25, 'color', [0.3 0.8 0.3])
grid on
xlabel('distance [m]')
ylabel('p [bar]')
title('CO_2 saturation vs distance from wellbore at t=15y')
xlim([0.1 10^5])
yticks(0:.05:.5)
ylim([0 0.5])

% P and viscosity
%rhoG = fluid.bG(pCells).*fluid.rhoGS;
%muG  = fluid.muG(pCells);
%rs   = fluid.rsSat(pCells);
%  rhoB = fluid.bO(p, rs, true(numel(p), 1)).*(rs.*fluid.rhoGS + fluid.rhoOS);
%  muB  = fluid.muO(p, rs, true(numel(p), 1));
rhos = [states{tstep}.FlowProps.phaseDensity{1}(monitoringCells2), ...
        states{tstep}.FlowProps.phaseDensity{2}(monitoringCells2)];
mus = [states{tstep}.FlowProps.Viscosity{1}(monitoringCells2), ...
       states{tstep}.FlowProps.Viscosity{2}(monitoringCells2)];
   
h = figure(9);
subplot(1,2,1)
semilogx(dist, rhos(:,1), '-', 'color', clrs(4,:), 'linewidth', 1.5);
hold on; 
semilogx(dist, rhos(:,2), '-', 'color', clrs(3,:), 'linewidth', 1.5);
hold off,  grid on; ax = gca; xlim([0.1 50000])
xlabel('distance [m]','fontsize', 14, latx{:})
ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
set(ax,'Xtick',[0.1 1 10 100 1000 10^4]);
set(ax,'Ytick',500:100:1100); ylim([500 1100]);
legend({'Brine', 'CO$_2$'}, latx{:}, 'fontSize', 10)
title('Density at t=15y')
   
subplot(1,2,2)
loglog(dist, mus(:,1)/(centi*poise), 'color', clrs(4,:), 'LineStyle', '-', 'linewidth', 1.5);
hold on; 
loglog(dist, mus(:,2)/(centi*poise), '-', 'color', clrs(3,:), 'linewidth', 1.5);
hold off;  %xlim([min(p) max(p)]./barsa)
ylabel('$\mu$ [cP]','fontsize', 14, latx{:}); ax = gca;
set(ax,'Ytick',[0 0.15 0.3 0.45 0.6]);
ylim([0.01 0.9])
set(ax,'Xtick',[0.1 1 10 100 1000 10^4]);
set(ax,'Ytick', [0.01:0.01:0.1 0.2:.1:.9]);
set(h, 'Position', [600, 600, 380, 320])
grid on
title('Viscosity at t=15y')


if saveData == 1
    disp('Data will be saved when any key is pressed. To cancel, press ctrl+c')
    %pause
    close all force
    if vaporizedOil == 1
        save('data_vap') 
    else
        save('data_noVap_triangleMesh.mat', '-v7.3')
    end
end

