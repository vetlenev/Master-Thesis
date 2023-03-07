%% Black Oil Test: relative permeability hysteresis
%
% 3D mesh with 2 fluids (water and CO2)
%
% Objective:
% Run a simulation with 1 region, no dissolution, with and without
% hysteresis. Assuming that the formation is strongly water-wet, kr
% hysteresis is only considered for the CO2 phase. Accounting for
% hysteresis is expected to greatly reduce the volume of free CO2 at the
% top of the injection layer, since residual trappping occurs as the
% plume migrates upwards due to secondary imbibition. PUNQ-S3 model, see
% Juanes et al., 2006.
%
% Note: The fn path (line 24) needs to be updated.
%
clear, clc, close all
mrstModule add ls-tests ls-utils ad-props deckformat mrst-gui ...
               ad-core ad-blackoil linearsolvers
gravity reset on;

% ------------- Input deck and generate grid, rock and fluid -------------
% We use the same input deck (CASE2.DATA) and properties that was used in 
% Juanes et al., WRR (2006), since MRST is able to import ECLIPSE decks.
folderName = 'Irate_0.15poreVol';
topDir = 'C:\Users\Lluis\matlab\sim_data\mrst\gcs3D\kr_hyst\';
fn   = [pwd '\mrst-dev\mrst-ls\ls-proj\gcs3D\kr-hysteresis\input_files\punq-s3\CASE2.DATA'];
deck = convertDeckUnits(readEclipseDeck(fn));
G = initEclipseGrid(deck);
shiftZup = 840 - min(G.nodes.coords(:, 3));
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + shiftZup;
G = computeGeometry(G);
%plotGrid(G); set(gca,'YDir','reverse');  view([-40 66])

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
%rock.poro = rock.poro + (0.2 - mean(rock.poro)); % mean poro = 0.2

fluid = initDeckADIFluid(deck);
fluid.bO = [];                           % why is this needed for mex?
rock.regions.saturation = ones(G.cells.num, 1); 
fluid.krHyst = 2;   % imbibition regions where hysteresis is active  
minSat = 0.02;      % saturation threshold for hysteresis activation
fluid = addScanKr(fluid, rock.regions.imbibition, minSat);

% --------------------- Plot hysteretic curve ----------------------------
% sgd = linspace(fluid.krPts.g(1,2),fluid.krPts.g(1,3),50)';
% krd = fluid.krG{1}(sgd);
% sgi = linspace(fluid.krPts.g(2,2),fluid.krPts.g(2,3),50)'; 
% kri = fluid.krG{2}(sgi);
% sgs = linspace(0.28,0.4,50)';
% krs = fluid.krGi{1}(sgs,repelem(0.4,50,1));
% figure(2)
% plot(sgd,krd,'-r','linewidth', 1.5, ...
%      'DisplayName', '$k_\mathrm{r,g}^\mathrm{d}$'), hold on
% plot(sgi,kri,'--r','linewidth', 1.5, ...
%      'DisplayName', '$k_\mathrm{r,g}^\mathrm{i}$')
% plot(sgs,krs,'-.r','linewidth', 1, ...
%      'DisplayName', '$k_\mathrm{r,g}^\mathrm{s}$'), hold off
% grid on
% xlabel('$S_\mathrm{g}$ [-]', 'interpreter', 'latex')
% ylabel('$k_\mathrm{r,g}$ [-]', 'interpreter', 'latex')
% title('Hysteretic gas relative permeability (Berea sst., gas-water, from Oak, 1990)')
% legend('interpreter', 'latex', 'fontsize', 12)
% xlim([0 1]); ylim([0 1])

% ----------------------- Model and ancceleration ------------------------
% Note that here, we multiply the pore volume of the boundary cells by
% 1000, following Juanes et al. (2006)
model = selectModelFromDeck(G, rock, fluid, deck); 
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
model = model.validateModel();
nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
                         'useCPR', true);
nls.LinearSolver = AMGCL_CPRSolverAD();
nls.useLinesearch = true;
nls.maxIterations = 15; 
nls.maxTimestepCuts = 12;
nls.acceptanceFactor = 2;

% idb = G.faces.normals(:, 3) == 0;
% cellsn0 = G.faces.neighbors(idb, :);
% idb1 = cellsn0(:,2) == 0; idb2 = cellsn0(:,1) == 0;
% cellsb = unique([cellsn0(idb1, 1); cellsn0(idb2, 2)]);
% cellsb(cellsb==72) = []; cellsb(cellsb==86) = [];
idAct = find(deck.GRID.ACTNUM);
idG = (1:prod(deck.GRID.cartDims))';
idG_mult = idG(deck.GRID.PORV==1000);
cellsb = ismember(idAct, idG_mult);
% figure(11); plotGrid(G, 'facecolor', 'none'); plotGrid(G,cellsb);
% set(gca,'YDir','reverse');  view([-40 66])
model.operators.pv(cellsb) = model.operators.pv(cellsb)*1e3;
model = model.validateModel(); % Set up the state function groups
model.FlowPropertyFunctions = ...
model.FlowPropertyFunctions.setStateFunction('RelativePermeability', ...
                                            MyRelativePermeability(model));
%model.OutputStateFunctions = {'ComponentTotalMass','RelativePermeability', 'Density', ...
%                              'Mobility', 'PhasePressures', 'PoreVolume', 'ShrinkageFactors', 'Viscosity'};


% -------------------------- Initialization -------------------------------
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
p_r = g*fluid.rhoWS*z_0;   
equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % fully saturated in water --> [1 0]
rs0 = 0;                                 % immiscible 
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% figure(11); plotCellData(G, p0/barsa); 
% set(gca,'YDir','reverse');  view([-40 66])
% colormap(jet), c = colorbar; c.Label.String = '$p_0$ [bar]'; 
% c.Label.Interpreter = 'latex'; c.Label.FontSize = 12;

% ------- Deck schedule into a MRST schedule by parsing the wells --------
schedule = convertDeckScheduleToMRST(model, deck);
pavg = mean(p0([schedule.control(1).W.cells]));
pvol = sum(rock.poro.*G.cells.volumes);
nwell = 8;
tinj  = 10*year;
fvol  = fluid.rhoGS*fluid.bG(pavg)/fluid.rhoGS; % m3 reservoir (0.15pv) to m3 surface
[schedule.control(1).W.val] = deal(0.15*pvol*fvol/(nwell*tinj));
[schedule.control(2).W.val] = deal(0);          % control 2 no injection
[schedule.control(1).W.type] = deal('rate');    % supported types are 'rate' or 'bhp'
[schedule.control(2).W.type] = deal('rate');
for n=1:nwell
    schedule.control(1).W(n).lims.bhp = 160*barsa;
end

% Timesteps
repTimes = [0.05:0.05:6 6.1:.1:10 10.5:0.5:20 21:50 52:2:100 ...
            105:5:200 210:10:500]'*year;
timesteps = [repTimes(1); diff(repTimes)];
schedule.step.val = timesteps;
% schedule.step.control = ones(numel(timesteps), 1)+1;
% schedule.step.control(1:20) = 1;
schedule.step.control = ones(numel(timesteps), 1);
schedule.step.control(161:end) = 2;

% Plot permeability and wells
% figure;
% plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
%              'FaceAlpha', 0.95, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
% plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
% set(gca,'YDir','reverse')
% title('Permeability (mD)')
% axis tight, view(-40, 70), c = colorbar('SouthOutside'); caxis([0.5 1000])


% ------------------------------ Simulation -------------------------------
N=1;
maxNumCompThreads(N)
nls.LinearSolver.amgcl_setup.nthreads = N; 
% [wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                                 'nonlinearsolver', nls, ...
%                                                 'Verbose', true);
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  

% ------------------------------- Results --------------------------------
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.pc = pc{1,2};
end

figure(10)
plotToolbar(G, states); view(-40, 70); set(gca,'YDir','reverse')

figure(11)
plotCellData(G, states{end}.s(:, 2), ...
             'FaceAlpha', 1, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
set(gca,'YDir','reverse')
title(['CO_2 saturation, t=' num2str(sum(schedule.step.val)/year) 'y']); 
colormap(jet)
axis tight, view(-40, 66), c = colorbar('SouthOutside'); caxis([0 1])

if strcmp(fn, 'mrst-2019b/mytests/punq-s3/CASE2_regi.DATA')
    figure(12)
    plotCellData(G, rock.regions.imbibition, ...
        'FaceAlpha', 0, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
    hold on
    allCells = 1:G.cells.num;
    plotCellData(G, rock.regions.imbibition, allCells(id), ...
        'FaceAlpha', 1, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
    set(gca,'YDir','reverse')
    title('Imb regions');
    colormap(jet)
    axis tight, view(-40, 66), c = colorbar('SouthOutside');
end

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
figure;
plot(s, kr, '.-');
title('kr_g')

figure;
stairs(s, kr_d);
title('d kr_g / d sg')