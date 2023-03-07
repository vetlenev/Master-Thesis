%% Black Oil Test: relative permeability hysteresis
%
% 3D mesh with 2 fluids (water and CO2)
%
% Goals:
% 1. Run a simulation with 1 region, no dissolution, with and without
%    hysteresis. Assuming that the formation is strongly water-wet, kr
%    hysteresis is only considered for the CO2 phase. Accounting for
%    hysteresis is expected to greatly reduce the volume of free CO2 at the
%    top of the injection layer, since residual trappping occurs as the
%    plume migrates upwards due to secondary imbibition. PUNQ-S3 model, see
%    Juanes et al., 2006.
%
% 2. Run a simulation with 2 regions and dissolution, with and without
%    hysteresis, to check that it works with other simulation complexities.
%    Use example bo_test_CO2brine_regs.m as guide.
%
%
clear, clc, close all
mrstModule add ls-tests ls-utils ad-props deckformat mrst-gui ad-core ad-blackoil linearsolvers
gravity reset on;

% Input deck and generate grid, rock and fluid.
fn   = 'mrst-dev/mrst-ls/ls-tests/hyst/punq-s3/CASE2_regi.DATA';
deck = convertDeckUnits(readEclipseDeck(fn));

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

id = G.cells.centroids(:,2) > 3500;
% 2 saturation regions
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.saturation(id) = 2;
deck.regions.SATNUM = rock.regions.saturation;
% 2 imbibition regions
rock.regions.imbibition(id) = 4;
deck.REGIONS.IMBNUM(id) = 4;

fluid = initDeckADIFluid(deck);
fluid.bO = [];                           % why is this needed for mex?
% Indicate imbibition regs with hysteresis if not all cells have hysteresis 
% active. Field should not exist if hysteresis is not used anywhere.
fluid.krHyst = 4;                        
fluid = addScanKr(fluid, rock.regions.imbibition);

% Modify to match Juanes et al. (2006) modifications.
%rock.poro = rock.poro + (0.2 - mean(rock.poro)); 
shiftZup = 840 - min(G.faces.centroids(:, 3));
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + shiftZup;
G.cells.centroids(:, 3) = G.cells.centroids(:, 3) + shiftZup;
G.faces.centroids(:, 3) = G.faces.centroids(:, 3) + shiftZup;

% Select model automatically
model = selectModelFromDeck(G, rock, fluid, deck);   % GenericBlackOilModel, initially error in relPerm evaluation

% Initialize reservoir
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
p_r = g*fluid.rhoWS*z_0;   
equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = 0;                                 % immiscible
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);
pavg = mean(p0([schedule.control(1).W.cells]));
pvol = sum(model.operators.pv);
nwell = 8;
tinj  = 10*year;
fvol  = fluid.rhoGS*fluid.bG(pavg)/fluid.rhoGS; % m3 reservoir (0.15pv) to m3 surface
[schedule.control(1).W.val] = deal(0.15*pvol*fvol/(nwell*tinj));
[schedule.control(2).W.val] = deal(0);          % control 2 no injection
[schedule.control(1).W.type] = deal('rate');    % supported types are 'rate' or 'bhp'
[schedule.control(2).W.type] = deal('rate');

% Timesteps
repTimes = [0.05:0.05:6 6.1:.1:10 10.5:0.5:20 21:50 52:2:100 ...
            105:5:200 210:10:500]'*year;
timesteps = [repTimes(1); diff(repTimes)];
schedule.step.val = timesteps;
schedule.step.control = ones(numel(timesteps), 1)+1;
schedule.step.control(1:20) = 1;

% Preliminary plots
figure;
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
             'FaceAlpha', 0.95, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
set(gca,'YDir','reverse')
title('Permeability (mD)')
axis tight, view(-40, 70), c = colorbar('SouthOutside'); caxis([0.5 1000])

% BC
bc = [];
%idext1 = G.faces.neighbors(:,2) == 0; idext2 = G.faces.neighbors(:,1) == 0;
%cellsext = [G.faces.neighbors(idext1, 1); G.faces.neighbors(idext2, 2)];
idb = G.faces.normals(:, 3) == 0;
cellsn0 = G.faces.neighbors(idb, :);
idb1 = cellsn0(:,2) == 0; idb2 = cellsn0(:,1) == 0;
cellsb = unique([cellsn0(idb1, 1); cellsn0(idb2, 2)]);
cellsb(cellsb==72) = []; cellsb(cellsb==86) = [];

% Changes to model - mex acceleration t = 83s (CASE2_regi.DATA, MRST 2019b).
model = model.validateModel();
ncomp  = model.getNumberOfComponents();
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
%model = model.validateModel();
solver = AMGCL_CPRSolverAD('tolerance', 1e-6, 'block_size', ncomp, ...
    'useSYMRCMOrdering', true, 'coarsening', 'aggregation', 'relaxation', 'ilu0');
nls = NonLinearSolver('LinearSolver', solver);
% nls = NonLinearSolver('LinearSolver', solver, 'maxIterations', 10, ...
%                       'maxTimestepCuts', 8);
%nls = [];

model.operators.pv(cellsb) = model.operators.pv(cellsb)*1e3;
%model.FlowPropertyFunctions = MyFlowPropertyFunctions(model); old
model = model.validateModel(); % Set up the state function groups
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
%model.OutputStateFunctions = {'ComponentTotalMass','RelativePermeability', 'Density', ...
%                              'Mobility', 'PhasePressures', 'PoreVolume', 'ShrinkageFactors', 'Viscosity'};

% Simulation
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
figure(13)
subplot(1,2,1)
plot(s, kr, '.-');
title('kr_g')

subplot(1,2,2)
stairs(s, kr_d);
title('d kr_g / d sg')