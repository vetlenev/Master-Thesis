%% Perform CO2 injection (blackoil module)
%--------------------------------------------------------------------------
% Mesh from simpleExtrudedMesh.
% 3 saturation regions, kr hysteresis, fault Pc
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
mrstVerbose on


%% Options
hyster = 'off';
mesh = 750;                                     % 150 (coarse), 250 (medium), 500 (fine)
wellno = 1;                                     % n injectors
rate = 1;                                       % mL/min should be ok up to a few mL/min


%% Save data
folderName = ['test_mesh', num2str(mesh), '_hyster', hyster, ...
              '_rate_mm_', num2str(rate)];
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
pth = fullfile(mrstPath('ls-proj'), 'fluidflower/test_upr_mesh/');
if mesh == 150
    G_dat  = load(fullfile(pth, 'simpleExtrudedG_nc150.mat'));
    sealID = [1 3 5 8];
    uGr = {5, [4 7], [1 8], [2 9], 3, 6};
elseif mesh == 250
    G_dat  = load(fullfile(pth, 'simpleExtrudedG_nc250.mat'));
    sealID = [2 3 5 8];  
    uGr = {2, [1, 7], [3 8], [4 9], 5, 6};
elseif mesh == 500
    G_dat  = load(fullfile(pth, 'simpleExtrudedG_nc500.mat'));
    sealID = [1 4 8 5];  
    uGr = {4, [3, 7], [1 8], [2 9], 5, 6};
elseif mesh == 750
    G_dat  = load(fullfile(pth, 'simpleExtrudedG_composite_nc750.mat'));
    sealID = [1 4 8 5];  
    uGr = {4, [3, 7], [1 8], [2 9], 5, 6};
end
G = G_dat.G;
faultID = 6;

% Plot
% fig3D();
% %plotCellData(G, G_dat.p)
% Ncompart = max(G_dat.p);
% cmap = copper(Ncompart);
% colr = [1, 7, 4, 8, 2, 5];
% for n=1:numel(uGr)
%     plotCellData(G, colr(n)*ones(sum(ismember(G_dat.p, uGr{n})), 1), ...
%                  ismember(G_dat.p, uGr{n}), 'edgealpha', 0.2)
% end
% outlineCoarseGrid(G, G_dat.compartID,'EdgeColor','w','LineWidth',2);
% setAxProps(gca), %camlight();
% colormap(copper); %c = colorbar; set(c, 'YTick', sort(colr));
% axis equal off
% ylim([0 1]), zlim([0 1])


%% Rock
% Layers follow uGr (bot to top, left to right for faulted, fault at the end)
rock.poro = nan(G.cells.num, 1);
rock.perm = nan(G.cells.num, 1);    % isotropic
poroMean = [0.25, 0.4, 0.3, 0.45, 0.27, 0.35];
permMean = [100, 6000, 1000, 10000, 500, 3000];   % [mD]
permtype = 'logn';
for n=1:numel(uGr)
    % Cell ids of corresponding unit
    cid = ismember(G_dat.compartID, uGr{n});
    nc = sum(cid);
    
    % Porosity
    r = 0.02; % half range
    rock.poro(cid) = (poroMean(n) - r) + rand(nc, 1)*(2*r);
    
    % Permeability
    if strcmp(permtype, 'randm')    % random
        logr = 0.2;
        rock.perm(cid) = 10.^(rand(nc, 1)*(2*logr) + ...
                              (log10(permMean(n)) - logr));
    elseif strcmp(permtype, 'logn') % lognormal
        N = [6, 150, 100];
        if n==faultID
           N = [150, 100, 1]; 
        end
        K = reshape(logNormLayers(N, 1, 'sigma', 0.5, 'a', 0.2), N);
        q = permMean(n)/7.3;
        rock.perm(cid) = sampleFromBox(G, q.*K, find(cid));
    end
end
rock.perm = rock.perm*(milli*darcy); % SI units

% Plot permeability
%fig3D(); plotCellData(G, log10(rock.perm/(milli*darcy)), 'edgealpha', 0.2)
%setAxProps(gca), colormap(copper), c = colorbar;
%axis equal off


%% Fluids
% 3 Saturation regions
fault.satRegNum = 3;
rock.regions.saturation = ones(G.cells.num, 1);                              % reservoir units
rock.regions.saturation(ismember(G_dat.compartID, sealID)) = 2;
rock.regions.saturation(ismember(G_dat.compartID, faultID)) = fault.satRegNum;

% 1 Regions for compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% load deck and initialize fluid
fn  = fullfile(pth, 'krHyst_fPcFault_co2brine_shallow.DATA');
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);

% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Scaled Pc in fault
refPc.val  = fluid.pcOG{fault.satRegNum};
refPc.poro = 0.07;                            % Table 3.3, sample 5 in Treviño & Meckel, 2017 (Pc curve is the green one in Fig. 3.4a)
refPc.perm = 0.00208*milli*darcy;             % "
fault.poro = rock.poro(G_dat.compartID==faultID);
fault.k.vals = rock.perm(G_dat.compartID==faultID);
fluid = assignFaultPc(fluid, fault, refPc);

% Hysteresis
if isfield(fluid, 'krHyst') && fluid.krHyst == 1
   numReg = max(rock.regions.saturation);
   fluid.krHyst = 1 + numReg;                % porous rock is first imb region (4), fault is third (6)
   rock.regions.imbibition = rock.regions.saturation + numReg;
   fluid = addScanKr(fluid, rock.regions.imbibition);
end


%% Initialize
gravity reset on
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
water_column = 0.3; % m
p_r = 1*barsa + g*rho_wr*water_column; % p at shallowest z
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% plot
% fig3D(); plotCellData(G, p0/barsa, 'edgealpha', 0.2)
% setAxProps(gca), colormap(jet), c = colorbar;
% axis equal off


%% Wells
t = [5*minute [0.5, 14]*day];     % rampup, injtime, simtime                          
injrate = rate*(milli*litre)/(minute*wellno);
injvol = injrate*t(2);
rhoInj = fluid.rhoGS;
mrate = injvol*rhoInj/(t(2)*wellno);                                      % [Sm^3/s]
wellInx = find(G_dat.wellNo==1);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', G.cells.centroids(wellInx, G.griddim), ...
            'Radius', 0.9*1e-3);             % refDepth must be included for 'bhp' wells.
reportTimes = [(1:t(1)/minute)*minute, ...                      % rampup inj rate
               [10:10:120 150:30:720]*minute, ...               % injection
               (722:1:740)*minute, ...                          % rampdown inj rate
               [12.5:.5:48 49:72 76:2:96 100:4:120]*hour, ...   % rest
               (5.5:.25:t(end)/day)*day];                       % rest
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==t(end), 'sum of timesteps must equal simTime')


%% Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'water', false);
                         
% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true);
model = model.validateModel();

nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, ...
                                      'targetChangeAbs', 0.25);
nls.timeStepSelector = stepSel;
nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');
% nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-3, 'Solver', 'gmres', ...
%                                           'preconditioner', 'amg', ...
%                                           'coarsening', 'smoothed_aggregation');
%nls.useRelaxation = 1;
%nls.relaxationType = 'sor';
nls.enforceResidualDecrease = 1;
nls.useLinesearch = 1;
nls.maxIterations = 12;
nls.maxTimestepCuts = 8;
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

%% BCs
bctype = 'pvm';
% Find top faces
f = false(G.faces.num, 1);
f(boundaryFaces(G)) = true;
zmean = mean(G.cells.centroids(:, 3));
f = all([G.faces.centroids(:, 3) < zmean, f], 2);
f = all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 2))], 2);
f = find(all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 1))], 2));

%fig3D(); plotFaces(G, f, 'edgecolor', 'r')
%setAxProps(gca), colormap(jet), c = colorbar;
%axis equal off

if strcmp(bctype, 'cp')
    [z_0, z_max] = deal(min(G.faces.centroids(f,3)), max(G.faces.centroids(f,3)));
    equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
    fp_val = reshape(deval(equil, G.faces.centroids(f,3)), [], 1);  clear equil
    bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0]);
elseif strcmp(bctype, 'pvm')
    cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
    cellsext(cellsext==0) = [];
    %fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
    %setAxProps(gca), axis equal off, view([90 0])
    model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e4; 
    bc = [];
end

%% Schedule
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
v = injrate;
injrates = [0.05*injrate 0.1*injrate 0.2*injrate 0.5*injrate injrate ...
            0.9*v 0.7*v 0.5*v 0.3*v 0.1*v 0.05*v 0.01*v ...
            0.005*v 0.001*v 0];
tmp = cell(numel(injrates), 1);                                             % create 2 schedules
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
for n=1:numel(injrates)
    schedule.control(n).W = W;                                              % field 1 used during injection
    schedule.control(n).W.val = injrates(n);                                % nr of wells must be the same for the entire simulation
    schedule.control(n).bc = bc;
end
% Smaller than rampup
idStep = find(cumsum(schedule.step.val) < t(1));
schedule.step.control(idStep) = 1:max(idStep);
schedule.step.control(idStep(end)+1:end) = max(idStep)+1;

% Larger than injtime
idStep2 = find(cumsum(schedule.step.val) > t(2), 1);
schedule.step.control(idStep2:idStep2+numel(injrates)-2-max(idStep)) = ...
    max(schedule.step.control)+1:numel(injrates); % set timesteps after injTime to use Well field 2
schedule.step.control(find(schedule.step.control > max(idStep)+1, 1, 'last')+1:end) = ...
                                                           numel(injrates); 

%% Simulation
N = 8;
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
bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    bhp(n) = wellSols{n}.bhp;
end
states{1}.reg = model.rock.regions.saturation;


% basic overview
rescells = find(ismember(G_dat.compartID, uGr{2}));
fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); 
plotGrid(G, rescells, 'edgecolor', [0.5 0.5 0.5]);
setAxProps(gca), colormap(jet), c = colorbar; %clim([0 40000])
axis equal off
view([90 0])