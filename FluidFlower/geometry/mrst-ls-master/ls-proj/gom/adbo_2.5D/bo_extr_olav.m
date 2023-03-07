%% Base case model (2D extruded) for GoM CO2 fault migration

%____________________________________________________________________
clc, clear, close all force;
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils

%% CHOOSE SIMULATION OPTIONS
fname     = 'sc2_hyst_62x23523';               % set fname = [] to not save 'case2_38x45'
figs      = [0, 0, 0, 0];                
scenario  = 2;
studyCase = 'sc2';
[mesh, opt, wells, t, resPlots] = simOpts(scenario, studyCase);


%% GET SIMULATION VARIABLES
% Domain grid & set gravity
[G, cellw] = getMesh_bo(mesh, figs(1));

% Simulation grid, Rock, unit cell ids, fault props
[G, rock, ucids, fault] = getRockParams_bo(mesh, G, opt, figs(2));

% Fluid & Regions
[fluid, rock, cP, fault, refPc] = getFluid_bo(mesh, G, rock, opt, ...
                                              ucids.unit_cell_ids, fault, figs(3));

% Initialize
state0 = init_bo(G, fluid, opt.watCol, ucids, figs(4));

% Wells & timesteps
[W, timesteps, wellInx] = getWell_bo(G, rock, fluid, wells, t);

%%
% Load module setup from Lluis here!
%% MODEL SETUP and BC
% modify pvMult to take into account:
%   - In rock compressibility: initial p at each cell
% modify fluid.Pc to take into account:
%   - In capillary pressure: pc is a fcn of poro and perm in the FZ (selected cases only)
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);
if strcmp(opt.faultPc, 'scaled')
    [fluid] = assignFaultPc(fluid, fault, refPc);
end

% Define base model class
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'water', false);

% Mex and Solver
[model, nls] = setAcceleration_bo(model, opt);
% model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'deferredAssembly', true, 'rowMajor', true);
% model = model.validateModel();
% nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
% nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-3);

% Changes to model. Note that this must be done after setting the autodiff
% backend (in setAcceleration_bo) since a copy of the backend will be made
% when setting the flow property functions.
model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Density', BlackOilDensity(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Mobility', Mobility(model));
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PressureReductionFactors', BlackOilPressureReductionFactors(model));
% model.FlowPropertyFunctions = model.FlowP ropertyFunctions.setStateFunction('ShrinkageFactors', BlackOilShrinkageFactors(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Viscosity', BlackOilViscosity(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RsMax', RsMax(model));
%model.FlowPropertyFunctions = MyFlowPropertyFunctions(model);
model.minimumPressure = 1*barsa;

% Boundary Conditions
[bc, model] = getBC_bo(G, fluid, model, opt, mesh, ucids);
schedule = getSchedule_bo(timesteps, W, bc, t);

%%
% for i = 1:numel(model.fluid.pvMultR)
%     model.fluid.pvMultR{i} = @(p) model.fluid.pvMultR{i}(p, model.operators.p0);
% end
%% SIMULATION
% Run a parallel simulation with N threads
% N = 12;
% maxNumCompThreads(N);
% nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
% [wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                                 'NonLinearSolver', nls, ...
%                                                 'Verbose', true);
%                                             
%%
nls.LinearSolver.tolerance = 1e-3;
nls.LinearSolver.setSolver('bicgstab')

problem = packSimulationProblem(state0, model, schedule, 'bo_extr', 'NonLinearSolver', nls);
simulatePackedProblem(problem);
%%
[ws, states, reports] = getPackedSimulatorOutput(problem, 'readStatesFromDisk', false);
plotWellSols(ws, cumsum(schedule.step.val))
%%
timing = getReportTimings(reports);
figure;
plot([timing.Iterations])
xlabel('Step')
ylabel('Iterations / step')
%%
figure;
plot(cumsum([timing.Total])/hour)
xlabel('Step')
ylabel('Elapsed runtime')


%% RESULTS, PLOTS AND SAVE DATA
model = model.validateModel();
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.pc = pc{1,2};
    states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
end
states{1}.reg = rock.regions.rocknum;

% Save data
if ~isempty(fname)
    disp(['ATTENTION: Data will be saved in: ' pwd])
    %disp('Press Ctrl+C to CANCEL or any other key to SAVE data.')
    
    save([fname '.mat'], '-v7.3')
    %save([fname '_states.mat'], 'states', '-v7.3') % larger than 2GB
end

% Plots
resPlots_bo(G, rock, fluid, model, ucids, W, states, state0, ncont, mesh, resPlots) 


