function [model, nls] = getModel(G, rock, fluid, state0, opt)
%
%
%
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'vapoil', opt.vapoil, 'water', false);
                         
% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true);
model.toleranceCNV = 5e-3;      
%model.toleranceMB = 3e-7;
model.dsMaxAbs = 0.05;
%model.useCNVConvergence = 0;
%model.nonlinearTolerance = 1e-6;
model = model.validateModel();

% Diffusion CO2 component
diffFlux = BlackOilComponentTotalFluxWithDiffusion(model);
diffFlux.componentDiffusion = [0 opt.D];  % opt.D = diffusion coeff.
diffFlux.faceAverage = true;
model.FlowDiscretization.ComponentTotalFlux = diffFlux;

nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
                          'useCPR', true);
%nls.LinearSolver.tolerance = 1e-6;
% nls = getNonLinearSolver(model, 'TimestepStrategy', 'none', ...
%                         'useCPR', true);
% stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, ...
%                                     'targetChangeAbs', 0.25);
% nls.timeStepSelector = stepSel;
if strcmp(opt.bctype, 'pvm')
    nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');
end
nls.LinearSolver.maxIterations = 25;
%nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'gmres', ...
%                                          'preconditioner', 'amg', ...
%                                          'coarsening', 'smoothed_aggregation');
%nls.useRelaxation = 1;
%nls.relaxationType = 'sor';
%nls.enforceResidualDecrease = 1;
% nls.useLinesearch = false;
nls.useLinesearch = true;
nls.useRelaxation = false;
nls.maxIterations = 12;      % 8, 6, 4, 2
%nls.minIterations = 0;     % Do NOT put minIterations to 0!
nls.maxTimestepCuts = 14;   % 10, 12, 14, 16
nls.acceptanceFactor = 5;  % 5, 10
%nls.alwaysUseStabilization=1;
%nls.convergenceIssues=true;
%nls = [];

% Model changes (after acceleration!)
model.operators.p0 = state0.pressure;                            % Needed for MyPvMult
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
if strcmp(opt.hyster, 'on')
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
    %model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));
end
model.minimumPressure = min(state0.pressure);

end
