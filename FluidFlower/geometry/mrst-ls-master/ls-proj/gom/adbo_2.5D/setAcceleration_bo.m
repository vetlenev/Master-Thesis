function [model, nls] = setAcceleration_bo(model, opt)
%
%
%
ncomp  = model.getNumberOfComponents();
if opt.mex == 1
    % Set up mex-accelerated backend and reduced variable set for wells
    % We use the Diagonal autodiff backend to calculate derivatives. By
    % default, this uses only Matlab, so we also set the "useMex" flag to be
    % true. It will then use C++ versions of most discrete operators during
    % assembly.
    %model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
    %model.FacilityModel.primaryVariableSet = 'bhp';
    
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                    'deferredAssembly', true);

    % Run assembly tests (will compile backends if required)
    % Note: Benchmarks are obviously not accurate if compilation is performed.
    %testMexDiagonalOperators(model, 'block_size', 2);
end
%model.toleranceCNV = 2e-3;      
%model.toleranceMB = 2e-7;
model.dsMaxAbs = 0.05;
model = model.validateModel();

if opt.compiled == 1
    % Set up a compiled linear solver
    % We use a AMGCL-based CPR solver to solve the problem. It is sufficient to
    % have a working C++ compiler, the AMGCL repository and BOOST available to
    % compile it. See the documentation of amgcl_matlab for more details on how
    % to set up these paths.
    %solver  = AMGCL_CPRSolverAD('tolerance', 1e-3, 'block_size', ncomp, ...
    %                            'useSYMRCMOrdering', true, ...
    %                            'coarsening', 'aggregation', 'relaxation', 'ilu0',...
    %                           'Solver', 'bicgstab');
    %nls     = NonLinearSolver('LinearSolver', solver);
    
    %nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
    %stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, 'targetChangeAbs', 0.25);
    %nls.timeStepSelector = stepSel;
    %nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');
    nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
                             'useCPR', true);
    nls.LinearSolver.maxIterations = 25;
    %nls.useRelaxation = 1;
    %nls.relaxationType = 'sor';
    %nls.enforceResidualDecrease = 1;
    nls.useLinesearch = 1;
    
%      nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-3, 'Solver', 'bicgstab', ...
%                                                'preconditioner', 'relaxation', 'coarsening', 'aggregation', ...
%                                                'relaxation', 'ilu0', 'maxIterations', 200, 'block_size', ncomp);
%     nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-3, 'Solver', 'gmres', ...
%                                               'preconditioner', 'amg', 'coarsening', 'smoothed_aggregation');
    nls.maxIterations = 12;
    nls.maxTimestepCuts = 8;
    %nls.acceptanceFactor = 2;
else
    nls     = [];
end

end