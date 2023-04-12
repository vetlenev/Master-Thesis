mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui

if false
    name = 'FluidFlowerRs';
    deck = readEclipseDeck('../deck/CSP11A_RS.DATA');
elseif true
    name = 'FluidFlowerImmiscible';
    deck = readEclipseDeck('deck/CSP11A_IMMISCIBLE.DATA');
else
    name = 'FluidFlowerRsRv';
    deck = readEclipseDeck('../deck/CSP11A_RSRV.DATA');
end
deck = convertDeckUnits(deck);
%% eclipse model
[state0, model, schedule, nls] = initEclipseProblemAD(deck);

%%
%model.OutputStateFunctions{end+1} = 'CapillaryPressure';
%model.outputFluxes = false;
% Need to add BC
G = model.G;
bf = boundaryFaces(G);
bf = bf(G.faces.centroids(bf, 3) < 1e-12);
if false
    bc = addBC([], bf, 'pressure', 1*atm, 'sat', [1, 0]);
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end
else
    bcells = sum(G.faces.neighbors(bf, :), 2);
    model.operators.pv(bcells) = model.operators.pv(bcells)*1000;
    state0.pressure(bcells) = 1*atm;
end

%%
nls.LinearSolver = BackslashSolverAD();
nls.maxTimestepCuts = 20;
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
% simulatePackedProblem(problem);
simulatePackedProblem(problem, 'restartStep', 1);
%%
[ws, states] = getPackedSimulatorOutput(problem);
%%
figure; plotToolbar(G, states); title(['MRST ', name])
%%
mrstModule add jutul
fn = writeJutulInput(state0, model, schedule, name);
%%
[wsj, statesj] = readJutulOutput(fn);
figure; plotToolbar(G, statesj); title(['Jutul ', name])

%%
inspectFluidModel(model, 'pressureRange', (0:0.1:2)*barsa)
