%% Hybrid model for diffuse leakage through thin horizontal shale
% Thin shale layer is either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layer is refined to fine-scale cells.

gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui test-suite

rootdir = strrep(ROOTDIR, '\', '/');
data_dir = strcat(rootdir, '../master_thesis/hybrid2D/data');
mkdir(data_dir);

my_seed = 5577;
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Setup simulation prooperties
% Using face constraints, fine cells only constitute the well
% Using cell constraints, fine cells also include sealing layer
useFaceConstraint = false;

if useFaceConstraint
    plot_dir = strcat(rootdir, '../Master-Thesis/hybrid2D/figs/face_constraint_pc/');   
else
    plot_dir = strcat(rootdir, '../Master-Thesis/hybrid2D/figs/cell_constraint_pc/');
end
mkdir(plot_dir);

nx = 200; ny = 1; nz = 50;
lx = 800; ly = 1; lz = 250;
trans_mult = useFaceConstraint*1e-4 + ~useFaceConstraint; % 1e-5
pe_sealing = 1*barsa; 
pe_rest = 0.05*barsa;
p_cap = 5*barsa;

%% Setup original fine-scale case
% The setupHorizontalGrid function generates the fine-scale grid and
% computes model, schedule and initial state, and adjusts model based on 
% given contraint

[state0, models, schedule, isFine, facesZeroTrans] = setupHorizontalGrid(useFaceConstraint, [nx,ny,nz], [lx,ly,lz], trans_mult);

sealingCells = isFine.sealing;
wellCells = isFine.well;
fineCells = sealingCells | wellCells;

model = models.original; % used as input to VE models
model_fine = models.fine;

swr = model_fine.fluid.krPts.w(1); % reisdual water sat
snr = model_fine.fluid.krPts.g(1); % residual CO2 sat
dummy_s = linspace(0, 1, 1000)';

sealingCellsFine = model_fine.G.cells.indexMap(sealingCells);
model_fine.fluid.pcWG = @(s) Capillary.runStandardPc(s, dummy_s, swr, snr, pe_sealing, pe_rest, p_cap, ...
                                            sealingCellsFine, model_fine.G);

G = model_fine.G;
W = schedule.control(1).W;
bc = schedule.control(1).bc;
T_all = model.operators.T_all;
T = model.operators.T;
[ii, jj, kk] = gridLogicalIndices(G);
nn = G.faces.neighbors;

% Plot fine grid
figure(1);
%plotOutlinedGrid(G, W, bc, facesZeroTrans);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')
plotFaces(G, facesZeroTrans, 'edgecolor', 'red', 'linewidth', 2);

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting

plotCellData(G, double(fineCells), fineCellsIdx, 'facecolor', 'green'); % requires numeric input, not logical
view(0, 0)
axis equal tight
daspect([1, 0.1, 0.5])
xlabel('Lateral position [m]');
ylabel('Depth [m]');
title({'Fine scale grid with flat sealing layer', 'Imposed fine cells highlighted in green'})

figure(2);
plotOutlinedGrid(G, W, bc, facesZeroTrans);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting

plotCellData(G, model.rock.perm); % requires numeric input, not logical
view(0, 0)
axis equal tight
daspect([1, 0.1, 0.5])
xlabel('Lateral position [m]');
ylabel('Depth [m]');
title('Permeability field')

%% Nonlinear solver
nls = NonLinearSolver('maxIterations', 70);

%% Pack and simulate problem
problem = packSimulationProblem(state0, model_fine, schedule, ...
    'horz_finescale', 'NonLinearSolver', nls);
 
% Simulate and get the output
simulatePackedProblem(problem);
[ws, states, report] = getPackedSimulatorOutput(problem);

%% Simple VE.
% Simulate standard VE model, not accounting for sealing faces.
model_ve = convertToMultiVEModel_test(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'horz_ve', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_ve);

[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);
states_ve_fs = convertMultiVEStates(model_ve, states_ve); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
[model_hybrid, model_coarse] = convertToMultiVEModel_test(model, fineCells, ...
                                            'sealingFaces', find(facesZeroTrans), ... % same as find(model.operators.T_all == 0)
                                            'sealingCells', sealingCells, ...
                                            'multiplier', trans_mult, ...
                                            'sumTrans', true);
                                        
%sealingCellsHybrid = model_hybrid.G.partition(sealingCells);
sealingCellsHybrid = model_hybrid.G.sealingCells;

model_hybrid.fluid.pcWG = @(s, n_sealing, isVE) Capillary.runHybridPc(s, swr, snr, pe_sealing, pe_rest, p_cap, n_sealing, isVE);                                      
                                        
schedule_hybrid = upscaleSchedule(model_hybrid, schedule);
state0_hybrid = upscaleState(model_hybrid, model, state0);

% Plot coarse grid
figure(2);
%plotOutlinedGrid(model_coarse.G, W, bc, T_all==0);
plotGrid(model_coarse.G, 'edgealpha', 0.2, 'facecolor', 'none')
view(0, 0)
axis equal tight
daspect([1, 1, 0.5])
title('Grid for coarse model')

% Plot grid 
figure(3);
plotOutlinedGrid(G, W, bc, facesZeroTrans);
plotGrid(model_hybrid.G, 'edgealpha', 0.2, 'facecolor', 'none')
idx = zeros(model_hybrid.G.cells.num, 1);
idx(sealingCellsHybrid) = 1;
plotCellData(model_hybrid.G, idx, 'edgealpha', 0.2)
view(0, 0)
axis equal tight
daspect([1, 1, 0.5])
title('Grid for hybrid VE model')

% Plot VE partitions on grid
% Here we can see how VE cells exist between sealing layers and can be
% stacked on top of each other.
figure(4);
plotOutlinedGrid(G, W, bc, facesZeroTrans);
plotCellData(G, model_hybrid.G.partition,'edgealpha', 0.2)
view(0, 0)
axis equal tight
daspect([1, 1, 0.5])
xlabel('Lateral position [m]');
ylabel('Depth [m]');
title('Partition of hybrid VE model')

figure(5)
vtc = model_hybrid.operators.connections.veTransitionHorizontalConn;
vtc_c1 = model_hybrid.operators.N(vtc, 1);
vtc_c2 = model_hybrid.operators.N(vtc, 2);
plotGrid(model_hybrid.G, 'edgealpha', 0.2, 'facecolor', 'none')
all_coarse_cells = zeros(model_hybrid.G.cells.num, 1);
all_coarse_cells(vtc_c1) = 1;
all_coarse_cells(vtc_c2) = 1;
plotCellData(model_hybrid.G, all_coarse_cells)
view(0, 0)
axis equal tight
daspect([1, 1, 0.5])
title('VE to Fine transition regions')

%% Simulate Hybrid model
nls.useRelaxation = true;

model_hybrid = model_hybrid.validateModel;
model_hybrid.FlowPropertyFunctions.CapillaryPressure = HybridCapillaryPressure(model_hybrid);

problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
    'horz_hybrid', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_hybrid);

[ws_hybrid, states_hybrid, report_hybrid] = getPackedSimulatorOutput(problem_hybrid);
states_hybrid_fs = convertMultiVEStates(model_hybrid, states_hybrid);

%% Separate simulateSchedule
%[ws_hybrid, states_hybrid] = simulateScheduleAD(state0_hybrid, model_hybrid, schedule_hybrid);

%% Plot CO2 saturation for each model
fafa = find(facesZeroTrans);

parg = {2, 1};
for i = 1:20:numel(states)
    f1 = figure(1); clf
    subplot(parg{:}, 1)
    plotCellData(model.G, states{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0); colormap(winter);
    axis tight off
    title(['Fine-scale saturation, step:', num2str(i)])
    
    subplot(parg{:}, 2)
    plotCellData(model.G, states_hybrid_fs{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    colormap(winter);
    axis tight off
    title(['Hybrid: VE reconstructed saturation, step:', num2str(i)])
    
%     subplot(parg{:}, 3)
%     plotCellData(model_hybrid.G, states_hybrid{i}.s(:, 2), 'edgec', 'none')
%     plotGrid(model_hybrid.G, 'facec', 'none', 'edgec', 'w', 'edgea', .3)
%     view(0, 0);
%     colormap(winter);
%     axis tight off
%     title(['Hybrid: Coarse saturation, step:', num2str(i)])
%     drawnow
    
    saveas(f1, sprintf(strcat(plot_dir, 'sat_%d'), i), 'png');  
end

%% Plot result at injection stop and end of simulation
nstep = numel(schedule.step.val);
end_inj = find(schedule.step.control == 1, 1, 'last');
end_mig = nstep;

substeps = [end_inj, end_mig];
%substeps = [nstep/10, nstep/5];
names = {'injected', 'migrated'};
model_names = {'fine_hybrid', 'coarse_ve'};
c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

fign = 5;
for i = 1:numel(substeps)
    ss = substeps(i);
    ff = UtilFunctions.fullsizeFig(fign);
    fign = fign + 1;
    for j = 1:4
        if j == 1
            g = G;
            s = states{ss}.s(:, 2);
            nm = 'Fine-scale';
        elseif j == 2
            g = G;
            s = states_hybrid_fs{ss}.s(:, 2);
            nm = 'Hybrid VE';
        elseif j == 3
            g = model_hybrid.G;
            s = states_hybrid{ss}.s(:, 2);
            nm = 'Coarse';
        else
            g = model_ve.G;
            s = states_ve{ss}.s(:, 2);
            nm = 'VE (no layers)';
        end
        subplot(1, 2, 2-mod(j,2));
        plotCellData(g, s, 'EdgeColor', 'none')
        plotOutlinedGrid(G, W, bc, facesZeroTrans);

        view(0, 0)
        axis equal tight
        daspect([1, 1, 0.5])        
        colormap(cc); caxis([0, 1])
        title(nm);
        
        if ~mod(j,2)
            set(gcf, 'Name', names{i});            
            saveas(ff, strcat(plot_dir, names{i}, '_', model_names{j/2}), 'png');
            ff = UtilFunctions.fullsizeFig(fign);
            fign = fign + 1;
        end
    end
    clf;
end