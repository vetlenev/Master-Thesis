%% Hybrid model for diffuse leakage through thin horizontal shale
% Thin shale layer is either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layer is refined to fine-scale cells.

gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui test-suite

rootdir = strrep(ROOTDIR, '\', '/');
data_dir = strcat(rootdir, '../master_thesis/hybrid2D/sloped/data');
mkdir(data_dir);

my_seed = 2233;
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Setup original fine-scale case
% The setupHorizontalGrid function generates the fine-scale grid and
% computes model, schedule and initial state, and adjusts model based on 
% given contraint

% Using face constraints, fine cells only constitute the well
% Using cell constraints, fine cells also include sealing layer
useFaceConstraint = false;

if useFaceConstraint
    plot_dir = strcat(rootdir, '../master_thesis/hybrid2D/sloped/figs/face_trans_low/');   
else
    plot_dir = strcat(rootdir, '../master_thesis/hybrid2D/sloped/figs/cell_trans_low/');
end
mkdir(plot_dir);

nx = 200; ny = 1; nz = 50;
lx = 750; ly = 1; lz = 250;
trans_mult = useFaceConstraint*1e-7 + ~useFaceConstraint; % 1e-4

[state0, models, schedule, isFine, facesZeroTrans] = setupSlopedGrid(useFaceConstraint, [nx,ny,nz], [lx,ly,lz], trans_mult);

sealingCells = isFine.sealing;
wellCells = isFine.well;
fineCells = sealingCells | wellCells;

model = models.original; % used as input to VE models
model_fine = models.fine;

G = model_fine.G;
W = schedule.control(1).W;
bc = schedule.control(1).bc;
T_all = model.operators.T_all;
T = model.operators.T;
[ii, jj, kk] = gridLogicalIndices(G);
nn = G.faces.neighbors;

% Plot fine grid
figure(1);
plotOutlinedGrid(G, W, bc, model.operators.T_all);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')
%plotFaceData(G, model.operators.T_all);

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

plotCellData(G, log10(convertTo(model.rock.perm, milli*darcy))); % requires numeric input, not logical
colorbar('southoutside')
if ~useFaceConstraint
    caxis([min(log10(convertTo(model.rock.perm, milli*darcy))), ...
            max(log10(convertTo(model.rock.perm, milli*darcy)))])
end
view(0, 0)
axis equal tight
daspect([1, 0.1, 0.5])
xlabel('Lateral position [m]');
ylabel('Depth [m]');
title('Permeability field')

%% Nonlinear solver
nls = NonLinearSolver('maxIterations', 50);

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
                                            'multiplier', trans_mult, ...
                                            'sumTrans', true);
                                        
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
%plotFaceData(model_hybrid.G, model_hybrid.operators.T_all);
view(0, 0)
axis equal tight
daspect([1, 0.1, 0.5])
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

%% Simulate Hybrid model
nls.useRelaxation = true;
problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
    'horz_hybrid', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_hybrid);

[ws_hybrid, states_hybrid, report_hybrid] = getPackedSimulatorOutput(problem_hybrid);
states_hybrid_fs = convertMultiVEStates(model_hybrid, states_hybrid);

%% Plot CO2 saturation for each model
fafa = find(facesZeroTrans);

for i = 1:10:numel(states)
    f1 = figure(1); clf    
    plotCellData(model.G, states{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0); colormap(winter);
    axis tight off
    title(['Fine-scale saturation, step:', num2str(i)])   
    
    saveas(f1, sprintf(strcat(plot_dir, 'fine_sat_%d'), i), 'png');  
end

for i = 1:10:numel(states)
    f1 = figure(1); clf    
    plotCellData(model.G, states_hybrid_fs{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    colormap(winter);
    axis tight off
    title(['Hybrid: VE reconstructed saturation, step:', num2str(i)])
    
    saveas(f1, sprintf(strcat(plot_dir, 'hybrid_sat_%d'), i), 'png');  
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
xl_min = -max(G.cells.centroids(:,1))/50;
xl_max = max(G.cells.centroids(:,1))*(1+1/50);
zl_min = -max(G.cells.centroids(:,3))/50;
zl_max = max(G.cells.centroids(:,3))*(1+1/50);

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
            nm = 'Hybrid VE reconstructed';
        elseif j == 3
            g = model_hybrid.G;
            s = states_hybrid{ss}.s(:, 2);
            nm = 'Coarse';
        else
            g = G;
            s = states_ve_fs{ss}.s(:, 2);
            nm = 'VE reconstructed (no layers)';
        end
        subplot(1, 2, 2-mod(j,2));
        plotCellData(g, s, 'EdgeColor', 'none')
        plotOutlinedGrid(G, W, bc, facesZeroTrans);

        view(0, 0)
        axis equal tight
        daspect([1, 1, 0.5]) 
        %xlim([xl_min, xl_max]);
        %zlim([zl_min, zl_max]);
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