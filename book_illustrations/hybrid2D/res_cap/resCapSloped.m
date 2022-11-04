%% Hybrid model for diffuse leakage through thin horizontal shale
% Thin shale layer is either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layer is refined to fine-scale cells.

gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui test-suite

rootdir = strrep(ROOTDIR, '\', '/');
data_dir = strcat(rootdir, '../Master-Thesis/book_illustrations/hybrid2D/res_cap/data');
mkdir(data_dir);

my_seed = 2234;
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

pe_sealing = 1*10^5*Pascal; %5*10^4*Pascal; % higher entry pressure in aquitards
pe_rest = 5*10^3*Pascal;%5*10^3*Pascal;
p_cap = 5*pe_sealing;

%% Setup original fine-scale case
% Using face constraints, fine cells only constitute the well
% Using cell constraints, fine cells also include sealing layer
useFaceConstraint = false;

if useFaceConstraint
    plot_dir = strcat(rootdir, '../Master-Thesis/book_illustrations/hybrid2D/res_cap/figs/face_orig_large/');   
else
    plot_dir = strcat(rootdir, '../Master-Thesis/book_illustrations/hybrid2D/res_cap/figs/cell_orig_large/');
end
mkdir(plot_dir);

nx = 200; ny = 1; nz = 50;
lx = 750; ly = 1; lz = 250;
trans_mult = useFaceConstraint*1e-4 + ~useFaceConstraint; % 1e-4

useAdaptive = false;

[state0, models, schedule, ...
    isFineCells, facesZeroTrans] = setupSlopedGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);

sealingCells = isFineCells.sealing;
sealingCells = any(cell2mat(sealingCells), 2);
wellCells = isFineCells.well;
fineCells = sealingCells | wellCells;

% Split into two models: only fine model has transmissibility multiplier
% applied
model = models.original; % used as input to VE models
model_fine = models.fine;

swr = model_fine.fluid.krPts.w(1); % reisdual water sat
snr = model_fine.fluid.krPts.g(1); % residual CO2 sat
dummy_s = linspace(0, 1, model_fine.G.cells.num)';

sealingCellsFine = model_fine.G.cells.indexMap(sealingCells);
model_fine.fluid.pcWG = @(s) Capillary.runStandardPcSharp(s, dummy_s, swr, snr, pe_sealing, pe_rest, ...
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
title({'Fine scale grid', 'Imposed fine cells highlighted in green'})

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

% Plot van Genuchten capillary pressure
swMin = 0.2; swr = 0.1; snr = 0.15;
alpha = -1/(10^4*Pascal); n = 1.5; gamma_D = 2; gamma_I = 1.5;
pc_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.pc_func(Sw, sne, swr, alpha, n, gamma);
sw_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.sw_func(Sw, sne, swr, alpha, n, gamma);
pc_scan = Hysteresis.Genuchten(dummy_s, swMin, snr, swr, pc_func, sw_func, alpha, n, gamma_D, gamma_I);

fig3 = figure(3);
plot(1-dummy_s, Capillary.PcGas(dummy_s, swr, snr, pe_rest, 4), 'LineWidth', 1.5);
xlabel('Water saturation');
xlim([swr, 1]);
title('Capillary pressure function');
saveas(fig3, strcat(plot_dir, '/cap_pres'), 'png');
hold off

%% Nonlinear solver
nls = NonLinearSolver('maxIterations', 70);

%% Pack and simulate problem
problem = packSimulationProblem(state0, model_fine, schedule, ...
    'horz_finescale', 'NonLinearSolver', nls);

%[states, wellSols] = simulateScheduleAD(state0, model_fine, schedule, 'NonLinearSolver', nls);
% Simulate and get the output
simulatePackedProblem(problem);
[ws, states, report] = getPackedSimulatorOutput(problem);

%% Simple VE.
% Simulate standard VE model, not accounting for sealing faces.
model_ve = convertToMultiVEModel_res(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'horz_ve', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_ve);

[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);
states_ve_fs = convertMultiVEStates_res(model_ve, states_ve); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
[model_hybrid, model_coarse] = convertToMultiVEModel_res(model, fineCells, ...
                                            'sealingFaces', find(facesZeroTrans), ... % same as find(model.operators.T_all == 0)
                                            'sealingCells', sealingCells, ...
                                            'multiplier', trans_mult, ...
                                            'sumTrans', true, ...
                                            'pe_rest', pe_rest);
                                       
model_hybrid.fluid.pcWG = @(s, n_sealing) ...
                            Capillary.runStandardPcSharp(s, dummy_s, swr, snr, pe_sealing, ...
                                                pe_rest, n_sealing, model_hybrid.G);                                                                                                      
                                        
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

figure(5)
vtc = model_hybrid.operators.connections.veToFineHorizontal;
vtc_c1 = model_hybrid.operators.N(vtc, 1);
vtc_c2 = model_hybrid.operators.N(vtc, 2);
plotGrid(model_hybrid.G, 'edgealpha', 0.2, 'facecolor', 'none')
all_coarse_cells = zeros(model_hybrid.G.cells.num, 1);
all_coarse_cells(vtc_c1) = 1;
all_coarse_cells(vtc_c2) = 1;
%all_coarse_cells(1:33) = 1;
plotCellData(model_hybrid.G, all_coarse_cells)
%plotCellData(model_hybrid.G, model_hybrid.G.cells.discretization);
view(0, 0)
axis equal tight
daspect([1, 1, 0.5])
title('VE to Fine transition regions')

%% Set capillary state function
nls.useRelaxation = true;

model_hybrid = model_hybrid.validateModel;
model_hybrid.FlowPropertyFunctions.CapillaryPressure = HybridCapillaryPressure(model_hybrid);

%% Simulate Hybrid model
problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
    'horz_hybrid', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_hybrid);

[ws_hybrid, states_hybrid, report_hybrid] = getPackedSimulatorOutput(problem_hybrid);

%% Reconstruct fine states
states_hybrid_fs = convertMultiVEStates_res(model_hybrid, states_hybrid);

%% Simulate schedule hybrid
[states, wellSols] = simulateScheduleAD(state0_hybrid, model_hybrid, schedule_hybrid, 'NonLinearSolver', nls);

%% Compute net CO2 volume in each VE column 
% and compare with fine solution
sn_f = states{end}.s(:,2); % fine saturations
sn_h = states_hybrid{end}.s(:,2); % hybrid saturations
p = model_hybrid.G.partition;
[vol_sort, idx_sort] = sort(model_hybrid.G.cells.volumes);

idx_new_vol = find(diff(vol_sort) > 1e-5);
idx_new_vol = [1; idx_new_vol];
idx_new_vol = [idx_new_vol; numel(vol_sort)];

pvh = poreVolume(model_hybrid.G, model_hybrid.rock);
pvf = poreVolume(model.G, model.rock);

%sn_f_net = accumarray(p, sn_f, size(unique(p)), @mean);
sn_f_net = accumarray(p, sn_f.*pvf, size(unique(p)));
sn_f_net = sn_f_net ./ pvh;
sn_f_net = sn_f_net(idx_sort);
%sn_h = sn_h.*pvh;
sn_h = sn_h(idx_sort);

%sn_f_net = sn_f;
diff_fh = zeros(numel(idx_new_vol)-1, 1);
var_fh = zeros(size(diff_fh));
for i=1:numel(idx_new_vol)-1
    ii = idx_new_vol(i);
    jj = idx_new_vol(i+1);
    diff_fh(i) = median(abs(sn_h(ii:jj) - sn_f_net(ii:jj)));
    var_fh(i) = var(abs(sn_h(ii:jj) - sn_f_net(ii:jj)));
end
%diff_fh = abs(sn_h - sn_f_net);

figure();
% Smallest volumes are fine cells, which obviously gives accurate
% reconstruction of saturation
plot(1:numel(diff_fh), diff_fh, 'b', 'DisplayName', 'Median')
hold on
plot(1:numel(var_fh), var_fh, '-r', 'DisplayName', 'Variance')
% SET VOLUME TICKS !
xlabel('Coarse cells (increasing volume)')
title('Absolute difference in CO2 saturation for hybrid and fine models')
legend();
drawnow;

%% Plot CO2 saturation for each model
fafa = find(facesZeroTrans);

c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

for i = 1:20:numel(states)
    f1 = figure(1); clf    
    plotCellData(model.G, states{i}.s(:,2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0); colormap(cc);
    axis tight off
    title(['Fine-scale saturation, step:', num2str(i)])   
    
    saveas(f1, sprintf(strcat(plot_dir, 'fine_sat_%d'), i), 'png'); 
    
    if 0
        f2 = figure(2); clf    
        plotCellData(model.G, model_fine.fluid.pcWG(states{i}.s(:,2)), 'edgec', 'none');
        plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
        view(0, 0); colormap(cc); 
        axis tight off
        title(['Fine-scale capillary pressure, step:', num2str(i)])   

        saveas(f2, sprintf(strcat(plot_dir, 'fine_pc_%d'), i), 'png');
    end
end

for i = 1:20:numel(states_hybrid)
    f1 = figure(1); clf    
    plotCellData(model.G, states_hybrid_fs{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    colormap(cc); caxis([0, 1]);
    axis tight off
    title(['Hybrid: VE reconstructed saturation, step:', num2str(i)])
    
    saveas(f1, sprintf(strcat(plot_dir, 'hybrid_sat_%d'), i), 'png');  
    
    if 0
        f2 = figure(2); clf    
        seal = ismember((1:model_hybrid.G.cells.num)', model_hybrid.G.sealingCells);
        ve = model_hybrid.G.cells.discretization > 1;
        plotCellData(model.G, model_hybrid.fluid.pcWG(states_hybrid_fs{i}.s(:,2), seal, ve), 'edgec', 'none');
        plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
        view(0, 0); colormap(cc); colorbar('location', 'southoutside');
        axis tight off
        title(['Hybrid capillary pressure, step:', num2str(i)])   

        saveas(f2, sprintf(strcat(plot_dir, 'hybrid_pc_%d'), i), 'png');
    end
end

%% Plot result at injection stop and end of simulation
nstep = numel(schedule.step.val);
end_inj = find(schedule.step.control == 1, 1, 'last');
end_mig = nstep;

substeps = [end_inj, end_mig];
%substeps = [nstep/10, nstep/5];
names = {'injected', 'migrated'};
model_names = {'fine_hybrid', 'coarse_ve'};

xl_min = -max(G.cells.centroids(:,1))/50;
xl_max = max(G.cells.centroids(:,1))*(1+1/50);
zl_min = -max(G.cells.centroids(:,3))/50;
zl_max = max(G.cells.centroids(:,3))*(1+1/50);

fign = 6;
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