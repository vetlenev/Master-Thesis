%% Hybrid model for diffuse leakage through thin horizontal shale
% Thin shale layer is either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layer is refined to fine-scale cells.

gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui test-suite

rootdir = strrep(ROOTDIR, '\', '/');
data_dir = strcat(rootdir, '../Master-Thesis/book_illustrations/hybrid3D/test/data');
mkdir(data_dir);

my_seed = 2234;
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

pe_sealing = 4*10^4*Pascal; %5*10^4*Pascal; % higher entry pressure in aquitards
pe_rest = 5*10^3*Pascal;%5*10^3*Pascal;
p_cap = 5*pe_sealing;

% Shut off ADI warning:
warning('off', 'Future:Deprecation');

%% Setup original fine-scale case
% Using face constraints, fine cells only constitute the well
% Using cell constraints, fine cells also include sealing layer
useFaceConstraint = false;
useAdaptive = false; % overrides useFaceConstraint in setupSlopedGrid
run3D = false;

trans_mult = 1e-5; % 1e-4
trans_mult = ~useAdaptive*(useFaceConstraint*trans_mult + ~useFaceConstraint) ...
                + useAdaptive*trans_mult;

if run3D
    nx = 100; ny = 8; nz = 50;
	lx = 500; ly = 25; lz = 250;
    [state0, models, schedule, ...
    isFineCells, sealingFaces] = setupSlopedGrid3D(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    vx = 30; vz = 20; % view angles
    hybrid_folder = 'hybrid3D';
else
    nx = 200; ny = 1; nz = 50;
    lx = 750; ly = 1; lz = 250;
    [state0, models, schedule, ...
    isFineCells, sealingFaces] = setupHorizontalGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    vx = 0; vz = 0;
    hybrid_folder = 'hybrid2D';
end

if useFaceConstraint
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/book_illustrations/%s/test/horizontal/face_lowperm_linrelperm/'), hybrid_folder);   
else
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/book_illustrations/%s/test/horizontal/cell_lowperm_linrelperm/'), hybrid_folder);
end
mkdir(plot_dir);

sealingCells = isFineCells.sealingCells;
sealingCells = any(cell2mat(sealingCells), 2);
% -----
sealingCells_faces = false(size(sealingFaces));
sealingCells_faces(isFineCells.sealingCells_faces) = true;
% -----
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

%% Plot fine grid
figure(1);
%plotOutlinedGrid(G, W, bc, model.operators.T_all);
plotWell(G, W, 'color', 'r')
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')
%plotFaceData(G, model.operators.T_all);

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting
%fineCellsIdx = G.faces.neighbors(sealingCells_faces);
%fineCellsIdx = fineCellsIdx(fineCellsIdx > 0)

plotCellData(G, double(fineCells), fineCellsIdx, 'facecolor', 'green'); % requires numeric input, not logical
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title({'Fine scale grid', 'Imposed fine cells highlighted in green'})

%% More plots
figure(2);
plotOutlinedGrid(G, W, bc, sealingFaces);
plotOutlinedGrid(G, W, bc, sealingCells_faces);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting

plotCellData(G, log10(convertTo(model.rock.perm, milli*darcy))); % requires numeric input, not logical
colorbar('southoutside')
if ~useFaceConstraint
    caxis([min(log10(convertTo(model.rock.perm, milli*darcy))), ...
            max(log10(convertTo(model.rock.perm, milli*darcy)))])
end
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title('Permeability field')

% Plot van Genuchten capillary pressure
swMin = 0.2; swr_pc = 0.1; snr_pc = 0.15;
alpha = -1/(10^4*Pascal); n = 1.5; gamma_D = 2; gamma_I = 1.5;
pc_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.pc_func(Sw, sne, swr, alpha, n, gamma);
sw_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.sw_func(Sw, sne, swr, alpha, n, gamma);
pc_scan = Hysteresis.Genuchten(dummy_s, swMin, snr_pc, swr_pc, pc_func, sw_func, alpha, n, gamma_D, gamma_I);

% fig3 = figure(3);
% plot(1-dummy_s, Capillary.PcGas(dummy_s, swr, snr, pe_rest, 4), 'LineWidth', 1.5);
% xlabel('Water saturation');
% xlim([swr, 1]);
% title('Capillary pressure function');
% saveas(fig3, strcat(plot_dir, '/cap_pres'), 'png');
% hold off

%% Nonlinear solver
nls = NonLinearSolver('maxIterations', 70);

%% Pack and simulate problem
problem = packSimulationProblem(state0, model_fine, schedule, ...
    'horz_finescale', 'NonLinearSolver', nls);

% Simulate and get the output
simulatePackedProblem(problem);
[ws, states, report] = getPackedSimulatorOutput(problem);
%[ws, states] = simulateScheduleAD(state0, model_fine, schedule, 'NonLinearSolver', nls);

%% Simple VE.
% Simulate standard VE model, not accounting for sealing faces.
model_ve = convertToMultiVEModel_res(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'horz_ve', 'NonLinearSolver', nls);

%[ok, status] = simulatePackedProblem(problem_ve);
%[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);

[ws_ve, states_ve] = simulateScheduleAD(state0_ve, model_ve, schedule_ve, 'NonLinearSolver', nls);

states_ve_fs = convertMultiVEStates_test(model_ve, states_ve); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
[model_hybrid, model_coarse] = convertToMultiVEModel_res(model, fineCells, ...
                                            'sealingFaces', find(sealingFaces), ... % same as find(model.operators.T_all == 0)
                                            'sealingCells', sealingCells, ...
                                            'sealingCells_faces', sealingCells_faces, ...
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
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
title('Grid for coarse model')

% Plot grid 
f3 = figure(3);
plotOutlinedGrid(G, W, bc, sealingFaces | sealingCells_faces);
%plotGrid(model_hybrid.G, 'edgealpha', 0.5, 'edgecolor', 'green', 'facecolor', 'blue')
%plotFaceData(model_hybrid.G, model_hybrid.operators.T_all);
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
title('Grid for hybrid VE model')
saveas(f3, strcat(plot_dir, 'sealing_faces'), 'png')

% Plot VE partitions on grid
% Here we can see how VE cells exist between sealing layers and can be
% stacked on top of each other.
f4 = figure(4);
plotOutlinedGrid(G, W, bc, sealingFaces | sealingCells_faces);
plotCellData(G, model_hybrid.G.partition,'edgealpha', 0.2)
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title('Partition of hybrid VE model')
saveas(f4, strcat(plot_dir, 'hybrid_discretization'), 'png')

figure(5)
vtcA = model_hybrid.operators.connections.veTransitionHorizontalConn;
vtcB = model_hybrid.operators.connections.veToFineVertical | model_hybrid.operators.connections.veTransitionVerticalConn;
vtc_cA1 = model_hybrid.operators.N(vtcA, 1);
vtc_cA2 = model_hybrid.operators.N(vtcA, 2);
vtc_cB1 = model_hybrid.operators.N(vtcB, 1);
vtc_cB2 = model_hybrid.operators.N(vtcB, 2);
vtc_A = union(vtc_cA1, vtc_cA2);
vtc_B = union(vtc_cB1, vtc_cB2);
vtc = intersect(vtc_A, vtc_B);

plotGrid(model_hybrid.G, 'edgealpha', 0.2, 'facecolor', 'none')
all_coarse_cells = zeros(model_hybrid.G.cells.num, 1);
all_coarse_cells(vtc_cB1) = 1;
all_coarse_cells(vtc_cB2) = 1;
%all_coarse_cells(1:200) = 1;
plotCellData(model_hybrid.G, all_coarse_cells)
%plotCellData(model_hybrid.G, model_hybrid.G.cells.discretization);
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])
title('Combined veToFineVertical and veHorizontal transition regions')

%% Set capillary state function
nls.useRelaxation = true;

model_hybrid = model_hybrid.validateModel;
model_hybrid.FlowPropertyFunctions.CapillaryPressure = HybridCapillaryPressure(model_hybrid);

%% Simulate Hybrid model
problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
    'horz_hybrid', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_hybrid);

[ws_hybrid, states_hybrid, report_hybrid] = getPackedSimulatorOutput(problem_hybrid);

%% Simulate schedule hybrid
[ws_hybrid, states_hybrid] = simulateScheduleAD(state0_hybrid, model_hybrid, schedule_hybrid, 'NonLinearSolver', nls);

%% Reconstruct fine states
states_hybrid_fs = convertMultiVEStates_res(model_hybrid, states_hybrid);
%states_hybrid_fs = convertMultiVEStates_test(model_hybrid, states_hybrid);

%% Compute net CO2 volume in each VE column 
% and compare with fine solution
sn_f = states{end}.s(:,2); % fine saturations
sn_h = states_hybrid{end}.s(:,2); % hybrid saturations
p = model_hybrid.G.partition;

[vol_unique, ~, idx_sort] = unique(model_hybrid.G.cells.volumes);

pvh = poreVolume(model_hybrid.G, model_hybrid.rock);
pvf = poreVolume(model.G, model.rock);

sn_f_net = accumarray(p, sn_f.*pvf);
sn_f_net = sn_f_net ./ pvh;

% VOLUME MISMATCH FOR COARSE CELLS
diff_fh = zeros(numel(vol_unique), 1);
var_fh = zeros(size(diff_fh));
for i=1:numel(vol_unique)    
    iv = idx_sort == i;
    diff_fh(i) = mean(abs(sn_h(iv) - sn_f_net(iv)));    
    var_fh(i) = var(abs(sn_h(iv) - sn_f_net(iv)));
end

f10 = figure(10);
plot(1:numel(diff_fh), diff_fh, 'b', 'DisplayName', 'Mean')
hold on
plot(1:numel(var_fh), var_fh, '-r', 'DisplayName', 'Variance')
% SET VOLUME TICKS !
xlabel('Coarse cells (increasing volume)')
title('Absolute difference in CO2 saturation: coarse')
legend();
drawnow;
saveas(f10, strcat(plot_dir, 'diff_sat_coarse'), 'png')

% VOLUME MISMATCH FOR FINE CELLS IN DIFFERENT DISCRETIZATION REGIONS
sn_hf = states_hybrid_fs{end}.s(:,2);
discr = model_hybrid.G.cells.discretization;
discr_u = unique(discr);

diff_fh = zeros(numel(discr_u), 1);
var_fh = zeros(size(diff_fh));
for i=1:numel(discr_u) % discretization starts at 1 for fine cells
    d = discr_u(i);
    dp = discr(p);
    iv = dp == d;
    diff_fh(i) = mean(abs(sn_hf(iv) - sn_f(iv)));    
    var_fh(i) = var(abs(sn_hf(iv) - sn_f(iv)));
end

f11 = figure(11);
plot(1:numel(diff_fh), diff_fh, 'b', 'DisplayName', 'Mean')
hold on
plot(1:numel(var_fh), var_fh, '-r', 'DisplayName', 'Variance')
% SET VOLUME TICKS !
xticklabels(discr_u);
xlabel('Discretization region')
title('Absolute difference in CO2 saturation: reconstruction')
legend();
drawnow;
saveas(f11, strcat(plot_dir, 'diff_sat_fine'), 'png')

%% Plot flux boundary
figure(20);
bfaces = zeros(G.faces.num, 1);
bf = boundaryFaces(G);
bfx = G.faces.centroids(bf, 1);
bfz = G.faces.centroids(bf, 3);
lz = max(bfz(bfx == min(bfx)));
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
top_left_faces = bf(bfx == min(bfx) & bfz <= lz/nz);
%plotOutlinedGrid(G, W, bc, top_left_faces);
nn = model_hybrid.G.faces.neighbors;
nn = find(~all(nn,2));
plotGrid(model_hybrid.G, 'faceColor', 'none', 'edgecolor', 'none');
plotFaces(model_hybrid.G, nn, 'r')
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])

%% Plot CO2 saturation for each model
fafa = find(sealingFaces | sealingCells_faces);

c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

for i = 1:20:numel(states)
    f1 = figure(1); clf
    plotCellData(model.G, states{i}.s(:,2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.3)
    view(vx, vz); colormap(cc); caxis([0, 1-swr]);
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
    f2 = figure(2); clf    
    plotCellData(model.G, states_hybrid_fs{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.3)
    view(vx, vz); colormap(cc); caxis([0, 1-swr]); 
    axis tight off
    title(['Hybrid: VE reconstructed saturation, step:', num2str(i)])
    %pause(0.5);
    saveas(f2, sprintf(strcat(plot_dir, 'hybrid_sat_%d'), i), 'png');  
    
    if 0
        f2 = figure(2); clf    
        seal = ismember((1:model_hybrid.G.cells.num)', model_hybrid.G.sealingCells);
        ve = model_hybrid.G.cells.discretization > 1;
        plotCellData(model.G, model_hybrid.fluid.pcWG(states_hybrid_fs{i}.s(:,2), seal, ve), 'edgec', 'none');
        plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2, 'facealpha', 0)
        view(0, 0); colormap(cc); colorbar('location', 'southoutside');
        axis tight off
        title(['Hybrid capillary pressure, step:', num2str(i)])   

        saveas(f2, sprintf(strcat(plot_dir, 'hybrid_pc_%d'), i), 'png');
    end
end
%% Plot 2D section
if run3D    
   section_2d = (jj == 6);  frf
   
   cells_2d = model.G.cells.indexMap(section_2d);
   n_2d = unique(model.G.faces.neighbors(fafa, :));
   cells_2d_faces = cells_2d(ismember(cells_2d, n_2d));     
   faces_2d = gridCellFaces(model.G, cells_2d_faces);
   faces_2d = fafa(ismember(fafa, faces_2d));
   
   for i=1:20:fix(2*numel(states)/3)
       f6 = figure(6); clf
       sf = states{i}.s(:,2);      
       sf = sf(section_2d);
       plotCellData(model.G, sf, cells_2d, 'edgec', 'none')
       plotFaces(G, faces_2d, 'facec', 'w', 'linewidth', 0.5, 'facealpha', 0)
       view(0, 0); colormap(cc); caxis([0,1-swr]);
       colorbar('location', 'southoutside');
       axis tight on    
       pause(0.2)
       title(['Fine saturation, step:', num2str(i)])
       
       saveas(f6, sprintf(strcat(plot_dir, 'jj6_fine_sat_%d'), i), 'png');
   end
   
   for i=1:20:fix(2*numel(states_hybrid)/3)
       f7 = figure(7); clf
       sh = states_hybrid_fs{i}.s(:,2);      
       sh = sh(section_2d);
       plotCellData(model.G, sh, cells_2d, 'edgec', 'none')
       plotFaces(G, faces_2d, 'facec', 'w', 'linewidth', 0.5, 'facealpha', 0)
       view(0, 0); colormap(cc); caxis([0,1-swr]);
       colorbar('location', 'southoutside');
       axis tight on
       pause(0.2)
       title(['Hybrid saturation, step:', num2str(i)])
       
       saveas(f7, sprintf(strcat(plot_dir, 'jj6_hybrid_sat_%d'), i), 'png');
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

fign = 12;
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
        %plotOutlinedGrid(G, W, bc, sealingCells_faces, 'facecolor', 'none');
        if ~isfield(g, 'parent')
            plotFaces(g, fafa, 'facec', 'w', 'facealpha', 0, 'linewidth', 0.3);
        end

        view(vx, vz)
        axis equal tight
        if run3D
            daspect([1, 0.1, 1]) 
        else
            daspect([1, 0.1, 0.5])
        end
        %xlim([xl_min, xl_max]);
        %zlim([zl_min, zl_max]);
        colormap(cc); caxis([0,1-swr]);
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