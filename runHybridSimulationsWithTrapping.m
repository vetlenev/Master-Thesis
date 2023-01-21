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

pe_sealing = 4*10^4*Pascal; %4*10^4*Pascal; % higher entry pressure in aquitards
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
sloped = true;
standard = false;

trans_mult = 1e-5; % 1e-4
trans_mult = ~useAdaptive*(useFaceConstraint*trans_mult + ~useFaceConstraint) ...
                + useAdaptive*trans_mult;

if run3D
    nx = 60; ny = 10; nz = 40; % 100, 8, 50
	lx = 500; ly = 50; lz = 250; % 500, 25, 250
    [state0, models, schedule, ...
    isFineCells, sealingFaces, allSealingFaces] = setupSlopedGrid3D(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    vx = 30; vz = 20; % view angles
    hybrid_folder = 'hybrid3D';
else
    nx = 75; ny = 1; nz = 50; % nz = 50
    lx = 750; ly = 1; lz = 250;
    if sloped
        [state0, models, schedule, ...
         isFineCells, sealingFaces, allSealingFaces] = setupSlopedGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult, standard);
    else
        [state0, models, schedule, ...
         isFineCells, sealingFaces] = setupHorizontalGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    end
    vx = 0; vz = 0;
    hybrid_folder = 'hybrid2D';
end

if sloped
    geometry_folder = 'caseMultilayered';
else
    geometry_folder = 'caseSimple';
end

if useAdaptive    
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/adaptive_test/'), hybrid_folder, geometry_folder);   
elseif useFaceConstraint
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/face_test/'), hybrid_folder, geometry_folder);       
else    
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/cell_test/'), hybrid_folder, geometry_folder);   
end
mkdir(plot_dir);

sealingCells = isFineCells.sealingCells;
sealingCells = any(cell2mat(sealingCells), 2);
% -----
sealingCells_faces = false(size(sealingFaces));
sealingCells_facesAll = vertcat(isFineCells.sealingCells_faces{:});
sealingCells_faces(sealingCells_facesAll) = true;
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
f1 = figure(1);
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
saveas(f1, strcat(plot_dir, '/fine_regions'), 'png');

%% More plots
figure(2);
plotOutlinedGrid(G, W, bc, sealingFaces);
plotOutlinedGrid(G, W, bc, sealingCells_faces);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting
log10perm = log10(convertTo(model.rock.perm, milli*darcy));

plotCellData(G, log10perm); % requires numeric input, not logical
colorbar('southoutside')
if ~useFaceConstraint
    caxis([min(log10perm)-0.1*min(log10perm), ...
            max(log10perm)+0.1*max(log10perm)])
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
model_ve = convertToMultiVEModel_test(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'horz_ve', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_ve);
[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);

%[ws_ve, states_ve] = simulateScheduleAD(state0_ve, model_ve, schedule_ve, 'NonLinearSolver', nls);

states_ve_fs = convertMultiVEStates_orig(model_ve, model_fine, states_ve, 'schedule', schedule, 'convert_flux', true); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
[model_hybrid, model_coarse] = convertToMultiVEModel_test(model, fineCells, ...
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
Gh = model_hybrid.G;

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
%[ws_hybrid, states_hybrid] = simulateScheduleAD(state0_hybrid, model_hybrid, schedule_hybrid, 'NonLinearSolver', nls);

%% Reconstruct fine states
%states_hybrid_fs = convertMultiVEStates_res(model_hybrid, model_fine, states_hybrid, states);
states_hybrid_fs = convertMultiVEStates_test(model_hybrid, model_fine, states_hybrid, 'schedule', schedule, 'convert_flux', true); % send in fine schedule to extract fine boundary faces

mh = model_hybrid;
oph = mh.operators;
n = oph.N;
shf = states_hybrid_fs;

%% Find traps for global top surface
Gsi = {}; Gfi = {}; Gts = {};
cmaps = {}; fmaps = {}; nmaps = {};
cmapf = {}; fmapf = {}; nmapf = {};
tas = {};

top_cells = G.cells.indexMap(kk == min(kk));
bf = boundaryFaces(G);
top_faces = zeros(numel(top_cells), 6);
for i=1:numel(top_cells)
    top_faces(i,:) = G.cells.faces(G.cells.facePos(top_cells(i)):G.cells.facePos(top_cells(i)+1)-1, 1);
end
top_faces = intersect(bf, top_faces);

[Gs_glob, cmap_glob, fmap_glob, ~] = ExtractLayerSubgrid(G, Gh, top_faces, ...
                                                     isFineCells.sealingCells, allSealingFaces); % cell constraints, face constraints

% OR: Do trap analysis for combined top surface
Gt = topSurfaceGrid(Gs_glob);
ta = trapAnalysis(Gt, true);
Gts = cat(1, Gts, Gt);
tas = cat(1, tas, ta);

%% Subgrids and traps for VE top surfaces

sealingLayers = [isFineCells.sealingBottom, allSealingFaces]; % merge all faces defining top surfaces

for i=1:numel(sealingLayers)
    subfaces = sealingLayers{i}; % bottom faces of sealing layer    
    [Gs, cmap, fmap, nmap] = ExtractLayerSubgrid(G, Gh, subfaces,  ...
                                                isFineCells.sealingCells, allSealingFaces); 
    
    if any(Gs.cells.num) % only append non-zero subgrids
        Gsi = cat(1, Gsi, Gs); % subgrid for VE regions
        Gtsi = topSurfaceGrid(Gs); % top surface of subgrid
        tai = trapAnalysis(Gtsi, true);
        Gts = cat(1, Gts, Gtsi);
        tas = cat(1, tas, tai);
        cmaps = cat(1, cmaps, cmap);
        fmaps = cat(1, fmaps, fmap);
        nmaps = cat(1, nmaps, nmap);
    end    
end

Gs_all = [{Gs_glob}; Gsi];
cmaps_all = [cmap_glob; cmaps];
fmaps_all = [fmap_glob; fmaps];

% Remaining fine cells not subject to trapping
top_surface_cells = vertcat(cmaps_all{:}); % should be unique
fine_rem = setdiff(G.cells.indexMap, top_surface_cells);

%% Plot trapping
figure(15)
%plotCellData(Gt, ones(Gt.cells.num,1), 'EdgeColor', 'none');
for i=1:numel(Gts)
    plotCellData(Gts{i}, tas{i}.traps, 'EdgeColor', 'k')
end
view(vx, vz+45)
axis equal tight
%light('Position',[-1 0 -1]);lighting phong
colorbar('horiz'); caxis([0 numel(unique(ta.traps))]);
title('Traps for sealing top surfaces')
daspect([2, 0.1, 1])

%% Plot subgrids
% Plot VE regions from top surfaces
figure()
subplot(1,2,1)
plotGrid(G, 'facecolor', 'none')
hold on

cm = jet(numel(Gs_all));
for i=1:numel(Gs_all)
    hold on
    if i == 1
        alpha = 0.2;
    else
        alpha = 1;
    end
    plotGrid(Gs_all{i}, 'facecolor', cm(numel(Gs_all)-i+1,:), 'facealpha', alpha)
end

view(vx, vz)
axis equal tight
title({'Top-surface subgrids under semi-perm layers', '(faded for caprock)'})
daspect([1, 0.1, 1])

% Plot remaining fine-perm regions
subplot(1,2,2)
plotGrid(G, 'facecolor', 'none')
hold on
plotGrid(G, fine_rem)
view(vx, vz)
axis equal tight
title('Remaining fine regions')
daspect([1, 0.1, 1])

%% Compute CO2 trapping
hybrid_reports = makeHybridReports(Gts, Gs_all, cmaps_all, fmaps_all, ... % top surface VE+fine regions under semi-perm layers, including caprock                                    
                                    fine_rem, ... % remaining fine regions not directly under semi-perm layer
                                    Gh, {state0_hybrid, states_hybrid{:}}, ...
                                    mh.rock, mh.fluid, schedule_hybrid, ...
                                    [swr, snr], tas, []);
                                                               
fine_reports = makeFineReports(Gts, Gs_all, cmaps_all, fmaps_all, ... % top surface VE+fine regions under semi-perm layers, including caprock                                                                    
                                    G, {state0, states{:}}, ...
                                    model_fine.rock, model_fine.fluid, schedule, ...
                                    [swr, snr], tas, []);                                
                                
%% Plot trapping inventory
fig10 = figure(10); plot(1); ax10 = get(fig10, 'currentaxes');
plotTrappingDistribution(ax10, hybrid_reports, 'legend_location', 'northwest', 'logScale', true)
title('Trapping inventory for hybrid model')
saveas(fig10, strcat(plot_dir, 'trapping_inventory_hybrid'), 'png');

fig11 = figure(11); plot(1); ax11 = get(fig11, 'currentaxes');
plotTrappingDistribution(ax11, fine_reports, 'legend_location', 'northwest', 'logScale', true)
title('Trapping inventory for full-dimensional model')
saveas(fig11, strcat(plot_dir, 'trapping_inventory_fine'), 'png');

%% Compute difference in CO2 vol between fine and hybrid models
ff = 20;
figure(ff);
%plotOutlinedGrid(G, W, bc, sealingFaces | sealingCells_faces);
p = model_hybrid.G.partition;
discr = model_hybrid.G.cells.discretization;
plotCellData(G, discr(p), 'edgealpha', 0.2)
xt = model_hybrid.G.cells.centroids(:,1);
yt = model_hybrid.G.cells.centroids(:,2);
zt = model_hybrid.G.cells.centroids(:,3);
discr_u = unique(discr);
for i=1:numel(discr_u)
    d = discr_u(i);   
    if d ~= 1 % only show discretization number for VE cells, to distinguish them
        xti = mean(xt(discr == d))-lx/nx;  
        yti = mean(yt(discr == d))-ly/ny;
        zti = mean(zt(discr == d))-lz/nz;    
        %text(xti, 0, zti, string(d), 'FontSize', 18, 'Color', 'black');
        text(xti, 0, zti, string(d), 'FontSize', 18, 'FontName', 'Castellar', 'Color', 'black');
    end
end
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])

sn_f = states{end}.s(:,2); % fine saturations
sn_h = states_hybrid{end}.s(:,2); % hybrid saturations
ff = ff + 1;
% VOLUME MISMATCH FOR COARSE REGIONS
[diff_coarse, var_coarse] = Volumes.diffCO2SatCoarse(model_hybrid, model_fine, sn_h, sn_f, ff, plot_dir);

sn_hf = states_hybrid_fs{end}.s(:,2);
ff = ff + 1;
% VOLUME MISMATCH FOR FINE CELLS IN DIFFERENT DISCRETIZATION REGIONS
[diff_fine, var_fine] = Volumes.diffCO2SatFine(model_hybrid, sn_hf, sn_f, ff, plot_dir);

% COMPARE EXITED VOLUMES
ff = ff + 1;
[Vh_exit, Vf_exit] = Volumes.plotExitedVolumes(model_fine, states_hybrid_fs, states, ...
                                            schedule, ff, plot_dir);

% Compare volumes in discretization regions connected to semi-perm layers
ff = ff + 1;
[vols_hybrid, vols_fine] = Volumes.plotSealingLayersVols(model_hybrid, model_fine, ...
                                                        states_hybrid_fs, states, ff);

%% Plot CO2 saturation for each model
fafa = find(sealingFaces | sealingCells_faces);

c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

for i = 1:20:numel(states)
    f1 = figure(1); clf
    plotCellData(model.G, states{i}.s(:,2), 'edgec', 'none');
    plotFaces(G, (1:G.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(G, fafa, 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
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
    plotFaces(G, (1:G.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(G, fafa, 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
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
   section_2d = (jj == 5);
   
   cells_2d = model.G.cells.indexMap(section_2d);
   n_2d = unique(model.G.faces.neighbors(fafa, :));
   cells_2d_faces = cells_2d(ismember(cells_2d, n_2d));     
   faces_2d = gridCellFaces(model.G, cells_2d_faces);
   faces_2d = fafa(ismember(fafa, faces_2d));
   
   for i=1:20:numel(states)%fix(2*numel(states)/3)
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
       
       saveas(f6, sprintf(strcat(plot_dir, 'jj5_fine_sat_%d'), i), 'png');
   end
   
   for i=1:20:numel(states_hybrid)%fix(2*numel(states_hybrid)/3)
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
       
       saveas(f7, sprintf(strcat(plot_dir, 'jj5_hybrid_sat_%d'), i), 'png');
   end
end



%% Functions:
function [v_hybrid, v_fine] = discrCO2Vol(d, model_hybrid,  model_fine, state_hybrid_fs, state_fine)
    discr = model_hybrid.G.cells.discretization == d;
    p = model_hybrid.G.partition;
    discr = discr(p);
    vol = model_fine.G.cells.volumes(discr);
    snh = state_hybrid_fs.s(:,2);    
    snh = snh(discr);
    snf = state_fine.s(:,2);
    snf = snf(discr);
    pv = poreVolume(model_fine.G, model_fine.rock);
    pv = pv(discr);
    
    v_hybrid = pv.*snh;
    v_fine = pv.*snf;
end
