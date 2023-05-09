%% Hybrid model for domains with sealing layers
% Thin shale layers either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layers are refined to fine-scale cells.

% 2D:
% - horizontal, rectilinear grid
% - long sloped, curvilinear grid
% 3D:
% - standard grid

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
useAdaptive = true; % overrides useFaceConstraint in setupSlopedGrid
run3D = false;
sloped = true;

trans_mult = 1e-5; % 1e-4
trans_mult = ~useAdaptive*(useFaceConstraint*trans_mult + ~useFaceConstraint) ...
                + useAdaptive*trans_mult;

if run3D
    nx = 80; ny = 8; nz = 40; % 100, 8, 50
	lx = 500; ly = 25; lz = 250;
    [state0, models, schedule, ...
    isFineCells, sealingFaces] = setupSlopedGrid3D(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    vx = 30; vz = 20; % view angles
    hybrid_folder = 'hybrid3D';
else
    nx = 75; ny = 1; nz = 50; % nz = 50
    lx = 750; ly = 1; lz = 250;
    if sloped
        [state0, models, schedule, ...
         isFineCells, sealingFaces] = setupSlopedGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
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
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/book_illustrations/%s/test/%s/adaptive_sealing_nz%d_quadratic/'), hybrid_folder, geometry_folder, nz);   
elseif useFaceConstraint
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/face_lowperm_nz%d/'), hybrid_folder, geometry_folder, nz);       
else    
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/cell_lowperm_nz%d/'), hybrid_folder, geometry_folder, nz);    
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
bcCells = isFineCells.bc;
fineCells = sealingCells | wellCells | bcCells;

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
% problem = packSimulationProblem(state0, model_fine, schedule, ...
%     'horz_finescale', 'NonLinearSolver', nls);

% Simulate and get the output
%simulatePackedProblem(problem);
%[ws, states, report] = getPackedSimulatorOutput(problem);
[ws, states] = simulateScheduleAD(state0, model_fine, schedule, 'NonLinearSolver', nls);

%% Simple VE.
% Simulate standard VE model, not accounting for sealing faces.
model_ve = convertToHybridModel(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'horz_ve', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_ve);
[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);

%[ws_ve, states_ve] = simulateScheduleAD(state0_ve, model_ve, schedule_ve, 'NonLinearSolver', nls);

states_ve_fs = convertHybridStates(model_ve, model_fine, states_ve, 'schedule', schedule, 'convert_flux', true); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
[model_hybrid, model_coarse] = convertToHybridModel(model, fineCells, ...
                                            'sealingFaces', find(sealingFaces), ... % same as find(model.operators.T_all == 0)
                                            'sealingCells', sealingCells, ...
                                            'sealingCells_faces', sealingCells_faces, ...
                                            'setSubColumnsFine', false, ...
                                            'multiplier', trans_mult, ...
                                            'sumTrans', true, ...
                                            'pe_rest', pe_rest);
                                       
model_hybrid.fluid.pcWG = @(s, n_sealing) ...
                            Capillary.runStandardPcSharp(s, dummy_s, swr, snr, pe_sealing, ...
                                                pe_rest, n_sealing, model_hybrid.G);                                                                                                      
                                        
schedule_hybrid = upscaleSchedule(model_hybrid, schedule);
state0_hybrid = upscaleStateHybrid(model_hybrid, model, state0);

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
%states_hybrid_fs = convertMultiVEStates_res(model_hybrid, model_fine, states_hybrid, states);
states_hybrid_fs = convertHybridStates(model_hybrid, model_fine, states_hybrid, 'schedule', schedule, 'convert_flux', true); % send in fine schedule to extract fine boundary faces

mh = model_hybrid;
oph = mh.operators;
n = oph.N;
shf = states_hybrid_fs;

%% Compute difference in CO2 vol between fine and hybrid models

% Plot discretization regions
ff = 20;
figure(ff);
%plotOutlinedGrid(G, W, bc, sealingFaces | sealingCells_faces);
p = model_hybrid.G.partition;
discr = model_hybrid.G.cells.discretization;
plotCellData(G, discr(p), 'edgealpha', 0.2)
xt = model_hybrid.G.cells.centroids(:,1);
zt = model_hybrid.G.cells.centroids(:,3);
discr_u = unique(discr);
for i=1:numel(discr_u)
    d = discr_u(i);   
    if d ~= 1 % only show discretization number for VE cells, to distinguish them
        xti = mean(xt(discr == d))-lx/nx;   
        zti = mean(zt(discr == d))-lz/nz;    
        %text(xti, 0, zti, string(d), 'FontSize', 18, 'Color', 'black');
        text(xti, 0, zti, string(d), 'FontSize', 18, 'FontName', 'Castellar', 'Color', 'black');
    end
end
view(vx, vz)
axis equal tight
daspect([1, 0.1, 1])

% Compare volumes in discretization regions connected to semi-perm layers
[vols_hybrid, vols_fine] = Volumes.plotSealingLayersVols(model_hybrid, model_fine, ...
                                                        states_hybrid_fs, states, ff);

sn_f = states{end}.s(:,2); % fine saturations
sn_h = states_hybrid{end}.s(:,2); % hybrid saturations

[diff_coarse, var_coarse] = Volumes.diffCO2SatCoarse(model_hybrid, model_fine, sn_h, sn_f, 10, plot_dir);

% VOLUME MISMATCH FOR FINE CELLS IN DIFFERENT DISCRETIZATION REGIONS
sn_hf = states_hybrid_fs{end}.s(:,2);

[diff_fine, var_fine] = Volumes.diffCO2SatFine(model_hybrid, sn_hf, sn_f, 11, plot_dir);

% COMPARE EXITED VOLUMES
[Vh_exit, Vf_exit] = Volumes.plotExitedVolumes(states_hybrid_fs, states, ...
                                            schedule, 12, plot_dir);

%% Plot CO2 saturation for each model
fafa = find(sealingFaces | sealingCells_faces);

% c1 = [48, 37, 255]/255;
% c2 = [0, 255, 0]/255;
% cc = interp1([0; 1], [c1; c2], (0:0.01:1)');
c1 = [255, 255, 255]/255;
c2 = [48, 37, 255]/255;
c3 = [0, 255, 0]/255;
cc = interp1([0; 1e-6; 1], [c1; c2; c3], (0:1e-6:1)');

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
   section_2d = (jj == 5);
   
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
cells_top = G.cells.indexMap(kk == min(kk));
cells_bottom = G.cells.indexMap(kk == max(kk)); 
faces_top = [];
for c=cells_top'
    f_top = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1, :);
    f_top = f_top(f_top(:,2) == 5, 1);
    faces_top = cat(1, faces_top, f_top);
end

sort_f = sort(unique(G.cells.faces(:,1)));
%faces_top = sort_f(1:numel(cells_top));
faces_bottom = sort_f(end-numel(cells_bottom)+1:end);

nstep = numel(schedule.step.val);
end_inj = find(schedule.step.control == 1, 1, 'last');
end_mig = nstep;

substeps = [end_inj, end_mig];
%substeps = [nstep/10, nstep/5];
names = {'injected', 'migrated'};
model_names = {'fine', 'hybrid', 'coarse', 've'};

xl_min = -max(G.cells.centroids(:,1))/50;
xl_max = max(G.cells.centroids(:,1))*(1+1/50);
zl_min = -max(G.cells.centroids(:,3))/50;
zl_max = max(G.cells.centroids(:,3))*(1+1/50);

fign = 13;
for i = 1:numel(substeps)
    ss = substeps(i);    
    for j = 1:2
        ff = UtilFunctions.fullsizeFig(fign);
        fign = fign + 1;
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
        %subplot(1, 2, 2-mod(j,2));
        %plotOutlinedGrid(g, W, bc, bcfaces | sealingCells_faces, 'lw', 0.8); 
        plotCellData(g, s, 'EdgeColor', 'none')
        if ~isfield(g, 'parent')                   
            plotFaces(g, fafa, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.8)                       
            hold on
            plotFaces(g, faces_top, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.1) 
            hold on
            plotFaces(g, faces_bottom, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.1) 
            %plotFaces(g, boundaryFaces(g), 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.7)
        else
            plotOutlinedGrid(g.parent, W, bc, bcfaces | sealingCells_faces, 'lw', 0.8);
        end

        view(vx, vz)
        axis equal tight
        if run3D
            daspect([1, 0.1, 1]) 
        else
            %daspect([1, 0.1, 0.5])
        end
        %xlim([xl_min, xl_max]);
        %zlim([zl_min, zl_max]);
        colormap(cc); caxis([0,1-swr]); colorbar('off')
        title(nm);
        
        if ~mod(j,1) % 2
            set(gcf, 'Name', names{i});            
            saveas(ff, strcat(plot_dir, names{i}, '_', model_names{j}));
            saveas(ff, strcat(plot_dir, names{i}, '_', model_names{j}), 'png');
            saveas(ff, strcat(plot_dir, names{i}, '_', model_names{j}), 'pdf');
            %ff = UtilFunctions.fullsizeFig(fign);
            %fign = fign + 1;
        end
    end
    %clf;
end
