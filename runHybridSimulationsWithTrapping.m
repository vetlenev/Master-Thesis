%% Hybrid model for diffuse leakage through thin horizontal shale
% Thin shale layer is either represented as face constraint or cell
% constraint. For face constraint, transmissibility is set to zero. For
% cell constraint, thin layer is refined to fine-scale cells.

gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui test-suite
    
rootdir = strrep(ROOTDIR, '\', '/');
%data_dir = strcat(rootdir, '../Master-Thesis/book_illustrations/hybrid3D/test/data');
data_dir = strcat(rootdir, '../Master-Thesis/hybrid2D/caseStochastic/data');
mkdir(data_dir);

my_seed = 662;
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%pe_sealing = 7*10^4*Pascal; %7*10^4*Pascal; % higher entry pressure in aquitards
pe_sealing = [2*10^5, 7*10^4, 3*10^4]*Pascal;
pe_rest = 5*10^3*Pascal;%5*10^3*Pascal;
% -----------------
pe_regions = [pe_sealing, pe_rest];
% -----------------
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
standard = true;

n_rel = 1.5;
n_layers = 10;
output_filename = 'stochastic';

if sloped && standard && strcmp(output_filename, 'stochastic')
    geometry_folder = 'caseStochastic';
elseif sloped && standard
    geometry_folder = 'caseMultilayered';
elseif sloped && ~standard
    geometry_folder = 'caseSlopedLong';
else
    geometry_folder = 'caseSimple';
end

trans_mult = 1e-6; % 1e-5
trans_mult = ~useAdaptive*(useFaceConstraint*trans_mult + ~useFaceConstraint) ...
                + useAdaptive*trans_mult;

if run3D    
    output_filename = strcat(output_filename, '_3D');
    if standard
        nx = 60; ny = 10; nz = 40; % 100, 8, 50
	    lx = 500; ly = 50; lz = 250; % 500, 25, 250
        [state0, models, schedule, ...
        isFineCells, sealingFaces, allSealingFaces] = setupSlopedGrid3D(useFaceConstraint, useAdaptive, ...
                                                        [nx,ny,nz], [lx,ly,lz], trans_mult);
    else
        nx = 16; ny = 16; nz = 70; % 100, 8, 50
	    lx = 750; ly = 500; lz = 400; % 500, 25, 250
        [state0, models, schedule, ...
        isFineCells, sealingFaces, allSealingFaces] = setupSlopedTallGrid3D(useFaceConstraint, useAdaptive, trans_mult, ...
                                                                                [nx,ny,nz], [lx,ly,lz]);
    end                                             
    vx = 30; vz = 20; % view angles
    hybrid_folder = 'hybrid3D';
else   
    output_filename = strcat(output_filename, '_2D');    
    nx = 50; ny = 1; nz = 100; % nx = 50, nz = 150
    lx = 300; ly = 10; lz = 1000;     

    if sloped
        [state0, models, schedule, ...
         isFineCells, sealingFaces, allSealingFaces] = setupStochasticGrid(useFaceConstraint, useAdaptive, ...
                                                                    [nx,ny,nz], [lx,ly,lz], trans_mult, standard, n_layers, n_rel);
    else
        [state0, models, schedule, ...
         isFineCells, sealingFaces] = setupHorizontalGrid(useFaceConstraint, useAdaptive, ...
                                                    [nx,ny,nz], [lx,ly,lz], trans_mult);
    end
    vx = 0; vz = 0;
    hybrid_folder = 'hybrid2D';
end

output_filename = strcat(output_filename, '_layers', string(n_layers), '_rel', strrep(string(n_rel),'.','p'));
output_filename = char(output_filename);

if useAdaptive    
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/adaptive_%s/'), hybrid_folder, geometry_folder, output_filename);   
elseif useFaceConstraint
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/face_%s/'), hybrid_folder, geometry_folder, output_filename);       
else    
    plot_dir = sprintf(strcat(rootdir, '../Master-Thesis/%s/%s/figs/cell_%s/'), hybrid_folder, geometry_folder, output_filename);   
end
mkdir(plot_dir);

% Split into two models: only fine model has transmissibility multiplier
% applied - for VE models it depends if we're using face or cell
% constraints
model = models.original; % used as input to VE models
model_fine = models.fine;
G = model_fine.G;
rock = model_fine.rock;

[ii, jj, kk] = gridLogicalIndices(G);

sealingCells = isFineCells.sealingCells;
sealingCells = any(cell2mat(sealingCells), 2);
% -----
sealingCells_faces = false(size(sealingFaces));
sealingCells_facesAll = vertcat(isFineCells.sealingCells_faces{:});
sealingCells_faces(sealingCells_facesAll) = true;
% -----
wellCells = isFineCells.well;
openBCCells = isFineCells.bc;

remFineCells = wellCells | openBCCells;
fineCells = sealingCells | remFineCells;

wellFaces = false(size(sealingFaces));
wellFacesAll = addConfiningLayers(G,'type','cells','full_dim',run3D,'cells',wellCells);
wellFaces(wellFacesAll) = true;

swr = model_fine.fluid.krPts.w(1); % reisdual water sat
snr = model_fine.fluid.krPts.g(1); % residual CO2 sat
dummy_s = linspace(0, 1, model_fine.G.cells.num)';

sealingCellsFine = model_fine.G.cells.indexMap(sealingCells);
% model_fine.fluid.pcWG = @(s, varargin) Capillary.runStandardPcSharp(s, dummy_s, swr, snr, pe_sealing, pe_rest, ...
%                                            sealingCellsFine, model_fine.G);
model_fine.fluid.pcWG = @(s, varargin) Capillary.runHybridPcSharp(s, dummy_s, swr, snr, pe_regions, [], model_fine);

W = schedule.control(1).W;
bc = schedule.control(1).bc;
T_all = model.operators.T_all;
T = model.operators.T;
nn = G.faces.neighbors;

perm_sealing = uniquetol(rock.perm(sealingCellsFine))'/(milli*darcy); % NB: assumes homogeneous sealing layers
perm_all = uniquetol(rock.perm)'/(milli*darcy);
perm_rest = uniquetol(setdiff(perm_all, perm_sealing));

pe_all = zeros(G.cells.num, 1);
for i=1:numel(perm_all)
    reg_i = rock.perm/(milli*darcy) == perm_all(i);
    pe_all(reg_i) = pe_regions(i);
end


%% Plot fine grid
f30 = figure(30);
%plotOutlinedGrid(G, W, bc, model.operators.T_all);
plotWell(G, W, 'color', 'r')
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')
%plotFaceData(G, model.operators.T_all);

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting

plotCellData(G, double(fineCells), fineCellsIdx, 'facecolor', 'green'); % requires numeric input, not logical
view(vx, vz)
axis equal tight
setDaspect(run3D, standard);
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title({'Fine scale grid', 'Imposed fine cells highlighted in green'})
%saveas(f30, strcat(plot_dir, '/fine_regions'), 'png');

%% More plots
f31 = figure(31);
plotOutlinedGrid(G, W, bc, sealingFaces);
plotOutlinedGrid(G, W, bc, sealingCells_faces);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')

fineCellsIdx = G.cells.indexMap(fineCells); % for plotting
log10perm = log10(convertTo(rock.perm, milli*darcy));

plotCellData(G, log10perm); % requires numeric input, not logical
if ~useFaceConstraint
    caxis([min(log10perm)-0.1*min(log10perm), ...
            max(log10perm)+0.1*max(log10perm)])
end
view(vx, vz)
axis equal tight
setDaspect(run3D, standard);
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title('Permeability field')
colormap(flipud(jet))
colorbar('east')
%saveas(f31, strcat(plot_dir, '/perm_field'), 'png');

% Plot van Genuchten capillary pressure
swMin = 0.2; swr_pc = 0.1; snr_pc = 0.15;
alpha = -1/(10^4*Pascal); n = 1.5; gamma_D = 2; gamma_I = 1.5;
pc_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.pc_func(Sw, sne, swr, alpha, n, gamma);
sw_func = @(Sw, sne, swr, alpha, n, gamma) Hysteresis.sw_func(Sw, sne, swr, alpha, n, gamma);
pc_scan = Hysteresis.Genuchten(dummy_s, swMin, snr_pc, swr_pc, pc_func, sw_func, alpha, n, gamma_D, gamma_I);

% fig3 = figure(3);
% plot(1-dummy_s, Capillary.PcGas(dummy_s, swr, snr, pe_rest, 4), 'LineWidth', 1.5);
% %plot(1-dummy_s, pc_scan, 'LineWidth', 1.5)
% xlabel('Water saturation');
% ylabel('[Pa]')
% xlim([swr, 1]);
% title('Capillary pressure function');
% saveas(fig3, strcat(plot_dir, '/cap_pres'), 'pdf');
% hold off

%% Nonlinear solver
nls = NonLinearSolver('maxIterations', 70);

%% Pack and simulate problem
problem = packSimulationProblem(state0, model_fine, schedule, ...
            strcat('finescale/',geometry_folder), 'NonLinearSolver', nls, 'Name', output_filename);

% Simulate and get the output
simulatePackedProblem(problem);
[ws, states, report] = getPackedSimulatorOutput(problem);
%[ws, states] = simulateScheduleAD(state0, model_fine, schedule, 'NonLinearSolver', nls);

%% Simple VE.
% Simulate standard VE model, not accounting for sealing faces.
model_ve = convertToMultiVEModel_test(model);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleStateHybrid(model_ve, model, state0);

problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
            strcat('ve/',geometry_folder), 'NonLinearSolver', nls, 'Name', output_filename);

[ok, status] = simulatePackedProblem(problem_ve);
[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);

%[ws_ve, states_ve] = simulateScheduleAD(state0_ve, model_ve, schedule_ve, 'NonLinearSolver', nls);

states_ve_fs = convertMultiVEStates_test(model_ve, model_fine, states_ve, 'schedule', schedule, 'convert_flux', true); % retrieve fine-scale states

%% Setup Hybrid model
% Simulate hybrid VE model, accounting for diffuse leakage at sealing face
%pe_all_hybrid = accumarray(ph, pe_all, [], @mode);
[model_hybrid, model_coarse] = convertToMultiVEModel_test(model, fineCells, 'remFineFaces', find(wellFaces), ...                                                    
                                            'sealingFaces', find(sealingFaces), ... % same as find(model.operators.T_all == 0)
                                            'sealingCells', sealingCells, ...
                                            'sealingCells_faces', sealingCells_faces, ...
                                            'setSubColumnsFine', false, ...
                                            'multiplier', trans_mult, ...
                                            'sumTrans', true, ...
                                            'pe_all', pe_all); % ! pe_rest !
                                       
% model_hybrid.fluid.pcWG = @(s, n_sealing) ...
%                             Capillary.runHybridPcSharp(s, dummy_s, swr, snr, pe_sealing, ...
%                                                 pe_rest, n_sealing, model_hybrid.G);        
model_hybrid.fluid.pcWG = @(s, subcells) Capillary.runHybridPcSharp(s, dummy_s, swr, snr, ...
                                                                    pe_regions, subcells, model_hybrid);
                                        
schedule_hybrid = upscaleSchedule(model_hybrid, schedule);
state0_hybrid = upscaleStateHybrid(model_hybrid, model, state0);

Gh = model_hybrid.G;
disc = Gh.cells.discretization;
ph = Gh.partition;

% Plot VE partitions on grid
% Here we can see how VE cells exist between sealing layers and can be
% stacked on top of each other.
f4 = figure(4);
plotOutlinedGrid(G, W, bc, sealingFaces | sealingCells_faces | wellFaces);
plotGrid(G, find(disc(ph) == 1), 'edgecolor', 'yellow')
%plotCellData(G, model_hybrid.G.partition,'edgealpha', 0.2)
plotCellData(G, disc(ph),'edgealpha', 0)
colormap("colorcube")
view(vx, vz)
axis equal tight
setDaspect(run3D, standard);
xlabel('Lateral position [m]');
zlabel('Depth [m]');
title('Partition of hybrid VE model')
%saveas(f4, strcat(plot_dir, 'hybrid_discretization'), 'png')

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
all_coarse_cells(vtc_A) = 1;
all_coarse_cells(vtc_B) = 1;
%all_coarse_cells(1:200) = 1;
plotCellData(model_hybrid.G, all_coarse_cells)
%plotCellData(model_hybrid.G, model_hybrid.G.cells.discretization);
view(vx, vz)
axis equal tight
setDaspect(run3D, standard);
title('Combined veToFineVertical and veHorizontal transition regions')

%% Set capillary state function
nls.useRelaxation = true;

model_hybrid = model_hybrid.validateModel;
model_hybrid.FlowPropertyFunctions.CapillaryPressure = HybridCapillaryPressure(model_hybrid);

%% Simulate Hybrid model
problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
            strcat('hybrid/',geometry_folder), 'NonLinearSolver', nls, 'Name', output_filename);

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

top_cells = G.cells.indexMap(kk == min(kk)); % global top surface
bf = boundaryFaces(G);
all_top_faces = zeros(numel(top_cells), 6);
for i=1:numel(top_cells)
    all_top_faces(i,:) = G.cells.faces(G.cells.facePos(top_cells(i)):G.cells.facePos(top_cells(i)+1)-1, 1);
end
top_faces = all_top_faces(:,5); % faces with normal vector pointing upwards
bf_bc = intersect(bf, bc.face); % boundary faces with open bc
top_faces_closed = top_faces(~ismember(top_faces, bf_bc)); % only extract global top surface faces NOT given as open boundary

if ~isempty(top_faces_closed)
    [Gs_glob, cmap_glob, fmap_glob, ~] = ExtractLayerSubgrid(G, Gh, top_faces_closed, ...
                                                        isFineCells.sealingCells, allSealingFaces); % cell constraints, face constraints
else
    Gs_glob = [];
    cmap_glob = [];
    fmap_glob = [];
end
% OR: Do trap analysis for combined top surface
if ~isempty(Gs_glob)
    Gt = topSurfaceGrid(Gs_glob);
    ta = trapAnalysis(Gt, true);
    Gts = cat(1, Gts, Gt);
    tas = cat(1, tas, ta);
end

%% Subgrids and traps for VE top surfaces

sealingLayers = [isFineCells.sealingBottom, allSealingFaces]; % merge all faces defining top surfaces

for i=1:numel(sealingLayers)
    subfaces = sealingLayers{i}; % bottom faces of sealing layer    
    [Gs, cmap, fmap, nmap] = ExtractLayerSubgrid(G, Gh, subfaces, ...
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

Gs_all = [Gs_glob; Gsi]; % [{Gs_glob}; Gsi];
cmaps_all = [cmap_glob; cmaps];
fmaps_all = [fmap_glob; fmaps];

%% Plot trapping
% Remaining fine cells not subject to trapping
top_surface_cells = vertcat(cmaps_all{:}); % should be unique
all_rem = setdiff(G.cells.indexMap, top_surface_cells);
all_rem_hybrid = Gh.partition(all_rem);
isVE_cells = Gh.cells.discretization(all_rem_hybrid);
ve_rem = all_rem(isVE_cells > 1);
fine_rem = all_rem(isVE_cells == 1);

figure(15)
%plotCellData(Gt, ones(Gt.cells.num,1), 'EdgeColor', 'none');
max_traps = 0;
for i=1:numel(Gts)
    plotCellData(Gts{i}, tas{i}.traps, 'EdgeColor', 'k')
    max_traps = max(max_traps, numel(unique(tas{i}.traps)));
end
view(vx, vz+45)
axis equal tight
%light('Position',[-1 0 -1]);lighting phong
colorbar('horiz'); clim([0 max_traps]);
title('Traps for sealing top surfaces')
setDaspect(run3D, standard);

%% Plot subgrids
% Plot VE regions from top surfaces
fig_topsurf = figure();
subplot(1,2,1)
plotGrid(G, 'facecolor', 'none')
hold on

cm = jet(numel(Gs_all));
for i=1:numel(Gs_all)
    hold on
    if i == 1 && ~isempty(Gs_glob)
        alpha = 0.2;
    else
        alpha = 1;
    end
    plotGrid(Gs_all{i}, 'facecolor', cm(numel(Gs_all)-i+1,:), 'facealpha', alpha)
end

view(vx, vz)
axis equal tight
title({'Top-surface subgrids under semi-perm layers', '(faded for caprock)'})
setDaspect(run3D, standard);

% Plot remaining fine-perm regions
subplot(1,2,2)
plotGrid(G, 'facecolor', 'none')
hold on
plotGrid(G, fine_rem, 'facecolor', 'green')
hold on
plotGrid(G, ve_rem, 'facecolor', 'red')
view(vx, vz)
axis equal tight
title({'Remaining regions.', 'Green: fine cells.', 'Red: ve cells.'})
setDaspect(run3D, standard);

%saveas(fig_topsurf, strcat(plot_dir, 'top_surface_subgrids'), 'png');

%% Compute CO2 trapping
hybrid_reports = makeHybridReports(Gts, Gs_all, cmaps_all, fmaps_all, ... % top surface VE+fine regions under semi-perm layers, including caprock                                    
                                    fine_rem, ve_rem, ... % remaining fine regions not directly under semi-perm layer
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
%saveas(fig10, strcat(plot_dir, 'trapping_inventory_hybrid'), 'pdf');

fig11 = figure(11); plot(1); ax11 = get(fig11, 'currentaxes');
plotTrappingDistribution(ax11, fine_reports, 'legend_location', 'northwest', 'logScale', true)
title('Trapping inventory for full-dimensional model')
%saveas(fig11, strcat(plot_dir, 'trapping_inventory_fine'), 'pdf');

%% Differences in trapping inventory
ff = 20;
fig_ff = figure(ff); plot(1); ax_ff = get(fig_ff, 'currentaxes');
[max_diff, max_diff_year, max_CO2, ...
    sum_diff, sum_CO2] = plotTrappingDifference(ax_ff, hybrid_reports, fine_reports, 'legend_location', 'northwest', 'logScale', true);
max_trap_diff = round(max_diff.tot/max_CO2.tot*100, 2);
title({'Normalized absolute difference in trapping inventory', ...
        sprintf('Max difference: %.2d MT, at year: %.1f,', max_diff.tot, max_diff_year.tot), ...
        sprintf('which equals %.1f %% of injected CO2.', max_trap_diff)})
%saveas(fig_ff, strcat(plot_dir, 'trapping_diff'), 'pdf');

%% Compute difference in CO2 vol between fine and hybrid models
control1 = schedule.step.control == 1;
inj_rate = schedule.control(1).W.val;
tot_inj = inj_rate .* cumsum(schedule.step.val) .* control1;
tot_inj = max(tot_inj);

ff = 22;
[vols_hybrid, vols_fine, discr_region, ff] = Volumes.plotSealingLayersVols(model_hybrid, model_fine, ...
                                                                        states_hybrid_fs, states, ff, false);

ff = ff + 1;
figure(ff);
plotFaces(model_fine.G, sealingCells_faces)
p = model_hybrid.G.partition;
discr = model_hybrid.G.cells.discretization;
plotCellData(G, discr(p), 'edgealpha', 0.2)
xt = model_hybrid.G.cells.centroids(:,1);
yt = model_hybrid.G.cells.centroids(:,2);
zt = model_hybrid.G.cells.centroids(:,3);
discr_u = unique(discr);
discr_used = [];
for i=1:numel(discr_u)
    d = discr_u(i);   
    rve_nonzero_flux = cellfun(@isempty, regexp(discr_region(:), string(d), 'match'));
    if d ~= 1 && ~all(rve_nonzero_flux) % only show discretization number for RVE cells with nonzero flux
        xti = mean(xt(discr == d))-lx/nx;  
        yti = mean(yt(discr == d))-ly/ny;
        zti = mean(zt(discr == d))-lz/nz;    
        %text(xti, 0, zti, string(d), 'FontSize', 18, 'Color', 'black');
        text(xti, 0, zti, string(d), 'FontSize', 18, 'FontName', 'Castellar', 'Color', 'black');
        discr_used = cat(1, discr_used, d);
    end
end
discr_plot = discr;
discr_plot(~ismember(discr_plot, discr_used)) = 0;
plotCellData(G, discr_plot(p), 'edgealpha', 0.2)
view(vx, vz)
axis equal tight
xlabel('Lateral position [m]')
zlabel('Depth [m]')
c1 = [255, 255, 255]/255;
c4 = [0, 255, 0]/255;
c3 = [255, 255, 0]/255;
c2 = [255, 165, 0]/255;
discr_plot_u = unique(discr_plot);
cmap_discr = interp1([0; linspace(discr_plot_u(2), discr_plot_u(end), 7)'], ...
                         [c1; c2; c3; c4; c3; c2; c3; c4], ...
                         (0:1:discr_plot_u(end))');
%ccube = flipud(colorcube(50));
colormap(cmap_discr)
setDaspect(run3D, standard);
%saveas(ff, strcat(plot_dir, 'partition_numbering'), 'png');

sn_f = states{end}.s(:,2); % fine saturations
sn_h = states_hybrid{end}.s(:,2); % hybrid saturations
sn_hf = states_hybrid_fs{end}.s(:,2);


% COMPARE EXITED VOLUMES
ff = ff + 1;
[Vh_exit, Vf_exit] = Volumes.plotExitedVolumes(model_fine, states_hybrid_fs, states, ...
                                            schedule, ff, plot_dir);

% Compare volumes in discretization regions connected to semi-perm layers
ff = ff + 1;
[vols_hybrid_CO2, vols_fine_CO2, ff] = Volumes.plotSealingLayersCO2(model_hybrid, model_fine, ...
                                                            states_hybrid_fs, states, schedule, plot_dir, ff);

%% Volumn difference - VE relaxed columns
ff = ff + 1;
vols_RVE_diff = abs(vols_hybrid - vols_fine);
t_cumsum = cumsum(schedule.step.val)./year();

fig_ff = figure(ff); plot(1); ax_ff = get(fig_ff, 'currentaxes');
[max_diff_net, year_net, max_CO2_net] = plotRVEDifference(ax_ff,vols_hybrid,vols_fine,t_cumsum,discr_region,'logScale',true);
max_RVE_diff = round(max_diff_net/max_CO2_net*100, 2);
title({'Normalized absolute difference in CO2 volume for RVE columns.', ...
        sprintf('Max difference: %.2d m^3, at year: %.1f,', max_diff_net, year_net), ...
        sprintf('which equals %.1f %% of CO2 in RVE columns.', max_RVE_diff)})
%saveas(fig_ff, strcat(plot_dir, 'vols_rve_diff'), 'pdf');

%% Store settings
tot_time = sum(schedule.step.val)/year();
inj_stop = round((sum(schedule.step.val.*(schedule.step.control == 1))/year())/tot_time, 1);
pv_rate = inj_rate/sum(poreVolume(G,rock))*(inj_stop*tot_time*year());
n_sealing = numel(isFineCells.sealingCells);

dirpath_fine = strcat(mrstOutputDirectory,'\',problem.BaseName,'\simulation_settings.txt');
dirpath_ve = strcat(mrstOutputDirectory,'\',problem_ve.BaseName,'\simulation_settings.txt');
dirpath_hybrid = strcat(mrstOutputDirectory,'\',problem_hybrid.BaseName,'\simulation_settings.txt');

time_fine = zeros(numel(report), 1);
time_ve = zeros('like', time_fine);
time_hybrid = zeros('like', time_fine);
for i=1:numel(report)
    time_fine(i) = report{i}.WallTime;
    time_ve(i) = report_ve{i}.WallTime;
    time_hybrid(i) = report_hybrid{i}.WallTime;
end
time_fine_tot = round(sum(time_fine),1);
time_fine_median = round(median(time_fine),3);
time_ve_tot = round(sum(time_ve),1);
time_ve_median = round(median(time_ve),3);
time_hybrid_tot = round(sum(time_hybrid),1);
time_hybrid_median = round(median(time_hybrid),3);

data = {pe_sealing,pe_rest,perm_sealing,perm_rest,trans_mult,...
        swr,snr,mode(rock.poro),tot_time,inj_stop,pv_rate,nx,ny,nz, ...
        n_sealing, n_rel, my_seed, max_trap_diff, max_RVE_diff};
other_info = 'Multiple lowperm vals';

UtilFunctions.storeProblemSettings(data, time_fine_tot, time_fine_median, ...
                                    dirpath_fine, problem.Name, other_info);
% UtilFunctions.storeProblemSettings(data, time_ve_tot, time_ve_median, ...
%                                     dirpath_ve, problem_ve.Name, other_info);
UtilFunctions.storeProblemSettings(data, time_hybrid_tot, time_hybrid_median, ...
                                    dirpath_hybrid, problem_hybrid.Name, other_info);

%% Store trapping data
dirpath_trap = strcat(mrstOutputDirectory,'\trapping\layers',string(n_layers));
mkdir(dirpath_trap);

fullfile_trap = strcat(dirpath_trap, '\trap_diff_max.txt');
other_info = 'max difference for each category';
UtilFunctions.storeTrappingData(max_diff, max_CO2, fullfile_trap, problem.Name, other_info);

fullfile_trap = strcat(dirpath_trap, '\trap_diff_sum.txt');
other_info = 'summed difference over all time steps';
UtilFunctions.storeTrappingData(sum_diff, sum_CO2, fullfile_trap, problem.Name, other_info);

%% Plot CO2 saturation for each model
fafa = find(sealingFaces | sealingCells_faces);

c1 = [255, 255, 255]/255;
c2 = [48, 37, 255]/255;
c3 = [0, 255, 0]/255;
cc = interp1([0; 1e-5; 1], [c1; c2; c3], (0:0.001:1)');
% c1 = [48, 37, 255]/255;
% c2 = [0, 255, 0]/255;
% cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

t = cumsum(schedule.step.val)/year();
% 
for i = 1:20:numel(states)
    f1 = UtilFunctions.fullsizeFig(1); clf
    %f1 = figure(1); clf
    plotCellData(model.G, states{i}.s(:,2), 'edgec', 'none');
    plotFaces(G, (1:G.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(G, fafa, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
    view(vx, vz); colormap(cc); caxis([0, 1-swr]); %colorbar('location','southoutside');
    axis tight off
    title({'Fine-scale saturation', sprintf('year: %.1f',t(i))})   
    setDaspect(run3D, standard);
    %saveas(f1, sprintf(strcat(plot_dir, 'fine_sat_%d'), i), 'png');       
end

for i = 1:20:numel(states_hybrid)
    %f2 = figure(2); clf
    f2 = UtilFunctions.fullsizeFig(2); clf    
    plotCellData(model.G, states_hybrid_fs{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, (1:G.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(G, fafa, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
    view(vx, vz); colormap(cc); caxis([0, 1-swr]); %colorbar('location','southoutside');
    axis tight off
    title({'Hybrid reconstructed saturation', sprintf('year: %.1f',t(i))})   
    %pause(0.5);
    setDaspect(run3D, standard);
    %saveas(f2, sprintf(strcat(plot_dir, 'hybrid_sat_%d'), i), 'png');         
end

%% Compute migration speed (m/year)
[tip_speed_fine, avg_speed_fine, z_well] = plumeLocation(model_fine, states, schedule, snr);
[tip_speed_hybrid, avg_speed_hybrid] = plumeLocation(model_hybrid, states_hybrid_fs, schedule_hybrid, snr);

fig_speed = figure();
tt = [0; t];
plot(tt(tt >= 0.01), tip_speed_fine(tt >= 0.01), '-b')
hold on
plot(tt(tt >= 0.01), tip_speed_hybrid(tt >= 0.01), '-r')
hold on
plot(tt(tt >= 0.01), avg_speed_fine(tt >= 0.01), '--b')
hold on
plot(tt(tt >= 0.01), avg_speed_hybrid(tt >= 0.01), '--r')
hold on
plot(tt(tt >= 0.01), repmat(z_well, nnz(tt >= 0.01), 1), '-k')
set(gca, 'XScale', 'log')
xlim([0.01, max(t)])
set(gca, 'YDir', 'reverse')
legend({'Tip: Fine', 'Tip: Hybrid', 'Average: Fine', 'Average: Hybrid'}, ...
            "Location","northwest");
xlabel('Years since simulation start')
ylabel('Depth (m)')
title('Depth of mobile plume')
%saveas(fig_speed, strcat(plot_dir, 'location_plume'), 'pdf');

%% Mean and var difference in depth
tip_plume_fine = tip_speed_fine(tt >= 0.01);
tip_plume_hybrid = tip_speed_hybrid(tt >= 0.01);
avg_plume_fine = avg_speed_fine(tt >= 0.01);
avg_plume_hybrid = avg_speed_hybrid(tt >= 0.01);
t_store = tt(tt >= 0.01);

plume_depth_data = [t_store, tip_plume_fine, tip_plume_hybrid, ...
                             avg_plume_fine, avg_plume_hybrid];

save(strcat(mrstOutputDirectory, '\plume_depth\', problem.Name, '.mat'), 'plume_depth_data')

%% Plot 2D section
if run3D    
   section_2d = (jj == fix(max(jj)/3));
   
   cells_2d = model.G.cells.indexMap(section_2d);
   n_2d = unique(model.G.faces.neighbors(fafa, :));
   cells_2d_faces = cells_2d(ismember(cells_2d, n_2d));     
   faces_2d = gridCellFaces(model.G, cells_2d_faces);
   faces_2d = fafa(ismember(fafa, faces_2d));
   faces_2d = faces_2d(model.G.faces.normals(faces_2d,3) > eps); % remove horizontal-directed faces
   
   for i=1:20:numel(states)
       f6 = figure(6); clf
       sf = states{i}.s(:,2);      
       sf = sf(section_2d);
       plotCellData(model.G, sf, cells_2d, 'edgec', 'none')
       plotFaces(G, faces_2d, 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
       view(0, 0); colormap(cc); caxis([0,1-swr]);
       %colorbar('location', 'southoutside');
       axis tight on    
       pause(0.2)
       title({'Fine saturation', sprintf('year: %.1f',t(i))}) 
       
       saveas(f6, sprintf(strcat(plot_dir, 'jj5_fine_sat_%d'), i), 'png');
   end
   
   for i=1:20:numel(states_hybrid)%fix(2*numel(states_hybrid)/3)
       f7 = figure(7); clf
       sh = states_hybrid_fs{i}.s(:,2);      
       sh = sh(section_2d);
       plotCellData(model.G, sh, cells_2d, 'edgec', 'none')
       plotFaces(G, faces_2d, 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
       view(0, 0); colormap(cc); caxis([0,1-swr]);
       %colorbar('location', 'southoutside');
       axis tight on
       pause(0.2)
       title({'Hybrid saturation', sprintf('year: %.1f',t(i))}) 
       
       saveas(f7, sprintf(strcat(plot_dir, 'jj5_hybrid_sat_%d'), i), 'png');
   end
end

%% Plot result at injection stop and end of simulation
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

fign = 12;
for i = 1:numel(substeps)
    ss = substeps(i);    
    for j = 1:4
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
            nm = 'Coarse solution';
        else
            g = G;
            s = states_ve_fs{ss}.s(:, 2);
            nm = 'VE reconstructed (no layers)';
        end
        %subplot(1, 2, 2-mod(j,2));       
        plotCellData(g, s, 'EdgeColor', 'none')           
        if ~isfield(g, 'parent')
            plotOutlinedGrid(g, W, bc, sealingCells_faces | sealingFaces, 'lw', 1.5); 
            plotFaces(g, fafa, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.7)
        else
            plotOutlinedGrid(g.parent, W, bc, sealingCells_faces | sealingFaces, 'lw', 1.5); 
        end

        view(vx, vz)
        axis equal tight        
        setDaspect(run3D, standard);         
        %xlim([xl_min, xl_max]);
        %zlim([zl_min, zl_max]);
        colormap(cc); caxis([0,1-swr]); colorbar('off')
        title(nm);
        
        if ~mod(j,1)
            set(gcf, 'Name', names{i});            
            saveas(ff, strcat(plot_dir, names{i}, '_', model_names{j}), 'png');            
            %fign = fign + 1;
        end
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

function setDaspect(run3D, standard, varargin)
    if run3D && ~standard
       daspect([1, 1, 0.5]); % 0.2, 0.25, 0.5
       if nargin > 2
           colorbar(varargin{1});
       end
    elseif standard
       daspect([0.5, 0.1, 1]);
       if nargin > 2
           colorbar(varargin{1});
       end
    else
       daspect([1, 0.1, 1]);
       if nargin > 2
           colorbar(varargin{1});
       end
    end
end
