%% Effect of number of lowperm layers
% Illustrate how different number of lowperm layers affect CO2 migration.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 8211;

%% Define 2D grid
nx = 60; ny = 1; nz = 40; % 100 50
lx = 800; ly = 1; lz = 400; % 1000 350
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);
x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);
dz = mean(diff(z)); % average cell spacing z-dir
dx = mode(diff(x));

%% Define rock and fluid objects
lowperm = 10*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [G.cells.num 1]);
perms = {lowperm, baseperm};

poro = 0.3;

%% Directories
n_layers = 30;
n_lowperm_layers = n_layers; % 13, 26, 40
n_imperm_layers = 0;

plot_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/plotsNumLayers');
data_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/dataNumLayers');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_layers);
data_dir = sprintf(strcat(data_base_dir, '/layers_%d'), n_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Set seed
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Set up grid with layers
grid = GridSetup(G, n_lowperm_layers, n_imperm_layers, perms);

[line_idx, anticline_idx] = grid.DefineLayers();

[added_layers, trapped_cells] = grid.GenerateLayers(line_idx, anticline_idx);

all_lowperm_layers = added_layers{2};
all_imperm_layers = added_layers{1};
all_added_layers = vertcat(added_layers{:});

s_trapped_imperm = trapped_cells{1};
s_trapped_lowperm = trapped_cells{2};

%% Compute rock+fluid objects
rock = makeRock(G, grid.Perm, poro);
T = computeTrans(G, rock, 'Verbose', true);
swr = 0.15;
snr = 0.2;

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
                           'smin', [swr, snr]);
             
s = linspace(0,1,100);                       
krW = fluid.krW(s).';
krO = fluid.krO(1-s).';
invalid_sw = find(s>1-snr);
invalid_so = find(s<swr);
krW(invalid_sw) = nan;
krO(invalid_so) = nan;

f1 = figure(1);
plot(s, krW, 'blue', s, krO, 'green', 'LineWidth', 1.5);
xlabel('Water saturation');
title('Relative permeability curves');
legend('krW', 'krO', 'Location', 'east');
saveas(f1, strcat(plot_dir, '/relperm'), 'png');
hold off
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, G.cells.num)';
p_e = 0.5*barsa;
p_cap = 3*barsa;

median_pc = 1*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 0;

if standard_pc
    plot_pc_lim = p_cap;
    fluid.pcOW = @(S, varargin) runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, G);
else
    plot_pc_lim = 5*median_pc;
    fluid.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, grid.Perm, baseperm, median_pc);
end

figure(1);
plot(dummy_Sw, fluid.pcOW(dummy_Sw), 'LineWidth', 1.5);
xlabel('Water saturation');
title('Capillary pressure function');
saveas(f1, strcat(plot_dir, '/cap_pres'), 'png');
hold off

%% Dummy well
well_h = 1; % cell perforations in vertical direction
perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
dummyW = addWell([], G, rock, perforation_idx, ...
        'Type', 'rate', 'Val', 1*meter^3/day(), ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');
    
%% Plot permeability field  
clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(grid.Perm, milli*darcy);
plotCellData(G, log10(perm_dummy), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn);
colorbarHist(log10(perm_dummy(all_added_layers)), [min(log10(perm_dummy)), max(log10(perm_dummy))], 'South', 51);
title('Log of permeability field');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
drawnow
hold off

saveas(f1, strcat(plot_dir, '/perm'), 'png');        

%% Plot structurally trapped cells
f2 = UtilFunctions.fullsizeFig(2);
plotGrid(G, all_lowperm_layers, 'FaceColor', 'yellow', 'EdgeColor', 'black', 'EdgeAlpha', 0.2, 'DisplayName', 'Lowperm layers');
plotGrid(G, all_imperm_layers, 'FaceColor', 'red', 'EdgeColor', 'black', 'EdgeAlpha', 0.2, 'DisplayName', 'Imperm layers');
plotGrid(G, s_trapped_lowperm, 'FaceColor', [0.25, 0.5, 0.25], 'EdgeColor', 'none', 'DisplayName', 'Lowperm trapped');
plotGrid(G, s_trapped_imperm, 'FaceColor', [0.5, 0, 0.5], 'EdgeColor', 'none', 'DisplayName', 'Imperm trapped');
plotGrid(G, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'none', 'DisplayName', 'Well');
title('Structurally trapped CO2', 'fontsize', 15);
axis equal tight
view([0, 0])
%zlim([min(z) max(z)]);
[hleg, hobj] = legend('Location', 'southoutside', 'Orientation', 'horizontal');
textobj = findobj(hobj, 'type', 'text');
set(hleg, 'position', [0.1, 0.06, 0.8, 0.05]);
set(textobj, 'fontsize', 12);
drawnow
hold off

saveas(f2, strcat(plot_dir, '/struct_trapped'), 'png');  


%% Set up solver
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
disp(model)

%% Boundary conditions, schedule
top_cells = G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = G.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, G, 'Right', pz, 'sat', [1 0]);

tot_time = 400*365*day();
%dt = rampupTimesteps(tot_time, 850*day(), 10);
dt = rampupTimesteps(round(tot_time/10), 500*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/10), 1000*day(), 0)']';
disp(numel(dt))

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(G, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

f4 = UtilFunctions.fullsizeFig(4); % to hold cap pressure

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, fluid.pcOW(state.s(:,1)), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colorbar; caxis([0, plot_pc_lim]);
title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f4, strcat(plot_dir, '/cap_pres_0'), 'png');

%% Run experiment
inj_stop_rate = 0.2;
% start with injection rate corresponding to 1/pv'th of total pore volume
pv = 30; % start with rate yielding total volume of 1/10th of pore volume in domain
rate = (sum(poreVolume(G, rock))/pv) / (inj_stop_rate*tot_time);

[max_volumes, leaked_boundary, structural_utilized, ...
    states, categorized_vols] = RunSimComputeTrapping(grid, rock, rate, state, model, ...
                                                                         s_trapped_imperm, s_trapped_lowperm, ...
                                                                         swr, snr, dt, bc, inj_stop_rate);


%% Show optimal solution found                                                       
disp(max_volumes)
dt_plot = cat(1, repmat(10, [fix(numel(states)/5), 1]), ...
                 repmat(20, [numel(states)-fix(numel(states)/5), 1]));
             
sat_CO2 = {}; % saturation in cells with non-zero CO2             

for i=1:numel(states)
    t = t + dt(i);
    
    sn = states{i}.s(:,2);
    sn = sn(sn > 1e-3);
    if ~isempty(sn)
        sat_CO2 = cat(1, sat_CO2, sn(sn > 1e-3)); % Buffer to avoid numerical artifacts
    else
        sat_CO2 = cat(1, sat_CO2, 0);
    end

    if ~mod(i, dt_plot(i))
        figure(3)
        set(f3, 'visible', 'off');
        plotCellData(G, states{i}.s(:,1), 'EdgeColor', 'none');
        plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
        colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
        title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
        axis equal tight;
        view([0, 0]);   

        set(f3, 'visible', 'on');
        filename_f3 = sprintf(strcat(plot_dir, '/sat_%d'), i);
        saveas(f3, filename_f3, 'png');
        
        figure(4);
        set(f4, 'visible', 'off');
        plotCellData(G, fluid.pcOW(states{i}.s(:,1)), 'EdgeColor', 'none');
        plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
        colorbar; caxis([0, plot_pc_lim]);
        title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
        axis equal tight
        view([0, 0])
       
        set(f4, 'visible', 'on');
        filename_f4 = sprintf(strcat(plot_dir, '/cap_pres_%d'), i);
        saveas(f4, filename_f4, 'png');
    end 
end

disp(leaked_boundary)


%% Store volume categories
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_seed_%d.mat'), seed.Seed);
struct_ratio_filename = sprintf(strcat(data_dir, '/struct_seed_%d.mat'), seed.Seed);
free_ratio_filename = sprintf(strcat(data_dir, '/free_seed_%d.mat'), seed.Seed);
sat_filename = sprintf(strcat(data_dir, '/sat_seed_%d.mat'), seed.Seed);

residual_vol = categorized_vols{1};
structural_vol = categorized_vols{2};
free_vol = categorized_vols{3};

save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
save(struct_ratio_filename, 'structural_vol');
save(free_ratio_filename, 'free_vol');

% STORE SATURATION FOR EACH FILE
save(sat_filename, 'sat_CO2');

%% Read data and store in structs
used_lab = {}; 

structural_struct = struct; % struct to hold structurally trapped values for each leaked percentage and num layer
residual_struct = struct;
free_struct = struct;
sat_struct = struct;
             
layer_folders = dir(strcat(data_base_dir, '/layers_*'));
layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
num_layers = numel(layer_folders);   

for j=1:num_layers
    lab_layers = regexp(layer_folders(j).name, 'layers_\d+', 'match'); 
    used_lab = cat(1, used_lab, lab_layers{1}); 

    structural_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/struct_*.mat'));
    free_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/free_*.mat'));
    residual_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/residual_*.mat'));
    sat_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/sat_*.mat'));

    structural_struct.imperm.(lab_layers{1}) = [];
    structural_struct.lowperm.(lab_layers{1}) = [];
    free_struct.(lab_layers{1}) = [];
    residual_struct.(lab_layers{1}) = [];
    sat_struct.(lab_layers{1}) = [];

    for k=1:numel(free_files)   
        get_seed = regexp(free_files(k).name, '\d+', 'match');
        get_seed = strcat('s', get_seed{1});
        
        load_structural = load(strcat(structural_files(k).folder, '\', structural_files(k).name), 'structural_vol');
        load_structural = load_structural.structural_vol;
        % NB: appended row-wise, so to access a specific seed do:
        % structural_struct.layers_X.lp_Y(i,:)
        structural_struct.imperm.(lab_layers{1}) = cat(2, structural_struct.imperm.(lab_layers{1}), load_structural.imperm);       
        structural_struct.lowperm.(lab_layers{1}) = cat(2, structural_struct.lowperm.(lab_layers{1}), load_structural.lowperm);
        load_free = load(strcat(free_files(k).folder, '\', free_files(k).name), 'free_vol');       
        free_struct.(lab_layers{1}) = cat(2, free_struct.(lab_layers{1}), load_free.free_vol);  
        load_residual = load(strcat(residual_files(k).folder, '\', residual_files(k).name), 'residual_vol');
        residual_struct.(lab_layers{1}) = cat(2, residual_struct.(lab_layers{1}), load_residual.residual_vol);
        
        % LOAD SATURATION FOR EACH FILE
        load_sat = load(strcat(sat_files(k).folder, '\', sat_files(k).name), 'sat_CO2'); 
        sat_struct.(lab_layers{1}).(get_seed) = load_sat.sat_CO2; % Append each seed for this num layer
    end

end


%% Plot ratio of total cells occupied by any CO2
f5 = figure(5);
clr = {'blue', 'red', 'green', 'magenta', 'orange'};
time = cumsum(dt)*second()/year();

for i=1:num_layers
    lab_layer = used_lab{i};
    sat_layer = struct2array(sat_struct.(lab_layer));
    n_states = numel(sat_layer);      

    sn_cells = cellfun(@numel, sat_layer) / G.cells.num; % ratio of total cells occupied by CO2     
    occupied_mean = mean(sn_cells, 2);
    occupied_std = std(sn_cells, 0, 2);   
   
    shadedErrorBar(time, sn_cells.', {@mean, @std}, 'lineprops', ...
                {'LineWidth', 1.5});
    
    hold on
end

xlabel('Time (years)')
ylabel('Ratio')
title('Ratio of total cells occupied by any CO2')
leg_names = strrep({used_lab{:}}, '_', ' ');
legend(leg_names, 'Location', 'southeast')
drawnow
saveas(f5, strcat(plot_base_dir, '/occupied_cells'), 'png');


%% TESTING


%% Functions
function pc = runStandardPc(S, dummy_S, swr, snr, p_e, p_cap, layers, G)
    pc_vals = Capillary.PcNew(dummy_S, swr, snr, p_e, p_cap, 2);

    region_table = {[dummy_S, zeros(numel(dummy_S), 1)], [dummy_S, pc_vals]}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg(region_table, S, region_idx);
end

function pc = runLeverettJ_2D(S, dummy_S, phi, K, dummy_K, K_base, layers, G)
    [grid_Sw, grid_K] = ndgrid(dummy_S, dummy_K);
    pc_vals = Capillary.LeverettJ(grid_Sw, phi, grid_K, K_base);
    
    region_table = {{grid_Sw, grid_K, zeros(size(grid_Sw))}, ...
                     {grid_Sw, grid_K,  pc_vals}}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg2D(region_table, S, K, region_idx);
end
