%% Effect of number of lowperm layers
% Illustrate how different number of lowperm layers affect CO2 migration.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite UPR
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 8441;
n_layers = 5;
n_lowperm_layers = ceil(n_layers/2);
n_imperm_layers = floor(n_layers/2);

%% Define 2D grid
nx = 60; ny = 1; nz = 40; % 100 50
lx = 600; ly = 1; lz = 400; % 1000 350
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);
x_dummy = G.cells.centroids(:,1);
z_dummy = G.cells.centroids(:,3);
dz = max(diff(z_dummy))/2;
dx = max(diff(x_dummy))/2;


%% Test PEBI
sin_arg = @(x, x_start, x_stop) 0.8*(x - x_start)/(x_stop - x_start) + 0.2*(x_stop - x)/(x_stop - x_start);
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr); 
rand_a_b = @(a, b, num) (b-a).*rand(num,1) + a; % random number between a and b

cons_faces = {};
cons_cells = {};
rand_x_mid = rand_a_b(min(x_dummy), max(x_dummy), n_layers); 
%rand_x_mid = rand_x_mid(randperm(length(rand_x_mid)));

layers_length = randi([round(perc_val(x_dummy, 0.10)), ...
                        round(perc_val(x_dummy, 0.35))], ...
                        [n_layers, 1]);
                    
x_start = rand_x_mid - round(layers_length/2);                 
x_stop = rand_x_mid + round(layers_length/2);
%rand_z = linspace(perc_val(z_dummy, 0.05), perc_val(z_dummy, 0.9), n_layers);
min_x = min(x_dummy)-dx; max_x = max(x_dummy)+dx;
min_z = min(z_dummy)-dz; max_z = max(z_dummy)+dz;

for i=1:n_layers 
    r = (x_stop(i) - x_start(i))/2;
    x0 = x_start(i):max(lx/nx, lz/nz):x_stop(i);
    x0 = x0(x0 >= min_x-dx & x0 <= max_x+dx); % prevent layers going out of x-dim. Extra dx added to ensure we are indeed out of bounds
    
    amp = sin_arg(x0, x_start(i), x_stop(i));   
    
    z_curve = r*sin(pi*amp); % anticline form
    z_shift = r/2*(1+cos(pi/2-pi*min(amp))); % shift anticline
    rand_z = max(0.01, min(rand(), 0.90)); % 0.90 to prevent overlaying the well

    z0 = perc_val(z_dummy, rand_z) - z_curve + z_shift;
   
    z0 = z0(z0 >= min_z & z0 <= max_z); % prevent layers going out of z-dim
    x0 = x0(z0 >= min_z & z0 <= max_z);
    
    if i <= n_imperm_layers
        cons_faces = cat(1, cons_faces, [x0(:), z0(:)]);
    else
        cons_cells = cat(1, cons_cells, [x0(:), z0(:)]);
    end
end

cons_faces = transpose(cons_faces);
cons_cells = transpose(cons_cells);

%PG = pebiGrid2D(20, [lx, lz], 'cellConstraints', cons);
PG = compositePebiGrid2D([lx/nx, lz/nz], [lx, lz], 'faceConstraints', cons_faces, ...
                        'cellConstraints', cons_cells, ...
                        'polyBdr', [min_x, min_z; min_x, max_z; max_x, max_z; max_x, min_z]);


PG = makeLayeredGrid(PG, 1); % Extrude pebi grid to 3D                    
PG = computeGeometry(PG); % to access centroids
PG.cells.indexMap = (1:PG.cells.num).';
%PG.nodes.coords = PG.nodes.coords(:, [3 1 2]);

x = PG.cells.centroids(:,1);
z = PG.cells.centroids(:,2); % z is second index since we extruded third dimension in y-direction
%z = PG.cells.centroids(:,3);

%% Plot grid
clf;
figure(1)
plotGrid(PG,'FaceColor','none');
% plotLinePath(cons_faces,'--o','color',[0 .447 .741], ...
%         'linewidth',2,'MarkerFaceColor','w');

imperm_layers = [];
for i=1:2
    cons_neighbors = PG.faces.neighbors(:,i);
    imperm_layers = cat(1, imperm_layers, cons_neighbors(PG.faces.tag == 1));
end
%faces_top_layer = find(PG.faces.normals(:,3) ~= 0); % & find(PG.faces.neighbors(:,:) == 0)
%imperm_faces = find(ismember(faces_top_layer, find(PG.faces.tag==1)));

plotGrid(PG, imperm_layers, 'facecolor', 'red')
%hold on
%plotFaces(PG, PG.faces.tag == 1, 'red')
plotGrid(PG, PG.cells.tag == 1, 'facecolor', 'green')
axis tight
%view([90,0])
%set(gca, 'YDir','reverse')

%% Define rock and fluid objects
lowperm = 10*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [PG.cells.num 1]);
perms = {lowperm, baseperm};

lowperm_layers = PG.cells.indexMap(PG.cells.tag == 1); % lowperm anticlines
all_added_layers = vertcat(imperm_layers, lowperm_layers);

perm(lowperm_layers) = lowperm;
perm(imperm_layers) = 1e-3*milli*darcy;

poro = 0.3;

%% Directories

plot_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/plotsAnticlines');
data_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/dataAnticlines');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_layers);
data_dir = sprintf(strcat(data_base_dir, '/layers_%d'), n_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Set seed
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Compute rock+fluid objects
rock = makeRock(PG, perm, poro);
T = computeTrans(PG, rock, 'Verbose', true);
swr = 0.15;
snr = 0.2;

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
                           'smin', [swr, snr]);
                       
%% Dummy well
well_h = 1; % cell perforations in vertical direction
perforation_idx = PG.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
dummyW = addWell([], PG, rock, perforation_idx, ...
        'Type', 'rate', 'Val', 1*meter^3/day(), ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');   
    
%% Plot permeability field  
clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(perm, milli*darcy);
plotCellData(PG, log10(perm_dummy), 'EdgeColor', 'black');
plotGrid(PG, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'black');
colormap(autumn);
colorbarHist(log10(perm_dummy(all_added_layers)), [min(log10(perm_dummy)), max(log10(perm_dummy))], 'South', 51);
title('Log of permeability field');
axis equal tight
drawnow
set(gca, 'YDir','reverse')

saveas(f1, strcat(plot_dir, '/perm'), 'png');        
                       

%% Plot relperm                       
s = linspace(0,1,100);                       
krW = fluid.krW(s).';
krO = fluid.krO(1-s).';
invalid_sw = find(s>1-snr);
invalid_so = find(s<swr);
krW(invalid_sw) = nan;
krO(invalid_so) = nan;

f2 = figure(2);
plot(s, krW, 'blue', s, krO, 'green', 'LineWidth', 1.5);
xlabel('Water saturation');
title('Relative permeability curves');
legend('krW', 'krO', 'Location', 'east');
saveas(f2, strcat(plot_dir, '/relperm'), 'png');
hold off
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, PG.cells.num)';
p_e = 0.5*barsa;
p_cap = 3*barsa;

median_pc = 1*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 0;

if standard_pc
    plot_pc_lim = p_cap;
    fluid.pcOW = @(S, varargin) runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, PG);
else
    plot_pc_lim = 5*median_pc;
    fluid.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, perm, baseperm, median_pc);
end

clf;
figure(2);
plot(dummy_Sw, fluid.pcOW(dummy_Sw), 'LineWidth', 1.5);
xlabel('Water saturation');
title('Capillary pressure function');
saveas(f2, strcat(plot_dir, '/cap_pres'), 'png');
hold off  


%% Set up solver
gravity reset on

model = TwoPhaseOilWaterModel(PG, rock, fluid);
disp(model)

%% Boundary conditions, schedule
top_cells = PG.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = PG.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid.rhoWS * norm(gravity) * min(z);
%bc = pside(bc, PG, 'Top', p_top, 'sat', [1 0]);
top_faces = find(PG.faces.neighbors(:,2) == 0); % only check second index since only these have normal out of top
bc = addBC(bc, top_faces, 'pressure', p_top);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = addBC(bc, PG, 'pressure', pz);

tot_time = 400*365*day();
%dt = rampupTimesteps(tot_time, 850*day(), 10);
dt = rampupTimesteps(round(tot_time/10), 500*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/10), 1000*day(), 0)']';
disp(numel(dt))

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(PG, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(PG, state.s(:,1), 'EdgeColor', 'none');
plotGrid(PG, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

f4 = UtilFunctions.fullsizeFig(4); % to hold cap pressure

plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(PG, fluid.pcOW(state.s(:,1)), 'EdgeColor', 'none');
plotGrid(PG, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
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
rate = (sum(poreVolume(PG, rock))/pv) / (inj_stop_rate*tot_time);

% [max_volumes, leaked_boundary, structural_utilized, ...
%     states, categorized_vols] = RunSimComputeTrapping(grid, rock, rate, state, model, ...
%                                                        s_trapped_imperm, s_trapped_lowperm, ...
%                                                        swr, snr, dt, bc, inj_stop_rate);
%       
schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

times = cumsum(dt)/year();
[~, inj_stop] = min(abs(times - inj_stop_ratio*times(end)));
schedule.control(2) = schedule.control(1); % create second well
schedule.control(2).W.status = 0; % shut off second well
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway      

[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);
                                                   


%% Show optimal solution found                                                       
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
        plotCellData(PG, states{i}.s(:,1), 'EdgeColor', 'none');
        plotGrid(PG, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
        colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
        title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
        axis equal tight;
        view([0, 0]);   

        set(f3, 'visible', 'on');
        filename_f3 = sprintf(strcat(plot_dir, '/sat_%d'), i);
        saveas(f3, filename_f3, 'png');
        
        figure(4);
        set(f4, 'visible', 'off');
        plotCellData(PG, fluid.pcOW(states{i}.s(:,1)), 'EdgeColor', 'none');
        plotGrid(PG, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
        colorbar; caxis([0, plot_pc_lim]);
        title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
        axis equal tight
        view([0, 0])
       
        set(f4, 'visible', 'on');
        filename_f4 = sprintf(strcat(plot_dir, '/cap_pres_%d'), i);
        saveas(f4, filename_f4, 'png');
    end 
end


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
