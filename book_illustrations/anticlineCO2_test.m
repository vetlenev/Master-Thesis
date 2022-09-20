%% Effect of number of lowperm layers
% Illustrate how different number of lowperm layers affect CO2 migration.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite UPR
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 7972;
% Parameters
n_layers = 30;
n_anticlines = ceil(n_layers/2);
n_straights = floor(n_layers/2);
corr_len_x = 5;

% Directories
plot_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/plotsAnticlines');
data_base_dir = strcat(ROOTDIR, '../master_thesis/book_illustrations/dataAnticlines');
data_layer_dir = sprintf(strcat(data_base_dir, '/layers_%d'), n_layers);
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d/corrlen_%d'), n_layers, corr_len_x);
data_dir = sprintf(strcat(data_base_dir, '/layers_%d/corrlen_%d'), n_layers, corr_len_x);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_layer_dir) ...
                & mkdir(data_dir) & mkdir(plot_dir);

% Set seed
seed = UtilFunctions.setSeed(data_layer_dir, my_seed);
rng(seed)

%% Define 2D grid
nx = 80; ny = 1; nz = 40; % 100 50
lx = 1000; ly = 1; lz = 400; % 1000 350
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

cons = {};

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

perm_layers = {};

% Make anticlines
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
        
    cons = cat(1, cons, [x0(:), z0(:)]);
end

% Make straights
% straights = {};
% unique_z = unique(z_dummy);
% for i=n_anticlines+1:n_anticlines+n_straights
%     x0 = x_start(i):max(lx/nx, lz/nz):x_stop(i);
%     x0 = x0(x0 >= min_x-dx & x0 <= max_x+dx);
%     z0 = perc_val(z_dummy, max(0.01, min(rand(), 0.90)));
%        
%     [min_z, min_zi] = min(abs(unique_z-z0)); 
%     z0 = unique_z(min_zi);
%     z0 = repmat(z0, numel(x0), 1).';
%  
%     cons = cat(1, cons, [x0(:), z0(:)]);
% end

cons = cons.';

%PG = pebiGrid2D(20, [lx, lz], 'cellConstraints', cons);
[PG, Pts, CCPts, F] = compositePebiGrid2D_vetle([lx/nx, lz/nz], [lx, lz], 'cellConstraints', cons);                        
                        %'polyBdr', [min_x, min_z; min_x, max_z; max_x, max_z; max_x, min_z]);

                     
PG = makeLayeredGrid(PG, 1); % Extrude pebi grid to 3D                    
PG = computeGeometry(PG); % to access centroids
PG.cells.indexMap = (1:PG.cells.num).';

disp(sum(cellfun(@(c) size(c,1), CCPts)))
disp(numel(find(PG.cells.tag == 1)))

x = PG.cells.centroids(:,1);
z = PG.cells.centroids(:,2); % z is second index since we extruded third dimension in y-direction
%z = PG.cells.centroids(:,3); % uncomment if rotated grid

%% Plot grid
clf;
figure(1)
plotGrid(PG,'FaceColor','none');
% plotLinePath(cons_faces,'--o','color',[0 .447 .741], ...
%         'linewidth',2,'MarkerFaceColor','w');

plotGrid(PG, PG.cells.tag == 1, 'facecolor', 'green')
%plotFaces(PG, PG.cells.faces(PG.faces.tag == 1))
axis tight

set(gca, 'YDir','reverse')

%% Define rock and fluid objects
poro = 0.3;

lowperm = 0.1*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [PG.cells.num 1]);

%corr_len_x = mean(x_stop - x_start) / 100;
corr_len_z = corr_len_x;

all_added_layers = PG.cells.indexMap(PG.cells.tag == 1); % lowperm anticlines

log_min = 0.01; log_max = 100;
layers_cells = cellfun(@(c) size(c,1), CCPts);
layers_totcells = sum(layers_cells);

log_field = exp(gaussianField([numel(all_added_layers) 1], ... % assume each layer only one cell thick in z-dir
                                   [log(log_min), log(log_max)], ...
                                   [fix(corr_len_x), 1, fix(corr_len_z)], 5));

[z_sorted_vals, z_sorted_idx] = sortrows(PG.cells.centroids(all_added_layers, :), 2);                             
perm(z_sorted_idx) = lowperm*log_field;                               

perm(~ismember(PG.cells.indexMap, all_added_layers)) = baseperm; % reset permeability for ambient rock
perm = perm.*sign(perm); % all positive perms
perm = max(perm, 1e-3*milli*darcy); % cap at non-zero to avoid singularity

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
rate = 0.3*meter^3/day;

perforation_idx = PG.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
W = addWell([], PG, rock, perforation_idx, ...
        'Type', 'rate', 'Val', rate, ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');   

% perforation_idx = PG.cells.indexMap(z <= max(z) & z >= max(z)-3*lz/nz & x == nx*dx);
% W = addWell([], PG, rock, perforation_idx, ...
%         'Type', 'rate', 'Val', rate, ...
%         'InnerProduct', 'ip_tpf', ...
%         'Radius', 0.1, ...
%         'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
%         'Name', 'P1'); 

% perforation_idx = PG.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x > max(x) - lx/5);    
% W = addWell(W, PG, rock, perforation_idx, ...
%         'Type', 'rate', 'Val', rate, ...
%         'InnerProduct', 'ip_tpf', ...
%         'Radius', 0.1, 'Dir', 'x', ...
%         'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
%         'Name', 'P1');   
%     

%% Plot permeability field  
clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(perm, milli*darcy);
plotCellData(PG, log10(perm_dummy), 'EdgeColor', 'black');
plotGrid(PG, W.cells, 'FaceColor', 'black');
%plotGrid(PG, W(1).cells, 'FaceColor', 'black')
%plotGrid(PG, W(2).cells, 'FaceColor', 'black')
colormap(flipud(jet));
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
p_e = 5*barsa;
p_cap = 15*barsa;

median_pc = 1*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 1;

if standard_pc
    plot_pc_lim = (p_cap+p_e)/2;
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

%% Boundary conditions
top_cells = PG.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = PG.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid.rhoWS * norm(gravity) * min(z);

bf = boundaryFaces(PG);
south_faces = bf(abs(PG.faces.normals(bf,1)) < eps & PG.faces.normals(bf,2) < 0); % NB: south-directed faces correspond to top since grid is simulated upside down 

bc = addBC(bc, south_faces, 'pressure', repmat(p_top, numel(south_faces), 1), 'sat', [1 0]);

east_faces = bf(PG.faces.normals(bf,1) < 0 & abs(PG.faces.normals(bf,2)) < eps);
east_neigh = PG.faces.neighbors(east_faces);
pz = fluid.rhoWS * norm(gravity) * z(east_neigh); % hydrostatic pressure at right boundary
%pz = flip(pz); % IMPORTANT: z-axis points upwards for PEBI-grid, so hydrostatic pressure must be reversed
bc = addBC(bc, east_faces, 'pressure', pz, 'sat', [1 0]);


clf;
figure(2)
plotFaces(PG, bf, 'edgecolor', 'red', 'linewidth', 1.5, 'DisplayName', 'No-flux')
plotFaces(PG, south_faces, 'edgecolor', 'blue', 'linewidth', 3, 'DisplayName', 'Top') 
plotFaces(PG, east_faces, 'edgecolor', 'green', 'linewidth', 3, 'DisplayName', 'Right')
title('Boundaries')
set(gca, 'YDir', 'reverse')
legend()

%% Schedule
tot_time = 600*365*day();
%dt = rampupTimesteps(tot_time, 850*day(), 10);
dt = rampupTimesteps(round(tot_time/10), 500*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/10), 1000*day(), 0)']';
disp(numel(dt))

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(PG, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotCellData(PG, state.s(:,1), 'EdgeColor', 'none');
plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.5);
plotGrid(PG, W.cells, 'FaceColor', 'black')
%plotGrid(PG, W(1).cells, 'FaceColor', 'black')
%plotGrid(PG, W(2).cells, 'FaceColor', 'black')
colormap(flipud(winter)); colorbar('southoutside'); caxis([swr, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
set(gca, 'YDir', 'reverse')
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

f4 = UtilFunctions.fullsizeFig(4); % to hold cap pressure

plotCellData(PG, fluid.pcOW(state.s(:,1)), 'EdgeColor', 'none');
plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.5);
plotGrid(PG, W.cells, 'FaceColor', 'black')
%plotGrid(PG, W(1).cells, 'FaceColor', 'black')
%plotGrid(PG, W(2).cells, 'FaceColor', 'black')
colorbar; caxis([0, plot_pc_lim]);
title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
set(gca, 'YDir', 'reverse')
drawnow

saveas(f4, strcat(plot_dir, '/cap_pres_0'), 'png');

%% Run experiment
inj_stop_rate = 0.2;
% start with injection rate corresponding to 1/pv'th of total pore volume
pv = 30; % start with rate yielding total volume of 1/30th of pore volume in domain
rate = (sum(poreVolume(PG, rock))/pv) / (inj_stop_rate*tot_time);

% [max_volumes, leaked_boundary, structural_utilized, ...
%     states, categorized_vols] = RunSimComputeTrapping(grid, rock, rate, state, model, ...
%                                                        s_trapped_imperm, s_trapped_lowperm, ...
%                                                        swr, snr, dt, bc, inj_stop_rate);
%       

times = cumsum(dt)/year();
n_steps = numel(dt);
[~, inj_stop] = min(abs(times - inj_stop_rate*times(end)));

schedule.control(2) = schedule.control(1); % new control for well shutoff
schedule.control(2).W.status = 0;
%schedule.control(2).W(1).status = 0; % shut off first well
%schedule.control(2).W(2).status = 0; % shut off second well
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway      

[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);                                                   


%% Show optimal solution found                                                       
dt_plot = cat(1, repmat(10, [fix(numel(states)/5), 1]), ...
                 repmat(20, [numel(states)-fix(numel(states)/5), 1]));
             
tot_vols = zeros(numel(states)+1, 1);
sim_vols = zeros(size(tot_vols));
leaked_vols = zeros(size(tot_vols));   
leaked_ratios = zeros(size(tot_vols));
sat_CO2 = zeros(size(tot_vols)); % saturation in cells with non-zero CO2

for i=1:numel(states)
    t = t + dt(i);
    
    sw = states{i}.s(:,1); 
    sn = 1 - sw; 
    
    sat_CO2(i+1) = numel(sn(sn > 1e-3));
    
    tot_vols(i+1) = sum(PG.cells.volumes.*sn)*mean(rock.poro);
    sim_vols(i+1) = min(rate*t, rate*sum(dt(1:inj_stop-1)));
      
    leaked_vols(i+1) = sim_vols(i+1) - tot_vols(i+1);
    leaked_ratios(i+1) = leaked_vols(i+1) / sim_vols(i+1);

    if ~mod(i, dt_plot(i))
        figure(3)
        set(f3, 'visible', 'off');
        plotCellData(PG, states{i}.s(:,1), 'EdgeColor', 'none');
        plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.5);
        plotGrid(PG, W.cells, 'FaceColor', 'black')
        %plotGrid(PG, W(1).cells, 'FaceColor', 'black')
        %plotGrid(PG, W(2).cells, 'FaceColor', 'black')
        colormap(flipud(winter)); colorbar('southoutside'); caxis([swr 1]);
        title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
        axis equal tight;
        set(gca, 'YDir', 'reverse')    

        set(f3, 'visible', 'on');
        filename_f3 = sprintf(strcat(plot_dir, '/sat_%d'), i);
        saveas(f3, filename_f3, 'png');
        
        figure(4);
        set(f4, 'visible', 'off');
        plotCellData(PG, fluid.pcOW(states{i}.s(:,1)), 'EdgeColor', 'none');
        plotGrid(PG, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.5);
        plotGrid(PG, W.cells, 'FaceColor', 'black')
        %plotGrid(PG, W(1).cells, 'FaceColor', 'black')
        %plotGrid(PG, W(2).cells, 'FaceColor', 'black')
        colorbar; caxis([0, plot_pc_lim]);
        title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
        axis equal tight
        set(gca, 'YDir', 'reverse')
        
        set(f4, 'visible', 'on');
        filename_f4 = sprintf(strcat(plot_dir, '/cap_pres_%d'), i);
        saveas(f4, filename_f4, 'png');
    end 
end


%% Store volume categories
leaked_ratio_filename = sprintf(strcat(data_dir, '/leaked_seed_%d.mat'), seed.Seed);
sat_filename = sprintf(strcat(data_dir, '/sat_seed_%d.mat'), seed.Seed);

save(leaked_ratio_filename, 'leaked_vols'); % save to compare for different nr low-perm layers
% STORE SATURATION FOR EACH FILE
save(sat_filename, 'sat_CO2');

%% Read data and store in structs
used_layer_lab = {}; 
used_corr_lab = struct; 

leaked_struct = struct;
sat_struct = struct;
             
layer_folders = dir(strcat(data_base_dir, '/layers_*'));
layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
num_layers = numel(layer_folders);   

for i=1:num_layers
    lab_layers = regexp(layer_folders(i).name, 'layers_\d+', 'match');
    corr_folders = dir(strcat(layer_folders(i).folder, '/', layer_folders(i).name, '/corrlen_*')); % all correlation folders for current num layer    
    used_layer_lab = cat(1, used_layer_lab, lab_layers{1}); 
    
    used_corr_lab.(lab_layers{1}) = {};
    
    for j=1:numel(corr_folders)
        lab_corr = regexp(corr_folders(j).name, 'corrlen_\d+', 'match');
        used_corr_lab.(lab_layers{1}) = cat(1, used_corr_lab.(lab_layers{1}), lab_corr{1});
    
        leaked_files = dir(strcat(corr_folders(j).folder, '/', corr_folders(j).name, '/leaked_*.mat'));
        sat_files = dir(strcat(corr_folders(j).folder, '/', corr_folders(j).name, '/sat_*.mat'));

        leaked_struct.(lab_layers{1}).(lab_corr{1}) = [];
        sat_struct.(lab_layers{1}).(lab_corr{1}) = [];

        for k=1:numel(leaked_files)   
            get_seed = regexp(leaked_files(k).name, '\d+', 'match');
            get_seed = strcat('s', get_seed{1});

            load_leaked = load(strcat(leaked_files(k).folder, '\', leaked_files(k).name), 'leaked_vols');
            leaked_struct.(lab_layers{1}).(lab_corr{1}) = cat(2, leaked_struct.(lab_layers{1}).(lab_corr{1}), ...
                                                                 load_leaked.leaked_vols);

            % LOAD SATURATION FOR EACH FILE
            load_sat = load(strcat(sat_files(k).folder, '\', sat_files(k).name), 'sat_CO2'); 
            sat_struct.(lab_layers{1}).(lab_corr{1}) = cat(2, sat_struct.(lab_layers{1}).(lab_corr{1}), ...
                                                              load_sat.sat_CO2); % Append each seed for this num layer
        end
    end
end


%% Plot ratio of total cells occupied by any CO2, for different number of layers
f5 = figure(5);
clr = {'blue', 'red', 'green', 'magenta', 'orange'};
time = cumsum(dt)*second()/year();

for i=1:num_layers
    lab_layer = used_layer_lab{i};
    sat_layer = struct2array(sat_struct.(lab_layer));
    %n_states = numel(sat_layer);      

    sn_cells = sat_layer / G.cells.num; % ratio of total cells occupied by CO2     
    occupied_mean = mean(sn_cells, 2);
    occupied_std = std(sn_cells, 0, 2);   
   
    shadedErrorBar([0; time], occupied_mean, occupied_std, 'lineprops', ...
                {'LineWidth', 1.5});
    
    hold on
end

xlabel('Time (years)')
ylabel('Ratio')
title('Ratio of total cells occupied by any CO2')
leg_names = strrep(used_layer_lab(:)', '_', ' ');
legend(leg_names, 'Location', 'southeast')
drawnow
saveas(f5, strcat(plot_base_dir, '/occupied_cells_all_layers'), 'png');


%% Same plot but for different correlation lengths, for each num layer
fi = 6;

for i=1:num_layers
    fig = figure(fi);
    lab_layer = used_layer_lab{i};
    corr_fields = fieldnames(sat_struct.(lab_layer));
    
    for k=1:numel(corr_fields)
        sat_corr = sat_struct.(lab_layer).(corr_fields{k});      

        sn_cells = sat_corr / G.cells.num; % occupied ratio for given correlation length    
        occupied_mean = mean(sn_cells, 2);
        occupied_std = std(sn_cells, 0, 2);   

        shadedErrorBar([0; time], occupied_mean, occupied_std, 'lineprops', ...
                    {'LineWidth', 1.5});

        hold on
    end
    
    xlabel('Time (years)')
    ylabel('Ratio')
    title({'CO2 occupation rate', strrep(lab_layer, '_', ': ')})
    leg_names = strrep(used_corr_lab.(lab_layer)(:)', '_', ' ');
    legend(leg_names, 'Location', 'southeast')
    drawnow
    saveas(fig, strcat(plot_base_dir, '/', lab_layer, '/occupied_cells'), 'png');
    
    fi = fi + 1;
end

%% Plot leaked volume for each num layer
fig = figure(fi);

for i=1:num_layers
    lab_layer = used_layer_lab{i};
    leaked_layer = struct2array(leaked_struct.(lab_layer));
     
    leaked_mean = mean(leaked_layer, 2);
    leaked_std = std(leaked_layer, 0, 2);   
   
    shadedErrorBar([0; time], leaked_mean, leaked_std, 'lineprops', ...
                {'LineWidth', 1.5});
    
    hold on
end

xlabel('Time (years)')
ylabel('Volume (m^3)')
title({'Leaked CO2 volume', sprintf('Injection rate: %d', round(convertTo(rate, meter^3/day), 2)) });
leg_names = strrep(used_layer_lab(:)', '_', ' ');
legend(leg_names, 'Location', 'southeast')
drawnow
saveas(fig, strcat(plot_base_dir, '/leaked_all_layers'), 'png');


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
