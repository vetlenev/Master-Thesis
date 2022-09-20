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
global G;
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

layers_length = randi([round(perc_val(x_dummy, 0.10)), ...
                        round(perc_val(x_dummy, 0.35))], ...
                        [n_layers, 1]);
                    
x_start = rand_x_mid - round(layers_length/2);                 
x_stop = rand_x_mid + round(layers_length/2);

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

PG = compositePebiGrid2D([lx/nx, lz/nz], [lx, lz], 'faceConstraints', cons_faces, ...
                        'cellConstraints', cons_cells, ...
                        'polyBdr', [min_x, min_z; min_x, max_z; max_x, max_z; max_x, min_z]);

% ------- UNCOMMENT THIS -----------
PG = makeLayeredGrid(PG, 1); % Extrude pebi grid to 3D   
% ----------------------------------

PG = computeGeometry(PG); % to access centroids
PG.cells.indexMap = (1:PG.cells.num).';

x = PG.cells.centroids(:,1);
z = PG.cells.centroids(:,2); % treat this as depth-dim

%% Plot grid
clf;
figure(1)
plotGrid(PG,'FaceColor','none');

imperm_layers = [];
for i=1:2
    cons_neighbors = PG.faces.neighbors(:,i);
    imperm_layers = cat(1, imperm_layers, cons_neighbors(PG.faces.tag == 1));
end

plotGrid(PG, imperm_layers, 'facecolor', 'red')
%plotFaces(PG, PG.cells.faces(PG.faces.tag == 1), 'red')
plotGrid(PG, PG.cells.tag == 1, 'facecolor', 'green')
axis tight
%set(gca, 'YDir','reverse')
