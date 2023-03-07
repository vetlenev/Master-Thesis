%% Create 2.5D extruded mesh
% Create a vertical cross section with a single fault and folded
% stratigraphy, extrude to have a 3D model with correct thickness.

clear, close all
mrstModule add upr
mrstVerbose on


%% Plotting functionality
fig2D = @() figure('Position', [0,0,800,500]);
fig3D = @() figure('Position', [0,0,1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );


%% Load data
% We start by loading a data structure that contains points describing the
% outline, faults, and well positions
pth = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/meshes/upr_mesh/');
meshpt  = load(fullfile(pth, 'test_datapoints.mat'));
meshpt  = meshpt.stratiPoints;

% rescale
l = max(meshpt.boundary(:,1));
meshpt.boundary = meshpt.boundary/l;
meshpt.wells{1} = meshpt.wells{1}/l;
for n=1:numel(meshpt.lines)
    meshpt.lines{n} = meshpt.lines{n}/l;
end

% plot
fig2D(), hold on
plot(meshpt.boundary(:,1), meshpt.boundary(:,2), 'k')
plotLinePath(meshpt.lines, 'b');
plotLinePath(meshpt.wells(1), '.r', 'markerSize', 10);
box on, axis equal tight, view([0, -90])


%% Generate 2D PEBI/composite grid
% We construct a 2D PEBI grid from the points using pebiGrid, adapted to
% the constraints imposed by the stratigraphy and fault
%--------------------------------------------------------------------------
% NOTE: For finer fault cells (between two relatively wide shear zone
% boundaries), add additional line(s) in between to constrain cell size!
%--------------------------------------------------------------------------
tic
% FCFactor,   CCFactor,   CCEps,   ~h (mm), ~telaps (s), N
% 0.25,       0.5,          -,       16,      83,        29095
L = max(abs(meshpt.boundary));
h = 0.005;    % m
cellSize = [h h];
nc_approx = round(L./cellSize);
G2D = compositePebiGrid2D(cellSize, [L(1) L(2)], ...
    'polybdr', meshpt.boundary, ...       % CellSize & Outline
    'faceConstraints', meshpt.lines, ...  % Fault lines
    'FCFactor'       , 0.5            , ... % Relative size of fault cells
    'cellConstraints', meshpt.wells , ... % Well coordinates
    'interpolateCC'  , false, ...
    'CCFactor'       , 0.5 ,  ...           % Relative size of well cells
    'interpolateFC'  , false);             % Interpolate along fault lines

G2D = removeShortEdges(G2D, 1e-05); % The grid may contain very short edges.
                                    % We remove these
toc
                  
% rescale
G2D.nodes.coords = G2D.nodes.coords*l;

%% Plot 2D grid
% Well cells are identified with G.cells.tag = true. We will need the well
% number later, so we find and store these.
wellNo2D                = nan(G2D.cells.num,1);
wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);
fig2D(), plotGrid(G2D); axis equal tight, box on                 % Grid
plotFaces(G2D, G2D.faces.tag, 'edgeColor', 'b', 'lineWidth', 2); % Lines
plotGrid(G2D, G2D.cells.tag , 'faceColor', 'r');                 % Wells
view([0 -90])


%% Find regions
% We identify these using functionality from the coarsegrid module
mrstModule add coarsegrid
p = ones(G2D.cells.num,1);
p = processPartition(G2D, p, find(G2D.faces.tag));

% Visualize regions
hold on
colormap(turbo)
for n=1:max(p)
   plotCellData(G2D, n*ones(sum(p==n), 1), p==n)
end
view([0 -90])

% Mesh 5.5 fault G add
% Load Mesh 6.5mm
% Find centroids of p=fault in 6.5mesh
% Find indices of minimum distance to cells in 5.5mesh
% If required, add remaining manually by plotting
% renumber p to 33 "regions"

%% Remove cells
% Run this section to remove silicone fault areas which are sealed
%remove_id = find(ismember(p, [21 26]));
remove_id = find(ismember(p, [22 27]));
[G2Drem, cellmap, facemap] = removeCells(G2D, remove_id);

% Update G2D to new mesh with correct cell tags for wells
wellid = find(G2D.cells.tag);
G2Drem_celltag = ismember(cellmap(1:G2Drem.cells.num), wellid);
G2Dfull = G2D;
G2D = G2Drem;
G2D.cells.tag = G2Drem_celltag;

% Find updated regions
p = ones(G2D.cells.num,1);
p = processPartition(G2D, p, find(G2D.faces.tag));

% Visualize regions
hold on
colormap(turbo)
for n=1:max(p)
   plotCellData(G2D, n*ones(sum(p==n), 1), p==n)
end
view([0 -90])

%wellNo2D
wellNo2D                = nan(G2D.cells.num,1);
wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);

% Plot grid with wells
fig2D(), plotGrid(G2D); axis equal tight, box on                 % Grid
plotFaces(G2D, G2D.faces.tag, 'edgeColor', 'b', 'lineWidth', 2); % Lines
plotGrid(G2D, G2D.cells.tag , 'faceColor', 'r');                 % Wells
view([0 -90])

%% Extrude (2.5D)
% We construct a volumetric reservoir model by extruding the 2D grid
thickness = 'variable';
% Make layered grid
nLayers = 1;   % Number of layers
layerThickness = 1e-3*round(mean([20 19 20 repelem(24, 4) 25 repelem(26, 5) ...
                            27 28]));
G = makeLayeredGrid(G2D, 1);
if strcmp(thickness, 'constant')
    G.nodes.coords(G.nodes.coords(:, 3) == 1, 3) = layerThickness;
elseif strcmp(thickness, 'variable')
    x = [0 0.03 0.73 1.33 2.13 2.83 2.86];
    ymax = max(G2D.nodes.coords(:,2));
    y = [ymax-0.03 ymax-0.33 ymax-0.63 ymax-0.93 ymax-1.23 ymax-1.53];
    [X, Y] = meshgrid(x, y);
    v = [0.019 0.019 0.02 0.019 0.02 0.019 0.019;
         0.019 0.019 0.024 0.024 0.024 0.019 0.019;
         0.019 0.019 0.026 0.028 0.026 0.019 0.019;
         0.019 0.019 0.027 0.025 0.026 0.019 0.019;
         0.019 0.019 0.026 0.026 0.024 0.019 0.019;
         0.019 0.019 0.019 0.023 0.02 0.019 0.019];
    h = h/1000;
    [Xq,Yq] = meshgrid(0:h:2.86,0:h:ymax);
    Vq = interp2(X,Y,v,Xq,Yq,'spline');
    surf(Xq,Yq,Vq,'edgecolor','none'); view([0 -90]);
    F = scatteredInterpolant(Xq(:),Yq(:),Vq(:));
    id1 = numel(G.nodes.coords(:,1))/2 + 1;
    vq = F(G.nodes.coords(id1:end,1), G.nodes.coords(id1:end,2));
    G.nodes.coords(G.nodes.coords(:, 3) == 1, 3) = vq;
end
G.nodes.coords = G.nodes.coords(:, [3 1 2]); % x -> y, y -> z, new dim = x
G           = computeGeometry(G);

% Cells tag
G.cells.tag = repmat(G2D.cells.tag, nLayers, 1); % hz well cells (all layers)

% Plot grid thickness
vqc = F(G.cells.centroids(:,2),G.cells.centroids(:,3));
fig3D(), plotCellData(G, vqc*1000,'linewidth', 0.1,'edgealpha',0.1); 
colormap(turbo)
axis equal tight, box on                   % Grid
c = colorbar; c.Label.String = '$T$ [mm]';
c.Label.FontSize = 12; c.Label.Interpreter = 'latex';
view([90 0])
ax = gca; ax.DataAspectRatio = [0.1 1 1];


% Faces tag
% This takes several minutes for fine grids (more than 30k cells), but the
% loop is necessary to avoid consuming all RAM with pdist. Find better way.

% G2D = computeGeometry(G2D);
% G0.faces.tag = false(G0.faces.num, 1);
% tagId = nan(G0.faces.num, 1);
% for n=1:G2D.faces.num 
%     disty = pdist2(G0.faces.centroids(:,2), G2D.faces.centroids(n, 1));
%     distz = pdist2( G0.faces.centroids(:,3), G2D.faces.centroids(n, 2));
%     dist = sqrt(disty.^2 + distz.^2);
%     tagId(n) = min(dist);
% end
% tagId = unique(tagId);
% G0.faces.tag(tagId) = repmat(G2D.faces.tag, nLayers, 1);

% IDs
wellNo       = repmat(wellNo2D, nLayers, 1);
layerID      = reshape(repmat(1:nLayers, G2D.cells.num, 1), [], 1);
compartID    = repmat(p,nLayers,1); 


%% Save data
% save(['G_benchmark_composite_cellSize' num2str(cellSize(1)*1e3) '_v2.mat'], 'G', 'G2D', 'p', ...
%       'compartID', 'wellNo', 'layerID')
save(['G_benchmark_composite_cellSize' num2str(cellSize(1)*1e3) ...
      '_Sremoved_ThickVar.mat'], 'G', 'G2D', 'p', 'compartID', 'wellNo', 'layerID')


%% Plot the resulting layered grid
% fig3D();
% Ncompart = max(p);
% cmap = copper(Ncompart);
% lyrs = {[1 9], [2 8], [4 7], 3, 5, 6}; % needs to be modified on a grid to grid basis
% colr = [8, 4, 7, 3, 1, 5];
% for n=1:numel(lyrs)
%     plotCellData(G, colr(n)*ones(sum(ismember(p, lyrs{n})), 1), ...
%                  ismember(p, lyrs{n}), 'edgealpha', 0.2)
% end
% outlineCoarseGrid(G, compartID,'EdgeColor','w','LineWidth',2);
% setAxProps(gca), %camlight();
% colormap(copper); c = colorbar; set(c, 'YTick', sort(colr));
% axis equal off
% %ylim([0 1]), zlim([0 1])

