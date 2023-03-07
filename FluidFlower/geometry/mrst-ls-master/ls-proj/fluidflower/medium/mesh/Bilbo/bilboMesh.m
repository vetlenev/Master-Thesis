%% Create 2.5D extruded mesh
% Create a vertical cross section with a single fault and folded
% stratigraphy, extrude to have a 3D model with correct thickness.

clear, close all
mrstModule add upr
mrstVerbose on

mesh_case = 'bilbo';

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
if strcmp(mesh_case, 'bilbo')
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh/Bilbo/');
    meshpt  = load(fullfile(pth, 'bilbo_datapoints.mat'));
elseif strcmp(mesh_case, 'albus')
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh/');
    meshpt  = load(fullfile(pth, 'mesh_datapoints.mat'));
end
meshpt  = meshpt.stratiPoints;

% Depths will be positive
if strcmp(mesh_case, 'albus')
    meshpt.boundary(:,2) = meshpt.boundary(:,2)*-1;
    for n=1:length(meshpt.lines)
        if n <= numel(meshpt.wells)
            meshpt.wells{n}(:, 2) = meshpt.wells{n}(:, 2)*-1;
        end
        meshpt.lines{n}(:, 2) = meshpt.lines{n}(:, 2)*-1;
    end

    % plot
    fig2D(), hold on
    plot(meshpt.boundary(:,1), meshpt.boundary(:,2), 'k')
    plotLinePath(meshpt.lines(5:end), 'b');
    plotLinePath(meshpt.lines(1:4), 'r');
    plotLinePath(meshpt.wells(1:2), '.c', 'markerSize', 10);
    plotLinePath(meshpt.wells(3:6), '.k', 'markerSize', 10);
    plotLinePath(meshpt.wells(7:end), '.g', 'markerSize', 10);
    box on, axis equal tight, view([0, -90])
else
    fig2D(), hold on
    plot(meshpt.boundary(:,1), meshpt.boundary(:,2), 'k')
    plotLinePath(meshpt.lines(1:10), 'b');
    plotLinePath(meshpt.lines(11), 'r');
    plotLinePath(meshpt.wells(1:2), '.c', 'markerSize', 10);
    box on, axis equal tight, view([0, -90])
end


%% Generate 2D PEBI/composite grid
% We construct a 2D PEBI grid from the points using pebiGrid, adapted to
% the constraints imposed by the stratigraphy and fault
%--------------------------------------------------------------------------
% NOTE: For finer fault cells (between two relatively wide shear zone
% boundaries), add additional line(s) in between to constrain cell size!
%--------------------------------------------------------------------------
tic
% rng,  n,  FCFactor,   CCFactor,   CCEps,   ~h (mm), ~telaps (s), N
% 862,  120,    0.85,        0.5,   0.15,    18,      30,         12,653
% 407,  150,    0.9,         0.5,   0.15,    14,      63,         19,030
% 119,  200,    0.9,         0.5,   0.15,    11,     146,         33,099
% 588,  250,    0.8,         0.5,   0.15,     8,     350,         51,601
% 327,  500,    0.8,         0.5,   0.15,     8,     3412,       201,046
% 327, 1000,    0.8,         0.5,   0.15,     8,     +24h,     1,000,000
%rng(2335);
%nc  = 350; % Approximate number of cells in x-direction
L   = max(abs(meshpt.boundary));
f = 1/L(1);
%cellSize  = [f f]*1e-3; % f f will give 1mm after rescaling [xSize, ySize], in m
cellSize  = [4 4]*1e-3; % [xSize, ySize], in m
%meshpt.boundary = meshpt.boundary.*f;
%meshpt.lines = cellfun(@(x) x.*f, meshpt.lines, 'UniformOutput', false);
%meshpt.wells = cellfun(@(x) x.*f, meshpt.wells, 'UniformOutput', false);
% G2D = pebiGrid2D(max(L)/nc, L          , ...
%                  'polybdr', meshpt.boundary, ...       % Outline
%                  'faceConstraints', meshpt.lines, ...  % Fault lines
%                  'FCFactor'       , 1          , ...   % Relative size of fault cells
%                  'cellConstraints', meshpt.wells , ... % Well coordinates
%                  'CCRefinement'   , true         , ... % Refine
%                  'CCFactor'       , 0.5         ,  ... % Relative size of well cells
%                  'interpolateFC'  , false        , ... % Interpolate along fault lines
%                  'CCEps'          , 0.15*max(L)  );    % Refinement transition
% G2D = compositePebiGrid2D([max(L)/nc max(L)/nc], [L(1) L(2)], ...
%                  'polybdr', meshpt.boundary, ...       % Outline
%                  'faceConstraints', meshpt.lines, ...  % Fault lines
%                  'FCFactor'       , 1          , ... % Relative size of fault cells
%                  'cellConstraints', meshpt.wells , ... % Well coordinates
%                  'interpolateCC'  , false, ...
%                  'CCFactor'       , 1         ,  ... % Relative size of well cells
%                  'interpolateFC'  , false);   
% G2D = compositePebiGrid2D([max(L)/nc max(L)/nc], [L(1) L(2)], ...
%                  'polybdr', meshpt.boundary, ...       % Outline
%                  'faceConstraints', meshpt.lines, ...  % Fault lines
%                  'FCFactor'       , 1            , ... % Relative size of fault cells
%                  'cellConstraints', meshpt.wells , ... % Well coordinates
%                  'interpolateCC'  , false, ...
%                  'CCFactor'       , 0.5         ,  ... % Relative size of well cells
%                  'mlqtMaxLevel', 1, ...   
%                  'mlqtLevelSteps', [.0075], ...
%                  'interpolateFC'  , false); % Interpolate along fault lines 
G2D = compositePebiGrid2D(cellSize, [L(1) L(2)], ...
                 'polybdr', meshpt.boundary, ...       % Outline
                 'faceConstraints', meshpt.lines, ...  % Fault lines
                 'FCFactor'       , 1            , ... % Relative size of fault cells
                 'cellConstraints', meshpt.wells , ... % Well coordinates
                 'interpolateCC'  , false, ...
                 'CCFactor'       , 0.5        ,  ... % Relative size of well cells
                 'interpolateFC'  , true); % Interpolate along fault lines
G2D = removeShortEdges(G2D, 1e-04); % The grid may contain very short edges.
                                    % We remove these
%G2D.nodes.coords = G2D.nodes.coords./f;
toc
                                

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
colormap(jet)
for n=1:max(p)
   plotCellData(G2D, n*ones(sum(p==n), 1), p==n)
end
view([0 -90])

%% Extrude (2.5D)
% We construct a volumetric reservoir model by extruding the 2D grid

% Make layered grid
nLayers = 1;   % Number of layers
layerThickness = 1e-3*mean([10 11]);
G = makeLayeredGrid(G2D, 1);
G.nodes.coords(G.nodes.coords(:, 3) == 1, 3) = layerThickness;
G.nodes.coords = G.nodes.coords(:, [3 1 2]); % x -> y, y -> z, new dim = x
G = computeGeometry(G);

% Cells tag
G.cells.tag = repmat(G2D.cells.tag, nLayers, 1); % hz well cells (all layers)

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
% save(['G_composite_cellSize' num2str(cellSize(1)/f*1e3) '.mat'], 'G', 'G2D', 'p', ...
%      'compartID', 'wellNo', 'layerID')
save(['G_composite_cellSize' num2str(cellSize(1)*1e3) '.mat'], 'G', 'G2D', 'p', ...
     'compartID', 'wellNo', 'layerID')


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

