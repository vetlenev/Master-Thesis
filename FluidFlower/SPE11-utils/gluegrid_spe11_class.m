%% Load geometry
filename = 'draft_spe11_with_facies_markers.geo';
[pts, loops, facies] = parse_spe11_geo(filename);
% pts (287 x 3): 3D coords of each Point
% loops (32 x variable): ordered array of points forming the surface of loop i
% facies (7 x variable): cell array of surfaces (set of loops) connected to each facie 

%% Other definitions
poly = struct; % to hold each polynomial
split_pts = struct; % to store split points of each poly
pts_overlap = struct;
nodes_overlap = struct;

logicalBottom = @(G) G.nodes.coords(:,2) == min(G.nodes.coords(:,2));
logicalTop = @(G) G.nodes.coords(:,2) == max(G.nodes.coords(:,2));
facesBottom = @(G) G.faces.centroids(:,2) == min(G.faces.centroids(:,2));
facesTop = @(G) G.faces.centroids(:,2) == max(G.faces.centroids(:,2));

gravity reset on

%% Global grid
[G_glob, x_glob, z_glob] = PolygonGrid.globalCartGrid(pts, 200, 100); % nx=200, nz=150

Lx = max(G_glob.faces.centroids(:,1));
Lz = max(G_glob.faces.centroids(:,2));
N = G_glob.cartDims;
Nx = N(1); Nz = N(2);

%% Organize polygons
polys = PolygonGrid.polyPoints(pts, loops, facies);
% Upper three surfaces: 28, 30, 32
poly.top = PolygonGrid(polys, 32);

ptop_orig = poly.top.p_orig;
ptop = poly.top.p;
poly.top.p_idx = 'p32';

%% Find top and bottom surfaces
[poly.top, split_pts.top] = reorderPts(poly.top, 32);

[top_side, bottom_side] = topAndBottomSurfaces(poly.top, split_pts.top, 32, pts_overlap);
poly.top.top_side = top_side;
poly.top.bottom_side = bottom_side;

poly.top = cartesianSubgrid(poly.top, Lx, Lz, Nx, Nz, []); % []: no conforming of horizontal grid size

%% dx-correction
% For each surface point, find closest x-coord in subgrid and change it to
% equal the surface coordinate
poly.top.bottom_mask = logicalBottom(poly.top.G); % only select nodes at bottom of grid
poly.top.top_mask = logicalTop(poly.top.G);
poly.top.G = computeGeometry(poly.top.G);
poly.top.faces_bottom = facesBottom(poly.top.G);
poly.top.faces_top = facesTop(poly.top.G);
% Bottom:
[poly.top, closest_bottom] = correct_dx_poly_new(poly.top, G_glob, poly.top.bottom_mask, 'bottom');
% Top:
[poly.top, closest_top] = correct_dx_poly_new(poly.top, G_glob, poly.top.top_mask, 'top');

%% Interpolation
% Interpolate z-values at surface of polygon
poly.top = interpolateZ_remaining_new(poly.top, closest_bottom, ...
                                poly.top.bottom_mask, 'bottom', 'linear'); % top face is the "remaining" parts of the polygon we want bottom surface
poly.top = interpolateZ_remaining_new(poly.top, closest_top, ...
                                poly.top.top_mask, 'top', 'linear');

% Interpolate x+z values in interior
poly.top = interpolateInternal(poly.top, poly.top.top_mask, poly.top.bottom_mask, []);

%% Finally, update grid
poly.top.G = computeGeometry(poly.top.G);
[ii, jj] = gridLogicalIndices(poly.top.G);
poly.top.G.i = ii;
poly.top.G.j = max(jj) - jj + 1; % To make lowest index start from top and increase downwards

Gtop = poly.top.G;
ptop_bottom = poly.top.bottom_side;
ptop_top = poly.top.top_side;

% Rename to polygon idx
poly = cell2struct(struct2cell(poly), {'p32'});

%% Create subgrid for neighboring polygon below
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 30, 32, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 18, 30, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 13, 18, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 28, 32, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 20, 28, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 15, 20, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 6, 15, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Create subgrids around pinch-outs
[poly, nodes_overlap, pts_overlap] = gluePinchOuts(polys, poly, 5, 13, 12, ...   
                                                     nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = gluePinchOuts(polys, poly, 29, 32, 25, ...
                                                     nodes_overlap, pts_overlap, G_glob);
% For polygon 19, upper neighbor contains pinch-out
[poly, nodes_overlap, pts_overlap] = gluePinchOuts(polys, poly, 19, 29, 25, ...
                                                     nodes_overlap, pts_overlap, G_glob);

%% Continue downwards until polygon 11 -> this one must be implemented manually
%For polygon 14, upper neighbor contains pinch-out -> handled internally in
%glueToUpperPolygon
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 14, 19, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 4, 14, ...
                                                        nodes_overlap, pts_overlap, G_glob, false); % Leftmost cells become very thin if inter_horz=false !

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 27, 6, ...
                                                        nodes_overlap, pts_overlap, G_glob, false);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 17, 6, ...
                                                        nodes_overlap, pts_overlap, G_glob, false);

%% Handle polygon 11
[poly, nodes_overlap, pts_overlap] = gluePolygon11(polys, poly, 11, [17,4,5], ...
                                                     nodes_overlap, pts_overlap, G_glob);

%% Continue downwards on both sides of bottom fault
[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 16, 27, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 10, 16, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 8, 10, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 3, 8, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Remaining polygons are glued to multiple sub-polygons
[poly, nodes_overlap, pts_overlap] = glueUpperPolygonMultiple(polys, poly, 7, 11, ...
                                                    nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = glueUpperPolygonMultiple(polys, poly, 2, 7, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Bottom polygon glued to structured polygon and fault
[poly, nodes_overlap, pts_overlap] = gluePolygon1(polys, poly, 1, [3,2], ...
                                                   nodes_overlap, pts_overlap, G_glob);


%% PEBI - top left fault
poly31_neighbors = cell(11, 3);
% List in counterclockwise order, starting from top left
poly31_neighbors{1,1} = poly.p32; poly31_neighbors{1,2} = 'top';
poly31_neighbors{2,1} = poly.p28; poly31_neighbors{2,2} = 'west';
poly31_neighbors{3,1} = poly.p20; poly31_neighbors{3,2} = 'west';
poly31_neighbors{4,1} = poly.p15; poly31_neighbors{4,2} = 'west';
poly31_neighbors{5,1} = poly.p6; poly31_neighbors{5,2} = 'west';
poly31_neighbors{6,1} = poly.p17; poly31_neighbors{6,2} = 'west';
poly31_neighbors{7,1} = poly.p11left; poly31_neighbors{7,2} = 'west';
poly31_neighbors{8,1} = poly.p4; poly31_neighbors{8,2} = 'east';
poly31_neighbors{9,1} = poly.p14; poly31_neighbors{9,2} = 'east';
poly31_neighbors{10,1} = poly.p19A; poly31_neighbors{10,2} = 'east';
poly31_neighbors{11,1} = poly.p29A; poly31_neighbors{11,2} = 'east';

[poly, nodes_overlap, pts_overlap] = glueFault_topleft(polys, poly, 31, poly31_neighbors, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Make PEBI subgrid
poly.p31PEBI = Faults.globalSubgrid(poly.p31PEBI, Lx, Lz, Nx, Nz, []);

%% In polygon
poly.p31PEBI = Faults.cellsInsideFault(poly.p31PEBI);

%% Set facies
[poly.p31PEBI, bneighbors] = Faults.assignFacies(poly.p31PEBI, poly31_neighbors);

%% Find discontinuous face transitions
p31GF = poly.p31PEBI.G_fault;
tip_faces = find(p31GF.faces.centroids(:,2) < 0.666 & ...
                   p31GF.faces.centroids(:,1) > 1.1065);
tip_faces = setdiff(tip_faces, boundaryFaces(p31GF));

edge_case_faces = find(p31GF.faces.centroids(:,1) > 1.087 & ...
                        p31GF.faces.centroids(:,1) < 1.09 & ...
                        p31GF.faces.centroids(:,2) > 0.735 & ...
                        p31GF.faces.centroids(:,2) < 0.737);  % somehow removeRedundantFaces does not work on this face ..., add it manually

tip_faces = [tip_faces; edge_case_faces];

poly.p31PEBI = Faults.removeRedundantFaces(poly.p31PEBI, bneighbors, tip_faces);

%% Fix faces unmatched with global grid
poly.p31PEBI = Faults.fixUnmatchedFaces(poly.p31PEBI);
p31GF = poly.p31PEBI.G_fault; % update

%% Plots / tests
figure()
xlim([0, 2.8])
ylim([0, 1.2])
colormap('jet')

poly_idxs = fieldnames(poly);
overlap_idxs = fieldnames(nodes_overlap);

if isfield(poly, 'p31PEBI')
   plotGrid(p31GF, 'EdgeAlpha', 0.1)
   plotCellData(p31GF, repmat(poly.p31PEBI.facies, p31GF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white')
   %plotFaces(p31G, poly.p31PEBI.G.faces.tag, 'EdgeColor', 'yellow', 'linewidth', 1.5)
   hold on
end

for i = 1:numel(poly_idxs) % plot subgrids first
    px = poly_idxs{i};
    if ~isempty(poly.(px).G) && ~strcmp(px, 'p31PEBI')
        plotGrid(poly.(px).G, 'EdgeAlpha', 0.3)
        plotCellData(poly.(px).G, repmat(poly.(px).facies, poly.(px).G.cells.num, 1), 'EdgeAlpha', 0.1)
        %plotCellData(poly.(px).G, poly.(px).G)
        hold on
    end    
    for j = 1:i-1%numel(overlap_idxs) % plot overlapping nodes
        px = overlap_idxs{j};
        plot(nodes_overlap.(px)(:,1), nodes_overlap.(px)(:,2), 'r.', 'markersize', 10)
        hold on
        plot(pts_overlap.(px)(:,1), pts_overlap.(px)(:,2), 'b.', 'markersize', 15)
        hold on
    %     plot(poly.(px).bottom_side(:,1), poly.(px).bottom_side(:,2), 'b.', 'markersize', 15)
    %     hold on
    end
    %pause(0.1)
end

%plot(all_nodes(:,1), all_nodes(:,2), 'g.', 'markersize', 10)

%% PEBI - top right fault
poly25_neighbors = cell(12, 2);
% List in counterclockwise order, starting from top left
poly25_neighbors{1,1} = poly.p32; poly25_neighbors{1,2} = 'top';
poly25_neighbors{2,1} = poly.p29B; poly25_neighbors{2,2} = {'west', 'top'};
poly25_neighbors{3,1} = poly.p29C; poly25_neighbors{3,2} = {'west', 'bottom'};
poly25_neighbors{4,1} = poly.p19B; poly25_neighbors{4,2} = {'west', 'top'};
poly25_neighbors{5,1} = poly.p19C; poly25_neighbors{5,2} = {'west', 'bottom'};
poly25_neighbors{6,1} = poly.p14; poly25_neighbors{6,2} = 'west';
poly25_neighbors{7,1} = poly.p4; poly25_neighbors{7,2} = 'west';
poly25_neighbors{8,1} = poly.p11right; poly25_neighbors{8,2} = 'east';
poly25_neighbors{9,1} = poly.p5A; poly25_neighbors{9,2} = 'east';
poly25_neighbors{10,1} = poly.p13; poly25_neighbors{10,2} = 'east';
poly25_neighbors{11,1} = poly.p18; poly25_neighbors{11,2} = 'east';
poly25_neighbors{12,1} = poly.p30; poly25_neighbors{12,2} = 'east';

[poly, nodes_overlap, pts_overlap] = glueFault_topright(polys, poly, 25, poly25_neighbors, ...
                                                        nodes_overlap, pts_overlap, G_glob);

%% Make PEBI subgrid
poly.p25PEBI = Faults.globalSubgrid(poly.p25PEBI, Lx, Lz, Nx, Nz, []);

%% In polygon
poly.p25PEBI = Faults.cellsInsideFault(poly.p25PEBI);

%% Plot grid
figure()
p25G = poly.p25PEBI.G;
bneighbors = unique(p25G.faces.neighbors(p25G.faces.tag, :));
plotGrid(p25G)
bcells = zeros(p25G.cells.num, 1);
bcells(bneighbors) = 1;
plotCellData(p25G, bcells)

figure()
p25GF = poly.p25PEBI.G_fault;
plotGrid(p25GF)

%% Set facies
[poly.p25PEBI, bneighbors] = Faults.assignFacies25(poly.p25PEBI, poly25_neighbors);

%% Find discontinuous face transitions
p25GF = poly.p25PEBI.G_fault;
tip_faces = find(p25GF.faces.centroids(:,2) < 0.673);
tip_faces = setdiff(tip_faces, boundaryFaces(p25GF));

edge_case_faces = find(p25GF.faces.centroids(:,1) > 1.433 & ...
                        p25GF.faces.centroids(:,1) < 1.44 & ...
                        p25GF.faces.centroids(:,2) > 0.988 & ...
                        p25GF.faces.centroids(:,2) < 0.9897);  % somehow removeRedundantFaces does not work on this face ..., add it manually
edge_case_faces = [edge_case_faces; find(p25GF.faces.centroids(:,1) > 1.245 & ...
                                        p25GF.faces.centroids(:,1) < 1.25 & ...
                                        p25GF.faces.centroids(:,2) > 0.904 & ...
                                        p25GF.faces.centroids(:,2) < 0.906)];  % somehow removeRedundantFaces does not work on this face ..., add it manually
edge_case_faces = [edge_case_faces; find(p25GF.faces.centroids(:,1) > 1.22 & ...
                                        p25GF.faces.centroids(:,1) < 1.225 & ...
                                        p25GF.faces.centroids(:,2) > 0.894 & ...
                                        p25GF.faces.centroids(:,2) < 0.896)];  % somehow removeRedundantFaces does not work on this face ..., add it manually


tip_faces = [tip_faces; edge_case_faces];
faces2keep = find(p25GF.faces.centroids(:,1) > 1.2085 & ...
                    p25GF.faces.centroids(:,1) < 1.2211 & ...
                    p25GF.faces.centroids(:,2) > 0.895 & ...
                    p25GF.faces.centroids(:,2) < 0.8975);

poly.p25PEBI = Faults.removeRedundantFaces(poly.p25PEBI, bneighbors, tip_faces, faces2keep);

%% Fix faces unmatched with global grid
point_set = {[1,2], [100,109], [109,115]};
poly.p25PEBI = Faults.fixUnmatchedFaces25(poly.p25PEBI, point_set);
p25GF = poly.p25PEBI.G_fault; % update

%% Change coordinate of degenerate face
bf = boundaryFaces(p25GF);
bf_mask = zeros(p25GF.faces.num, 1);
bf_mask(bf) = 1;

old_face = find(p25GF.faces.centroids(:,2) == min(p25GF.faces.centroids(~bf_mask,2)));
true_coord_left = poly.p25PEBI.bnodes(108,:);
true_coord_right = poly.p25PEBI.bnodes(110,:);

poly.p25PEBI = Faults.changeNodesForFaces(poly.p25PEBI, old_face, true_coord_left, true_coord_right);

%% Plots / tests
figure()
xlim([0, 2.8])
ylim([0, 1.2])
colormap('jet')

if isfield(poly, 'p31PEBI') && isfield(poly, 'p25PEBI')
   plotGrid(p31GF, 'EdgeAlpha', 0.1)
   plotCellData(p31GF, repmat(poly.p31PEBI.facies, p31GF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white')
   %plotFaces(p31G, poly.p31PEBI.G.faces.tag, 'EdgeColor', 'yellow', 'linewidth', 1.5)
   hold on
   plotGrid(p25GF, 'EdgeAlpha', 0.1)
   plotCellData(p25GF, repmat(poly.p25PEBI.facies, p25GF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white') 
   hold on
end

for i = 1:numel(poly_idxs) % plot subgrids first
    px = poly_idxs{i};
    if ~isempty(poly.(px).G) && ~strcmp(px, 'p31PEBI') && ~strcmp(px, 'p25PEBI')
        plotGrid(poly.(px).G, 'EdgeAlpha', 0.3)
        plotCellData(poly.(px).G, repmat(poly.(px).facies, poly.(px).G.cells.num, 1), 'EdgeAlpha', 0.1)
        %plotCellData(poly.(px).G, poly.(px).G)
        hold on
    end    
    for j = 1:i-1%numel(overlap_idxs) % plot overlapping nodes
        px = overlap_idxs{j};
        plot(nodes_overlap.(px)(:,1), nodes_overlap.(px)(:,2), 'r.', 'markersize', 10)
        hold on
        plot(pts_overlap.(px)(:,1), pts_overlap.(px)(:,2), 'b.', 'markersize', 15)
        hold on
    end
end


%% PEBI - bottom fault
polyBF_nums = [24, 9, 21, 26, 22, 23]; % polyBF = poly bottom fault

polyBF_neighbors.p24 = cell(6,2);
polyBF_neighbors.p9 = cell(1,2);
polyBF_neighbors.p21 = cell(5,2);
polyBF_neighbors.p26 = cell(2,2);
polyBF_neighbors.p22 = cell(1,2);
polyBF_neighbors.p23 = cell(3,2);
% poly 24
polyBF_neighbors.p24{1,1} = poly.p6; polyBF_neighbors.p24{1,2} = 'top';
polyBF_neighbors.p24{2,1} = poly.p27; polyBF_neighbors.p24{2,2} = 'west';
polyBF_neighbors.p24{3,1} = poly.p16; polyBF_neighbors.p24{3,2} = 'west';
polyBF_neighbors.p24{4,1} = poly.p10; polyBF_neighbors.p24{4,2} = 'west';
polyBF_neighbors.p24{5,1} = poly.p11left; polyBF_neighbors.p24{5,2} = 'east';
polyBF_neighbors.p24{6,1} = poly.p17; polyBF_neighbors.p24{6,2} = 'east';
% poly 9
polyBF_neighbors.p9{1,1} = poly.p11left; polyBF_neighbors.p9{1,2} = 'east';
% poly 21
polyBF_neighbors.p21{1,1} = poly.p8; polyBF_neighbors.p21{1,2} = 'west';
polyBF_neighbors.p21{2,1} = poly.p2left; polyBF_neighbors.p21{2,2} = 'east';
polyBF_neighbors.p21{3,1} = poly.p7small; polyBF_neighbors.p21{3,2} = {'top', 'east', 'bottom'};
polyBF_neighbors.p21{4,1} = poly.p7left; polyBF_neighbors.p21{4,2} = 'east';
polyBF_neighbors.p21{5,1} = poly.p11left; polyBF_neighbors.p21{5,2} = 'east';
% poly 26
polyBF_neighbors.p26{1,1} = poly.p8; polyBF_neighbors.p26{1,2} = 'west';
polyBF_neighbors.p26{2,1} = poly.p3; polyBF_neighbors.p26{2,2} = 'west';
% poly 22
polyBF_neighbors.p22{1,1} = poly.p2left; polyBF_neighbors.p22{1,2} = 'east';
% poly 23
polyBF_neighbors.p23{1,1} = poly.p3; polyBF_neighbors.p23{1,2} = 'west';
polyBF_neighbors.p23{2,1} = poly.p1mid; polyBF_neighbors.p23{2,2} = 'bottom';
polyBF_neighbors.p23{3,1} = poly.p2left; polyBF_neighbors.p23{3,2} = 'east';

polyBF_norder = cell(17,3);
for i=1:4
    polyBF_norder{i,1} = polyBF_neighbors.p24{i,1};
    polyBF_norder{i,2} = polyBF_neighbors.p24{i,2};
    polyBF_norder{i,3} = 'p24';
end
polyBF_norder{5,1} = polyBF_neighbors.p21{1,1}; 
polyBF_norder{5,2} = polyBF_neighbors.p21{1,2};
polyBF_norder{5,3} = 'p21';
for i=6:7
    polyBF_norder{i,1} = polyBF_neighbors.p26{i-5,1};
    polyBF_norder{i,2} = polyBF_neighbors.p26{i-5,2};
    polyBF_norder{i,3} = 'p26';
end
for i=8:10
    polyBF_norder{i,1} = polyBF_neighbors.p23{i-7,1};
    polyBF_norder{i,2} = polyBF_neighbors.p23{i-7,2};
    polyBF_norder{i,3} = 'p23';
end
polyBF_norder{11,1} = polyBF_neighbors.p22{1,1}; 
polyBF_norder{11,2} = polyBF_neighbors.p22{1,2};
polyBF_norder{11,3} = 'p22';
for i=12:15
    polyBF_norder{i,1} = polyBF_neighbors.p21{i-10,1}; 
    polyBF_norder{i,2} = polyBF_neighbors.p21{i-10,2};
    polyBF_norder{i,3} = 'p21';
end
polyBF_norder{16,1} = polyBF_neighbors.p9{1,1}; 
polyBF_norder{16,2} = polyBF_neighbors.p9{1,2};
polyBF_norder{16,3} = 'p9';
for i=17:18
    polyBF_norder{i,1} = polyBF_neighbors.p24{i-12,1}; 
    polyBF_norder{i,2} = polyBF_neighbors.p24{i-12,2};
    polyBF_norder{i,3} = 'p24';
end

[poly, external_nodes] = bottomFaultBoundary(polys, poly, polyBF_nums, polyBF_norder, ...
                                                 nodes_overlap, pts_overlap, G_glob);

%% Then separate out individual polygons in fault
poly.pBFPEBI = Faults.globalSubgrid(poly.pBFPEBI, Lx, Lz, Nx, Nz, []);

%% Plot bottom fault
figure()
pBFG = poly.pBFPEBI.G;
bneighbors = unique(pBFG.faces.neighbors(pBFG.faces.tag, :));
bcells = zeros(pBFG.cells.num, 1);
bcells(bneighbors) = 1;
plotGrid(pBFG);
plotCellData(pBFG, bcells)

%% Fixing ...
poly.pBFPEBI = Faults.cellsInsideFault(poly.pBFPEBI);
[poly.pBFPEBI, bneighbors] = Faults.assignFacies(poly.pBFPEBI, polyBF_norder);
pBFGF = poly.pBFPEBI.G_fault;

%% New
tip_faces = find(pBFGF.faces.centroids(:,1) > 0.53);
tip_faces = setdiff(tip_faces, boundaryFaces(pBFGF));

edge_case_faces = find(pBFGF.faces.centroids(:,1) > 0.485 & ...
                        pBFGF.faces.centroids(:,1) < 0.490 & ...
                        pBFGF.faces.centroids(:,2) > 0.746);

edge_case_faces = [edge_case_faces; find(pBFGF.faces.centroids(:,1) > 0.478 & ...
                                        pBFGF.faces.centroids(:,2) > 0.615 & ...
                                        pBFGF.faces.centroids(:,2) < 0.618)];

tip_faces = [tip_faces; edge_case_faces];

poly.pBFPEBI = Faults.removeRedundantFaces(poly.pBFPEBI, bneighbors, tip_faces);
pBFGF = poly.pBFPEBI.G_fault;

%% Plots / tests
figure()
xlim([0, 2.8])
ylim([0, 1.2])
colormap('jet')

if isfield(poly, 'p31PEBI') && isfield(poly, 'p25PEBI') && isfield(poly, 'pBFPEBI')
   plotGrid(p31GF, 'EdgeAlpha', 0.1)
   plotCellData(p31GF, repmat(poly.p31PEBI.facies, p31GF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white')
   %plotFaces(p31G, poly.p31PEBI.G.faces.tag, 'EdgeColor', 'yellow', 'linewidth', 1.5)
   hold on
   plotGrid(p25GF, 'EdgeAlpha', 0.1)
   plotCellData(p25GF, repmat(poly.p25PEBI.facies, p25GF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white') 
   hold on
   plotGrid(pBFGF, 'EdgeAlpha', 0.1)
   plotCellData(pBFGF, repmat(poly.pBFPEBI.facies, pBFGF.cells.num, 1), 'EdgeAlpha', 0.7, 'EdgeColor', 'white')
   hold on
end

for i = 1:numel(poly_idxs) % plot subgrids first
    px = poly_idxs{i};
    if ~isempty(poly.(px).G) && ~strcmp(px, 'p31PEBI') && ~strcmp(px, 'p25PEBI')
        plotGrid(poly.(px).G, 'EdgeAlpha', 0.3)
        plotCellData(poly.(px).G, repmat(poly.(px).facies, poly.(px).G.cells.num, 1), 'EdgeAlpha', 0.1)
        %plotCellData(poly.(px).G, poly.(px).G)
        hold on
    end    
    for j = 1:i-1%numel(overlap_idxs) % plot overlapping nodes
        px = overlap_idxs{j};
        plot(nodes_overlap.(px)(:,1), nodes_overlap.(px)(:,2), 'r.', 'markersize', 10)
        hold on
        plot(pts_overlap.(px)(:,1), pts_overlap.(px)(:,2), 'b.', 'markersize', 15)
        hold on
    end
end




%% Make top polygon basis for gluing
assert(strcmp(poly_idxs{1}, 'p32')) % make sure top subgrid is first element
poly_idxs = poly_idxs(2:end);

%% Nested gluing
poly.glued = poly.p32; % top grid is basis for gluing

% NB: Gluing doesn't work until faults are set up!
for i=1:numel(poly_idxs)    
    px = poly_idxs{i};
    disp(px)
    G_new = poly.(px).G;
    G_glued = glue2DGrid(poly.glued.G, G_new);

    G_glued.cells.indexMap = [poly.glued.G.cells.indexMap; G_new.cells.indexMap + poly.glued.G.cells.num];
    G_glued.i = [poly.glued.G.i; G_new.i];
    G_glued.j = [poly.glued.G.j; G_new.j];

    poly.glued.G = G_glued;  
end


%% FUNTCION DEFINTIONS
function [poly_obj, nodes_overlap, pts_overlap] = glueToUpperPolygon(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)
    
    if nargin > 7
        inter_horz = varargin{1};
    else        
        inter_horz = true;
    end

    p_idx = strcat('p', string(poly_num));
    
    p_name = fieldnames(poly_obj); 
    is_pinch = 0;
    for i=1:numel(p_name)
        is_pinch = is_pinch + ~isempty(regexp(p_name{i}, strcat('p', string(poly_num_upper), '[ABC]'), 'match'));
    end

    if poly_num == 2 % handle this separately
        p_idx_upper = {strcat('p', string(poly_num_upper), 'left'), ...
                        strcat('p', string(poly_num_upper), 'mid'), ...
                        strcat('p', string(poly_num_upper), 'right')};
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2}), ...
                       poly_obj.(p_idx_upper{3})};

    elseif is_pinch
        p_idx_upper = {strcat('p', string(poly_num_upper), 'A'), ...
                        strcat('p', string(poly_num_upper), 'C')}; % subgrid A and C together comprise the bottom side of upper polygon
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2})};
    else
        p_idx_upper = strcat('p', string(poly_num_upper));
        poly_upper = poly_obj.(p_idx_upper);
    end

    poly_obj.(p_idx) = PolygonGrid(all_polys, poly_num);           
    poly = poly_obj.(p_idx); % dummy variable for readability    
    poly.p_idx = p_idx;    

    for i=1:numel(poly_upper)
        % if points have been added to neighboring polygon, make sure to
        % add these to current polygon
        if ~iscell(poly_upper)
            poly_i = poly_upper;
        else
            poly_i = poly_upper{i};
        end
        if ~isempty(poly_i.p_added)
            n_added = size(poly_i.p_added, 1);
            poly.p(end+1:end+n_added, :) = poly_i.p_added;
        end
    end

    Lx = max(G_glob.faces.centroids(:,1));
    Lz = max(G_glob.faces.centroids(:,2));
    N = G_glob.cartDims;
    Nx = N(1); Nz = N(2);   
    
    if numel(poly_upper) > 1 % Upper neighbor contains a pinch-out
        %[nodes_overlap.(p_idx), pts_overlap.(p_idx)] = PinchOuts.findOverlappingNodesPinch(poly, poly_upper{1}, poly_upper{2}, 'top');
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');
    else
        [nodes_overlap.(p_idx), pts_overlap.(p_idx)] = findOverlappingNodes(poly, poly_upper, 'top'); 
    end
    
    num_x_overlap = size(nodes_overlap.(p_idx), 1);

    [poly, split_pts] = reorderPts(poly, poly_num); % two last args required for edge-cases
    
    [top_side, bottom_side] = topAndBottomSurfaces(poly, split_pts, poly_num, pts_overlap); % poly_num required for edge-cases
    
    poly.top_side = unique(top_side, 'rows');
    poly.bottom_side = unique(bottom_side, 'rows');
    
    poly = cartesianSubgrid(poly, Lx, Lz, Nx, Nz, num_x_overlap);                
    
    poly.G.nodes.coords(poly.top_mask, :) = nodes_overlap.(p_idx); % glue to overlapping nodes on top side
    
    [poly, closest_bottom] = correct_dx_poly_new(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
    
    poly = interpolateZ_remaining_new(poly, closest_bottom, ...
                                    poly.bottom_mask, 'bottom', 'linear');
           
    poly = interpolateInternal(poly, poly.top_mask, poly.bottom_mask, []);

    if ~isempty(poly.p_we)
        poly = interpolateSide(poly);
    end
    % Interpolate x-points in internal to conform with shifted boundary                         
    if inter_horz
        poly = interpolateHorizontal(poly);
    else
        poly = interpolatePartlyHorizontal(poly, 0.2); % on west and east side, interpolate 20% of horizontal extent
    end

    % Get logical indices for new polygon   
    poly = logicalIndicesUpperOverlap(poly, poly_upper);  

    % --- Is this needed? ---
    poly.faces_bottom = poly.G.faces.centroids(:,2) == min(poly.G.faces.centroids(:,2)); 
    poly.faces_top = poly.G.faces.centroids(:,2) == max(poly.G.faces.centroids(:,2));
    % ---
    
    poly_obj.(p_idx) = poly;
end


function [poly_obj, nodes_overlap, pts_overlap] = gluePinchOuts(all_polys, poly_obj, poly_num, poly_num_upper, poly_num_pinch, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)   
    p_idx = strcat('p', string(poly_num));   

    p_name = fieldnames(poly_obj); 
    is_pinch = 0;
    for i=1:numel(p_name)
        is_pinch = is_pinch + ~isempty(regexp(p_name{i}, strcat('p', string(poly_num_upper), '[ABC]'), 'match'));
    end

    if is_pinch
        p_idx_upper = {strcat('p', string(poly_num_upper), 'A'), ...
                        strcat('p', string(poly_num_upper), 'C')}; % subgrid A and C together comprise the bottom side of upper polygon
        poly_upper = {poly_obj.(p_idx_upper{1}), ...
                       poly_obj.(p_idx_upper{2})};
    else
        p_idx_upper = strcat('p', string(poly_num_upper));
        poly_upper = poly_obj.(p_idx_upper);
    end
    p_idx_pinch = strcat('p', string(poly_num_pinch));

    poly = PolygonGrid(all_polys, poly_num); % poly_obj.(p_idx) = ...
    poly_pinch = PolygonGrid(all_polys, poly_num_pinch);
    %poly = poly_obj.(p_idx); % dummy variable for readability    
    p = poly.p;
    pp = poly_pinch.p;

    % Create separation points at end of pinch-out   
    switch poly_num
        case 5
            pin_idx = 13;
        case 29
            pin_idx = 12;
        case 19
            pin_idx = 15;
    end
    pinch = p(pin_idx, :);
    z_pinch = p(pin_idx, 2);
    x_pinch = p(pin_idx, 1);

    Lx_glob = max(G_glob.faces.centroids(:,1));
    Lz_glob = max(G_glob.faces.centroids(:,2));
    N = G_glob.cartDims;
    Nx_glob = N(1); Nz_glob = N(2);          

    if numel(poly_upper) > 1 % Upper neighbor contains a pinch-out
        [nodes_overlap_pX, pts_overlap_pX] = PinchOuts.findOverlappingNodesMultiple(poly, poly_upper, 'top');
    else
        [nodes_overlap_pX, pts_overlap_pX] = findOverlappingNodes(poly, poly_upper, 'top');
    end
    top_nodes = nodes_overlap_pX;
    top_pts = pts_overlap_pX;
    num_x_overlap = size(top_nodes, 1);
        
    % Find corresponding node in bottom side     
    switch poly_num
        case 5
            bottom_pts = [p(26:end,:); p(1:8,:)];
        case 29
            bottom_pts = p(2:6, :);
        case 19
            bottom_pts = p(8:11, :);
    end
    poly.bottom_side = bottom_pts;
    poly.top_side = top_pts;
    % Create (dummy) subgrid for full polygon
    poly = cartesianSubgrid(poly, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_overlap);           

    [bottom_nodes, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(poly, poly.bottom_side, G_glob);
    poly.bottom_side_new = poly.bottom_side;
    poly.top_side_new = poly.top_side;
    poly.bottom_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(poly.bottom_side_new, bottom_nodes, closest_mask, 'linear'); % changed first argument from poly.bottom_side to poly.bottom_side_new
    bottom_nodes(rem_idx,2) = z_new;

    % Find top/bottom-node closest to tip of pinch
    % --- changeTPL ---
    if true
        dist_pinch = @(node) sqrt((node(:,1)-x_pinch).^2 + (node(:,2)-z_pinch).^2);
        top2pinch = dist_pinch(top_nodes);
        bottom2pinch = dist_pinch(bottom_nodes);
        net_dist = top2pinch + bottom2pinch;
        [~, pin_idx] = min(net_dist);
    else
        [~, pin_idx] = min(abs(top_nodes(:,1) - x_pinch));        
        pin_idx = nnz(top_nodes(top_nodes(:,1) <= sep_point_top(:,1), 1)); % this is index for separation point at bottom            
    end
    sep_point_top = top_nodes(pin_idx,:);
    sep_point_bottom = bottom_nodes(pin_idx,:);
    % -----------------    

    % Using separation points, divide polygon into pXA, pXB and pXC (for
    % polygon of number X):
    x_top = top_nodes(:,1);
    x_bottom = bottom_nodes(:,1);    


    % --- A: Make catesian subgrid of pXA based on nodes and points to left
    % of separation point ---
    top_nodesA = top_nodes(x_top <= sep_point_top(:,1), :);
    bottom_nodesA = bottom_nodes(x_bottom <= sep_point_bottom(:,1), :);        
    num_x_pXA = size(top_nodesA, 1);

    top_ptsA = [top_pts(top_pts(:,1) <= sep_point_top(:,1), :); sep_point_top]; % polygon data points on top side LEFT of separation point, including top separation point
    bottom_ptsA = [bottom_pts(bottom_pts(:,1) <= sep_point_bottom(:,1), :); sep_point_bottom]; % (:,1) -> compare x-value
    side_ptsA = [top_ptsA; bottom_ptsA];    

    polyA = PinchOuts(all_polys, poly_num, side_ptsA);
    polyA.bottom_side = bottom_ptsA;
    polyA.top_side = top_ptsA;
    polyA = cartesianSubgrid(polyA, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXA);     

    % CHANGE COORDS (pXA)
    [polyA, z_sep_idx] = PinchOuts.coordCorrectionSubgridA(polyA, top_nodesA, bottom_nodesA, pinch);
    %polyA = interpolateHorizontal(polyA);
   
    % Get logical indices for new polygon
    polyA = logicalIndicesUpperOverlap(polyA, poly_upper);     

    p_idxA = strcat(p_idx, 'A');
    polyA.p_idx = p_idxA;
    poly_obj.(p_idxA) = polyA;
    nodes_overlap.(p_idxA) = top_nodesA; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxA) = top_ptsA;

    % Get nodes of east side (needed when gluing pXB and pXC)
    NxA = polyA.G.cartDims(1)+1; NzA = polyA.G.cartDims(2)+1; % +1 since we are indexing NODES not CELLS
    east_nodesA = polyA.G.nodes.coords(NxA:NxA:NxA*NzA, :); % correct indexing ???
    east_nodesAC = east_nodesA(1:z_sep_idx, :); % from bottom to pinch
    east_nodesAB = east_nodesA(z_sep_idx:end, :); % from pinch to top


    % --- B: Make cartesian subgrid of pXB. ---
    top_nodesB = top_nodes(x_top >= sep_point_top(:,1), :);       
    [bottom_ptsB, top_ptsC] = PinchOuts.getSides_pinchout(poly_pinch, poly_num_pinch, poly_num); % bottom points of polygon pXB are top points of fault, and top points of pXC are bottom point of fault
    num_x_pXB = size(top_nodesB, 1);  

    top_ptsB = [sep_point_top; top_pts(top_pts(:,1) >= sep_point_top(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsB = [top_ptsB; bottom_ptsB];

    polyB = PinchOuts(all_polys, poly_num, side_ptsB);       
    polyB.bottom_side = bottom_ptsB;
    polyB.top_side = top_ptsB;
    NzB = NzA - z_sep_idx; % NB: number of CELLS not NODES
    polyB = cartesianSubgrid(polyB, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXB, NzB);
   
    [bottom_nodesB, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(polyB, bottom_ptsB, G_glob);
    polyB.bottom_side_new = polyB.bottom_side;
    polyB.top_side_new = polyB.top_side;
    polyB.bottom_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(polyB.bottom_side_new, bottom_nodesB, closest_mask, 'linear');
    bottom_nodesB(rem_idx,2) = z_new;

    polyB = PinchOuts.coordCorrectionSubgridBC(polyB, top_nodesB, bottom_nodesB, east_nodesAB);    
    %polyB = interpolateHorizontal(polyB);

    polyB = logicalIndicesUpperOverlap(polyB, poly_upper);      

    p_idxB = strcat(p_idx, 'B');
    polyB.p_idx = p_idxB;
    poly_obj.(p_idxB) = polyB;
    nodes_overlap.(p_idxB) = top_nodesB; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxB) = top_ptsB;


    % --- C: Make cartesian subgrid of pXC. ---
    bottom_nodesC = bottom_nodes(x_bottom >= sep_point_bottom(:,1), :);
    num_x_pXC = num_x_pXB; % enforce same number of columns (not strictly necessary, but ensures same number of nodes for bottom and top side of full polygon subgrid
    NzC = z_sep_idx - 1; % -1 to remove end-node

    bottom_ptsC = [sep_point_bottom; bottom_pts(bottom_pts(:,1) >= sep_point_bottom(:,1), :)]; % polygon data points on top side RIGHT of separation point, including top separation point   
    side_ptsC = [top_ptsC; bottom_ptsC];

    polyC = PinchOuts(all_polys, poly_num, side_ptsC);            
    polyC.bottom_side = bottom_ptsC;
    polyC.top_side = top_ptsC;
    polyC = cartesianSubgrid(polyC, Lx_glob, Lz_glob, Nx_glob, Nz_glob, num_x_pXC, NzC);

    [top_nodesC, closest_mask, xs_new] = PinchOuts.separationPoint_pinchout(polyC, top_ptsC, G_glob);
    polyC.bottom_side_new = polyC.bottom_side;
    polyC.top_side_new = polyC.top_side;
    polyC.top_side_new(:,1) = xs_new;

    [z_new, rem_idx] = PinchOuts.interpolateZSide(polyC.top_side_new, top_nodesC, closest_mask, 'linear');
    top_nodesC(rem_idx,2) = z_new;

    polyC = PinchOuts.coordCorrectionSubgridBC(polyC, top_nodesC, bottom_nodesC, east_nodesAC);    
    %polyC = interpolateHorizontal(polyC);

    polyC = logicalIndicesWestOverlap(polyC, polyA);           

    p_idxC = strcat(p_idx, 'C');
    polyC.p_idx = p_idxC;
    poly_obj.(p_idxC) = polyC;
    nodes_overlap.(p_idxC) = top_nodesC; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxC) = top_ptsC;
end
