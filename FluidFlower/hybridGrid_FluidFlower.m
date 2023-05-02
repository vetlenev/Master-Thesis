%% Load modules and geometry
mrstModule add matlab_bgl coarsegrid

filename = 'draft_spe11_with_facies_markers.geo';
[pts, loops, facies] = parse_spe11_geo(filename);
% pts (287 x 3): 3D coords of each Point
% loops (32 x variable): ordered array of points forming the surface of loop i
% facies (7 x variable): cell array of surfaces (set of loops) connected to each facie 

%% Directory and seed
rootdir = strrep(ROOTDIR, '\', '/');
make_figs_thesis = false;
make_pinch_fig = false;
make_fig_bottomfault = false;
% NB: change this to your local path to FluidFlower folder
simtype_dir = 'hybrid/extended/';
simtype_dir_fine = 'fine/';
output_filename = 'sealing_1_2_3_regions_muchCO2';
data_dir = strcat(rootdir, '../Master-Thesis/FluidFlower/data');
plot_dir = strcat(rootdir, '../Master-Thesis/FluidFlower/output/', simtype_dir, output_filename, '/');
plot_dir_fine = strcat(rootdir, '../Master-Thesis/FluidFlower/output/', simtype_dir_fine, output_filename, '/');
mkdir(data_dir);
mkdir(plot_dir);
mkdir(plot_dir_fine);

my_seed = 295;
seed = UtilFunctionsFF.setSeed(data_dir, my_seed);
rng(seed)

%% SET PARAMETERS
node_density = 0.5; % The density of nodes inside triangulated parts of faults. Should be < 1 to not get to dense triangles (standard: 0.3)
nx_glob = 200; % 200! % standard: nx_glob=250, nz_glob=150, works for this setting...
nz_glob = 150;

%% Other definitions
poly = struct; % to hold each polynomial
pts_overlap = struct; % to hold overlapping polygonal points (provided from data file)
nodes_overlap = struct; % to hold overlapping nodes from generated subgrids

logicalBottom = @(G) G.nodes.coords(:,2) == min(G.nodes.coords(:,2));
logicalTop = @(G) G.nodes.coords(:,2) == max(G.nodes.coords(:,2));

gravity reset on

%% Global grid
[G_glob, x_glob, z_glob] = PolygonGrid.globalCartGrid(pts, nx_glob, nz_glob);

Lx = max(G_glob.faces.centroids(:,1));
Lz = max(G_glob.faces.centroids(:,2));
N = G_glob.cartDims;
Nx = N(1); Nz = N(2);

%% Organize polygons
polys = PolygonGrid.polyPoints(pts, loops, facies);
% Upper, fully horizontally extending surface: 32
poly.top = PolygonGrid(polys, 32);

ptop_orig = poly.top.p_orig;
ptop = poly.top.p;
poly.top.p_idx = 'p32';

%% Find top and bottom surfaces
[poly.top, split_pts] = reorderPts(poly.top, 32);

[top_side, bottom_side] = topAndBottomSurfaces(poly.top, split_pts, 32, pts_overlap);
poly.top.top_side = top_side;
poly.top.bottom_side = bottom_side;

poly.top = cartesianSubgrid(poly.top, Lx, Lz, Nx, Nz, []); % []: no conforming of horizontal grid size

%% dx-correction
% For each surface point, find closest x-coord in subgrid and change it to
% equal the surface coordinate
poly.top.bottom_mask = logicalBottom(poly.top.G); % only select nodes at bottom of grid
poly.top.top_mask = logicalTop(poly.top.G);
poly.top.G = computeGeometry(poly.top.G);
% Bottom:
[poly.top, closest_bottom] = correct_dx_poly(poly.top, G_glob, poly.top.bottom_mask, 'bottom');
% Top:
[poly.top, closest_top] = correct_dx_poly(poly.top, G_glob, poly.top.top_mask, 'top');

if make_figs_thesis
    fig_top = figure();
    plotGrid(poly.top.G, 'facecolor', 'blue');
    hold on
    scatter(poly.top.bottom_side(:,1), poly.top.bottom_side(:,2), 60, 'MarkerFaceColor', 'yellow', 'LineWidth', 1, 'MarkerEdgeColor', 'black');
    hold on
    scatter(poly.top.bottom_side_new(:,1), poly.top.bottom_side_new(:,2), 10, 'MarkerFaceColor', 'red');
    xlim([min(poly.top.bottom_side(:,1)), max(poly.top.bottom_side(:,1))])
    ylim([1, 1.2])
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_0'));
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_0'), 'pdf');
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_0'), 'png');
end

%% Interpolation
% Interpolate z-values at surface of polygon
poly.top = interpolateZ_remaining(poly.top, closest_bottom, ...
                                poly.top.bottom_mask, 'bottom', 'linear'); % top face is the "remaining" parts of the polygon we want bottom surface
poly.top = interpolateZ_remaining(poly.top, closest_top, ...
                                poly.top.top_mask, 'top', 'linear');

if make_figs_thesis
    fig_top = figure();
    plotGrid(poly.top.G, 'facecolor', 'blue', 'facealpha', 0.3);
    hold on
    scatter(poly.top.G.nodes.coords(poly.top.top_mask, 1), poly.top.G.nodes.coords(poly.top.top_mask, 2), 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none') %[50, 205, 50]/255
    hold on
    scatter(poly.top.G.nodes.coords(poly.top.bottom_mask, 1), poly.top.G.nodes.coords(poly.top.bottom_mask, 2), 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none')
    hold on
    scatter(poly.top.bottom_side(:,1), poly.top.bottom_side(:,2), 60, 'MarkerFaceColor', 'yellow', 'LineWidth', 1, 'MarkerEdgeColor', 'black');
    hold on
    scatter(poly.top.bottom_side_new(:,1), poly.top.bottom_side_new(:,2), 10, 'MarkerFaceColor', 'red');  
    xlim([min(poly.top.bottom_side(:,1)), max(poly.top.bottom_side(:,1))])
    ylim([1, 1.2])
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_a'));
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_a'), 'pdf');
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_a'), 'png');
end

% Interpolate x+z values in interior
poly.top = interpolateInternal(poly.top, poly.top.top_mask, poly.top.bottom_mask, []);

if make_figs_thesis
    fig_top = figure();
    plotGrid(poly.top.G, 'facecolor', 'blue', 'facealpha', 0.3, 'edgealpha', 0.2);
    hold on
    %scatter(poly.top.G.nodes.coords(:, 1), poly.top.G.nodes.coords(:, 2), 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none')
    %hold on    
    %scatter(poly.top.bottom_side(:,1), poly.top.bottom_side(:,2), 60, 'MarkerFaceColor', 'yellow', 'LineWidth', 1, 'MarkerEdgeColor', 'black');
    %hold on
    scatter(poly.top.bottom_side_new(:,1), poly.top.bottom_side_new(:,2), 40, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor','k');  
    xlim([min(poly.top.bottom_side(:,1)), max(poly.top.bottom_side(:,1))])
    ylim([1, 1.2])
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_b'));
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_b'), 'pdf');
    saveas(fig_top, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_b'), 'png');
end

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

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 28, 32, ...
                                                        nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap, pts_overlap] = gluePinchOuts(polys, poly, 29, 32, 25, ...
                                                     nodes_overlap, pts_overlap, G_glob);

if make_figs_thesis
    fig_top_glued = figure();
    plotGrid(poly.p32.G, 'facecolor', 'blue', 'facealpha', 0.3, 'edgealpha', 0.1);
    plotGrid(poly.p30.G, 'facecolor', 'red', 'facealpha', 0.3, 'edgealpha', 0.1);
    plotGrid(poly.p28.G, 'facecolor', 'red', 'facealpha', 0.3, 'edgealpha', 0.1);
    plotGrid(poly.p29A.G, 'facecolor', 'red', 'facealpha', 0.3, 'edgealpha', 0.1);
    plotGrid(poly.p29B.G, 'facecolor', 'red', 'facealpha', 0.3, 'edgealpha', 0.1);
    plotGrid(poly.p29C.G, 'facecolor', 'red', 'facealpha', 0.3, 'edgealpha', 0.1);
    hold on
    %scatter(poly.p30.top_side(:,1), poly.p30.top_side(:,2), 50, 'MarkerFaceColor','yellow', 'MarkerEdgeColor','black');
    %hold on
    scatter(nodes_overlap.p30(:,1), nodes_overlap.p30(:,2), 10, 'MarkerFaceColor','red', 'MarkerEdgeColor','none');
    hold on
    scatter(nodes_overlap.p28(:,1), nodes_overlap.p28(:,2), 10, 'MarkerFaceColor','red', 'MarkerEdgeColor','none');
    hold on
    scatter(nodes_overlap.p29A(:,1), nodes_overlap.p29A(:,2), 10, 'MarkerFaceColor','red', 'MarkerEdgeColor','none');
    hold on
    scatter(nodes_overlap.p29B(:,1), nodes_overlap.p29B(:,2), 10, 'MarkerFaceColor','red', 'MarkerEdgeColor','none');
    %xlim([0, max(poly.p32.G.nodes.coords(:,1))])
    xlim([0.7, 1.8])
    %ylim([min(poly.p30.G.nodes.coords(:,2)), 1.1])
    ylim([0.94, 1.12])
    daspect([1 0.25 1])
    saveas(fig_top_glued, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_c'));
    saveas(fig_top_glued, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_c'), 'pdf');
    saveas(fig_top_glued, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_top_c'), 'png');
end

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 18, 30, ...
                                                        nodes_overlap, pts_overlap, G_glob, true, 1.5);

[poly, nodes_overlap, pts_overlap] = glueToUpperPolygon(polys, poly, 13, 18, ...
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
[poly, nodes_overlap, pts_overlap, ...
    top_nodes, bottom_nodes, west_nodes] = gluePolygon1(polys, poly, 1, [3,2], ...
                                                   nodes_overlap, pts_overlap, G_glob, node_density);

colors = [[0 0.4470 0.7410]; 
          [0.8500 0.3250 0.0980];
          [0.9290 0.6940 0.1250];
          [0.4940 0.1840 0.5560];
          [0.4660 0.6740 0.1880];
          [0.3010 0.7450 0.9330];
          [0.6350 0.0780 0.1840]];

%% Triangulation of Top Left Fault
poly31_neighbors = cell(11, 2);
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

% Quadrilateration
poly = Faults.makeQuadrilaterals(polys, poly, 31, poly31_neighbors, ...
                                    2, Lx, Lz, Nx, Nz);

[poly, p31At] = Faults.QuasiRandomPointDistribution(poly, "p31fAt", G_glob, node_density, false, colors);

%% Triangulation of Top Right Fault
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

% Quadrilateration
poly = Faults.makeQuadrilaterals(polys, poly, 25, poly25_neighbors, ...
                                    2, Lx, Lz, Nx, Nz);
% Triangulate bottom part of fault and end-point of pinch-outs
[poly, p25E] = Faults.QuasiRandomPointDistribution(poly, "p25fE", G_glob, node_density, false, colors);
[poly, p25Ft] = Faults.QuasiRandomPointDistribution(poly, "p25fFt", G_glob, 2*node_density, false, colors);
[poly, p25Gt] = Faults.QuasiRandomPointDistribution(poly, "p25fGt", G_glob, 2*node_density, false, colors);

%% Save pinch-out fig
if make_pinch_fig
    fig_pinch = figure();
    plotGrid(poly.p29A.G, 'facecolor', 'red', 'facealpha', 0.7)
    plotGrid(poly.p29B.G, 'facecolor', 'red', 'facealpha', 0.5)
    plotGrid(poly.p29C.G, 'facecolor', 'red', 'facealpha', 0.3)    
    plotGrid(poly.p25fFq.G, 'facecolor', [0, 191, 255]/255)
    plotGrid(poly.p25fFt.G, 'facecolor', [0, 191, 255]/255)
    hold on
    scatter(nodes_overlap.p29C(:,1), nodes_overlap.p29C(:,2), 23, 'MarkerFaceColor', 'yellow')
    hold on
    scatter(poly.p29B.G.nodes.coords(poly.p29B.bottom_mask, 1), poly.p29B.G.nodes.coords(poly.p29B.bottom_mask, 2), 23, 'MarkerFaceColor', 'yellow')
end

%% Triangulation of Right Fault
poly12_neighbors = cell(2, 2);
poly12_neighbors{1,1} = poly.p5B; poly12_neighbors{1,2} = 'top';
poly12_neighbors{2,1} = poly.p5C; poly12_neighbors{2,2} = 'bottom';

poly = Faults.makeQuadrilaterals(polys, poly, 12, poly12_neighbors, ...
                                    1.5, Lx, Lz, Nx, Nz);

[poly, p12At] = Faults.QuasiRandomPointDistribution(poly, "p12fAt", G_glob, node_density, false, colors);

%% Triangulation of Bottom Fault
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

[poly, nodes_overlap] = triangulateBottomFault(polys, poly, polyBF_nums, polyBF_norder, ...
                                                   nodes_overlap, pts_overlap, G_glob, 2);

%% Triangulation with added internal points (for internal faults)
[poly, nodes_overlap] = triangulateInternalFaults(poly, 24, ...
                                                 nodes_overlap, pts_overlap, G_glob);
   
[poly, nodes_overlap] = triangulateInternalFaults(poly, 9, ...
                                                  nodes_overlap, pts_overlap, G_glob);
 
[poly, nodes_overlap] = triangulateInternalFaults(poly, 21, ...
                                                   nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap] = triangulateInternalFaults(poly, 26, ...
                                                    nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap] = triangulateInternalFaults(poly, 22, ...
                                                    nodes_overlap, pts_overlap, G_glob);

[poly, nodes_overlap] = triangulateInternalFaults(poly, 23, ...
                                                    nodes_overlap, pts_overlap, G_glob);

%% Triangulate internal parts separately
%bf_polys = [24,9,21,26,22,23]; % polygons constituting bottom fault
bf_polys = [26,22]; % other polygons are cartesian quadrilaterals

% Triangulate irregular twin polygons
[poly, p26_22] = Faults.QuasiRandomPointDistribution(poly, bf_polys, G_glob, node_density, false, colors);
% Quadrilaterate the more regular polygons
poly = Faults.heterogeneousFault(polys, poly, 1.5, Lx, Lz, Nx, Nz);
% Triangulate bottom-most part
[poly, p23B] = Faults.QuasiRandomPointDistribution(poly, "p23fB", G_glob, node_density, false, colors);

%% Plot part of bottom fault
if make_fig_bottomfault
    fig_bfault = figure();
    plotGrid(poly.p23fA.G, 'facecolor', [1, 191, 255]/255)
    plotGrid(poly.p26f.G, 'facecolor', 'red')
    plotGrid(poly.p22f.G, 'facecolor', [128, 0, 128]/255)
    hold on
    %scatter(poly.p26f.internal_east(:,1), poly.p26f.internal_east(:,2), 25, 'MarkerFaceColor','yellow','MarkerEdgeColor','black')
    %hold on
    w23 = poly.p23fA.G.nodes.coords(poly.p23fA.west_mask, :);
    [uwe23, ~, ic] = uniquetol(we, 'ByRows', true);
    counts = accumarray(ic, 1);
    collapsed_nodes = unique(we(counts(ic) == 2, :), 'rows');
    scatter(collapsed_nodes(:,1), collapsed_nodes(:,2), 25, 'MarkerFaceColor','yellow','MarkerEdgeColor','black')
    xlim([min(poly.p23fA.G.nodes.coords(:,1)), max(poly.p22f.G.nodes.coords(:,1))])  
    ylim([min(poly.p23fA.G.nodes.coords(:,2)), max(poly.p26f.G.nodes.coords(:,2))])
    saveas(fig_bfault, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_bottomfault'));
    saveas(fig_bfault, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_bottomfault'), 'pdf');
    saveas(fig_bfault, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_bottomfault'), 'png');
end

%% Plots / tests
fig_grid = UtilFunctions.fullsizeFig(2);
xlim([0, 2.8])
ylim([0, 1.2])
%colormap('jet')

poly_idxs = fieldnames(poly);
overlap_idxs = fieldnames(nodes_overlap);
for i = 1:numel(poly_idxs) % plot subgrids first
    px = poly_idxs{i};
    if ~isempty(poly.(px).G) && isempty(regexp(px, 'p\w*PEBI', 'once'))
        plotGrid(poly.(px).G, 'EdgeAlpha', 0.3, 'FaceColor', colors(poly.(px).facies,:))
        %plotCellData(poly.(px).G, repmat(poly.(px).facies, poly.(px).G.cells.num, 1), 'EdgeAlpha', 0.1)        
        hold on
    end    
    for j = 1:numel(overlap_idxs) % plot overlapping nodes
        px = overlap_idxs{j};
        plot(nodes_overlap.(px)(:,1), nodes_overlap.(px)(:,2), 'k.', 'markersize', 6)
        hold on
%         plot(pts_overlap.(px)(:,1), pts_overlap.(px)(:,2), 'b.', 'markersize', 15)
%         hold on
    end
end

saveas(fig_grid, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_fullgrid'));
saveas(fig_grid, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_fullgrid'), 'png');
saveas(fig_grid, strcat(rootdir, '../Master-Thesis/FluidFlower/geometry/figs/FF_fullgrid'), 'pdf');
%plotGrid(poly.p25pA.G, 'EdgeAlpha', 0.3, 'FaceColor', colors(poly.p25pA.facies,:))

%% Glue-preparation
% Remove all PEBI-polygons to check if triangulation causes problems
             
poly_idxs = {poly.p31fAq.p_idx; poly.p31fAt.p_idx; ...
             poly.p25fA.p_idx; poly.p25fB.p_idx; poly.p25fC.p_idx; ...             
             poly.p29A.p_idx; poly.p29B.p_idx; poly.p25fFq.p_idx; ...
             poly.p25fFt.p_idx; poly.p29C.p_idx; poly.p19A.p_idx; ...
             poly.p19B.p_idx; poly.p25fGq.p_idx; poly.p25fGt.p_idx; ...
             poly.p25fD.p_idx; poly.p25fE.p_idx; poly.p28.p_idx; ...
             poly.p20.p_idx; poly.p15.p_idx; poly.p6.p_idx; ...
             poly.p19C.p_idx; poly.p14.p_idx; poly.p4.p_idx; ...
             poly.p30.p_idx; poly.p18.p_idx; poly.p13.p_idx; ...
             poly.p5A.p_idx; poly.p5B.p_idx; poly.p12fAq.p_idx; ...
             poly.p12fAt.p_idx; poly.p5C.p_idx; poly.p27.p_idx; ...
             poly.p17.p_idx; poly.p11left.p_idx; poly.p11mid.p_idx; ...
             poly.p11right.p_idx; poly.p16.p_idx; poly.p10.p_idx; ...
             poly.p8.p_idx; poly.p3.p_idx; poly.p7left.p_idx; ...
             poly.p7mid.p_idx; poly.p7right.p_idx; poly.p7small.p_idx; ...
             poly.p2left.p_idx; poly.p2mid.p_idx; poly.p2right.p_idx; ...
             poly.p24fA.p_idx; poly.p24fB.p_idx; poly.p24fC.p_idx; ...
             poly.p9fA.p_idx; poly.p21fA.p_idx; poly.p21fB.p_idx; ...
             poly.p21fC.p_idx; poly.p21fD.p_idx; poly.p22f.p_idx; ...
             poly.p26f.p_idx; poly.p23fA.p_idx; poly.p23fB.p_idx; ...
             poly.p1left.p_idx; poly.p1mid.p_idx; poly.p1rightA.p_idx; ...
             poly.p1rightB.p_idx; poly.p1rightC.p_idx; poly.p1small.p_idx};
             

%% Glue together subgrids
poly.glued = poly.p32; % top grid is basis for gluing
poly.glued.G.facies = repmat(poly.glued.facies, poly.glued.G.cells.num, 1);
poly.glued.cell_range.p32 = [1, poly.p32.G.cells.num];

for i=1:numel(poly_idxs)    
    px = poly_idxs{i};
    disp(px)
    G_sub = poly.(px).G; 
    G_glued = glue2DGrid_FF(poly.glued.G, G_sub, 'poly_idx', px);    

    new_cells = G_sub.cells.indexMap + poly.glued.G.cells.num;
    poly.glued.cell_range.(px) = [poly.glued.G.cells.num+1, poly.glued.G.cells.num+max(G_sub.cells.indexMap)];

    G_glued.cells.indexMap = [poly.glued.G.cells.indexMap; new_cells];
    G_glued.facies = [poly.glued.G.facies; repmat(poly.(px).facies, G_sub.cells.num, 1)];
    G_glued.i = [poly.glued.G.i; G_sub.i];
    G_glued.j = [poly.glued.G.j; G_sub.j];

    poly.glued.G = G_glued;  
end

poly.glued.G = computeGeometry(poly.glued.G);

figure(4)
plotGrid(poly.glued.G)
plotCellData(poly.glued.G, double(isnan(poly.glued.G.i)))

%% Remove impermeable layers from glued grid
Gg = poly.glued.G;
poly.sim = poly.glued;
[G_sim, gc, gf] = extractSubgrid(poly.glued.G, poly.glued.G.facies ~= 7);
% gc == G_sim.cells.indexMap 
G_sim.i = Gg.i(gc);
G_sim.j = Gg.j(gc);
G_sim.facies = Gg.facies(gc);
poly.sim.G = G_sim;

for i=1:numel(poly_idxs)    
    px = poly_idxs{i};
    G_sub = poly.(px).G; 
    pg = poly.glued;
    % --------------------
    % TODO: Find mapping from glued facies to sim facies!
    % --------------------
    if all(pg.G.facies(pg.cell_range.(px)) == 7) % not included in new grid
        poly.sim.cell_range = rmfield(poly.sim.cell_range, px);
    else
        poly.sim.cell_range.(px) = find(ismember(gc, pg.cell_range.(px)));
    end
end

Gg = poly.sim.G; % THIS IS THE COMPLETE GLUED GRID

%% Find neighboring structured cells to fault
nan_cells = find(isnan(Gg.i)); % all nan logical indices are fine-scale triangle cells
gN = Gg.faces.neighbors;
gN1 = gN(:,1);
gN2 = gN(:,2);

nan_N = ismember(gN, nan_cells);
fault_N = gN(sum(nan_N, 2) == 1, :); % if only one neighbor of a face is nan, this face is at boundary of fault
fault_idx = nan_N(sum(nan_N, 2) == 1, :); % index of nan cell for each neighbor pair

% select cartesian neighbor, not triangle cell
fault_buffer = zeros(size(fault_N,1),1);
for i=1:size(fault_N,1)
    fault_buffer(i) = fault_N(i, ~fault_idx(i,:)); % ~nn will select no-nan neighbor, i.e., the cartesian cell
end
fault_buffer = fault_buffer(fault_buffer ~= 0);

% Add impermeable fault (by extracting cells adjacent to "internal"
% boundary) -> to give correct VE discretization
bcells = Gg.faces.neighbors(boundaryFaces(Gg),:);
bcells_all = bcells(:);
bcells_all = bcells_all(bcells_all ~= 0);
imperm_fault = bcells_all(Gg.cells.centroids(bcells_all,1) > 0.5 & ...
                            Gg.cells.centroids(bcells_all,1) < 1.5 & ...
                            Gg.cells.centroids(bcells_all,2) > 0.5 & ...
                            Gg.cells.centroids(bcells_all,2) < 1.15);

fault_buffer = [fault_buffer; imperm_fault];

fault_buffer = fault_buffer(~ismember(fault_buffer, nan_cells));

%% Some plotting
if false
    figure(3)
    %colormap('jet')
    plotGrid(Gg)    
    colorbar('southoutside')
    %plotCellData(Gg, double(Gg.i == 46))
    ss = zeros(Gg.cells.num,1);
    ss(fault_buffer) = 1;    
    plotCellData(Gg, ss)
    hold on
    plotFaces(Gg, boundaryFaces(Gg), 'EdgeColor', 'red', 'linewidth', 2)
end

%% Read deck -> assign fluids and rock
deck = readEclipseDeck('deck/CSP11A.DATA');
deck = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);

% Change from oil-gas system to water-gas system
[fluid.krO] = fluid.krOG;
fluid = rmfield(fluid, 'krOG');
[fluid.krPts.o] = fluid.krPts.og;
fluid.krPts = rmfield(fluid.krPts, 'og');
[fluid.krW] = fluid.krO;
fluid = rmfield(fluid, 'krO');
[fluid.krPts.w] = fluid.krPts.o;
fluid.krPts = rmfield(fluid.krPts, 'o');
[fluid.pcWG] = fluid.pcOG;
fluid = rmfield(fluid, 'pcOG');

b = 1; % unit formation volume factor
fluid.bW = @(p, varargin) 1 + 0.*p; %b*constantReciprocalFVF(p, varargin{:});
fluid.bG = @(p, varargin) 1 + 0.*p; %b*constantReciprocalFVF(p, varargin{:});
fluid.muW = @(p, varargin) 3e-5*Pascal*second + 0.*p; %3e-5*Pascal*second + 0.*p;
fluid.muG = @(p, varargin) 8e-4*Pascal*second + 0.*p; %8e-4*Pascal*second + 0.*p;
%fluid = assignRelPerm(fluid);

facies_all = deck.REGIONS.SATNUM;

perm = zeros(Gg.cells.num, 1);
poro = zeros(Gg.cells.num, 1);
for i=1:numel(unique(facies_all))
    perm(Gg.facies == i) = unique(deck.GRID.PERMX(facies_all == i));
    poro(Gg.facies == i) = unique(deck.GRID.PORO(facies_all == i));
end
%poro(Gg.facies == 7) = 0.5; % just set to a positive value to avoid errors

rock = makeRock(Gg, perm, poro);
rock.regions.saturation = Gg.facies;

%% Setup model, schedule, state ...
isFine = struct;
isFine.bufferCells = false(Gg.cells.num, 1);
isFine.bufferCells(fault_buffer) = true;
isFine.faultCells = false(Gg.cells.num, 1);
isFine.faultCells(nan_cells) = true;

isFine.faultFaces = false(Gg.faces.num, 1);
for i=find(isFine.faultCells)'
    fault_face = Gg.cells.faces(Gg.cells.facePos(i):Gg.cells.facePos(i+1)-1, 1);
    isFine.faultFaces(fault_face) = true;
end
isFine.bufferFaces = false(Gg.faces.num, 1);
for i=find(isFine.bufferCells)'
    buffer_face = Gg.cells.faces(Gg.cells.facePos(i):Gg.cells.facePos(i+1)-1, 1);
    isFine.bufferFaces(buffer_face) = true;
end

allSealingCells = {}; % cell-indices for sealing cells
allSealingBFaces = {}; % boundary faces of sealing cells
allSealingBottom = {}; % bottom faces of sealing layer (needed for trap analysis)
allLayersBottom = {};
allLayersCells = {};

lowperm_facie = [1,2,3]; % [1,4,7]

%lowperm_cells = find(Gg.facies == lowperm_facie(1) | Gg.facies == lowperm_facie(2));
lowperm_cells = find(ismember(Gg.facies, lowperm_facie));
if ismember(lowperm_facie, 1) % Small p7 part is included in isFine.fault -> remove from sealing cells in structured regions
    p7small_cells = poly.sim.cell_range.p7small(1):poly.sim.cell_range.p7small(2); % poly.glued -> poly.sim
    lowperm_p7small = ismember(lowperm_cells, p7small_cells);
    lowperm_cells(lowperm_p7small) = [];
end
%facie_polys = struct; % append polys part of this facies
lowperm_cells_split = struct; % lowperm cells in separate arrays for each polygon
poly_idxs_new = fieldnames(poly.sim.cell_range);
for i=1:numel(poly_idxs_new)
    pidx = string(poly_idxs_new{i});
    pcells = (poly.sim.cell_range.(pidx)(1):poly.sim.cell_range.(pidx)(2))';

    [sealing_faces, bottom_faces, sealing_cells] = addConfiningLayersFF(Gg, 'cells', pcells);

    if isempty(sealing_faces) || isempty(bottom_faces) || isempty(sealing_cells)
        continue; % don't add sealing layer if triangulated cells
    end
    if all(ismember(pcells, lowperm_cells))
        %facie_polys.(pidx) = poly.(pidx);
        lowperm_cells_split.(pidx) = pcells;
                
        allSealingCells = cat(2, allSealingCells, sealing_cells);
        allSealingBFaces = cat(2, allSealingBFaces, sealing_faces);
        allSealingBottom = cat(2, allSealingBottom, bottom_faces);
    end
    allLayersCells = cat(2, allLayersCells, sealing_cells);
    allLayersBottom = cat(2, allLayersBottom, bottom_faces);    
    %allLayersBottom.(pidx) = bottom_faces;
end

isFine.sealingCells = false(Gg.cells.num, 1);
isFine.sealingCells(vertcat(allSealingCells{:})) = true;

allSealingBFaces = UtilFunctionsFF.removeOverlappingElements(vertcat(allSealingBFaces{:}));
allLayersBottom_rem = UtilFunctionsFF.removeOverlappingElements(vertcat(allLayersBottom{:}));

isFine.sealingBFaces = false(Gg.faces.num, 1);
isFine.sealingBFaces(allSealingBFaces) = true; % bounding faces of sealing cells

isFine.layersBFaces = false(Gg.faces.num, 1);
isFine.layersBFaces(allLayersBottom_rem) = true; % bounding faces of each layer

isFine.sealingBottom = allSealingBottom;

%% Operational schedule
bf = boundaryFaces(Gg);
bfz = Gg.faces.centroids(bf, 2);
top_faces = bf(abs(bfz - max(bfz)) < 1e-10);
bc = addBC([], top_faces, 'pressure', 1e5*Pascal, 'sat', [1,0]);

ii = Gg.i;
jj = Gg.j;
% ds = deck.SCHEDULE;
% dt = ds.step.val;
% t = cumsum(dt);
% tot_time = t(end);
tot_time = 2*day; % 5*day % 18*hour
t_control2 = 5*hour;% 2.5*hour; % 2*hour
t_control3 = 10*hour; % 5*hour; % 4*hour
dt1 = rampupTimesteps(t_control2, t_control2/120, 10); % t_control2/70 % 100
dt2 = rampupTimesteps(t_control3-t_control2, (t_control3-t_control2)/120, 12); % (t_control3-t_control2)/120 % 100
dt3 = rampupTimesteps(tot_time-t_control3, (tot_time-t_control3)/100, 12); % (tot_time-t_control3)/120 % 70
dt = [dt1; dt2; dt3];
t = cumsum(dt);

inj_rate_W1 = 3.5*10^(-7) * kilogram/second; % 1.6*10^(-7) * kilogram/second
inj_rate_W2 = 3.5*10^(-7) * kilogram/second; % 1.6*10^(-7) * kilogram/second

xc = Gg.cells.centroids(:,1);
zc = Gg.cells.centroids(:,2);
min_x = @(x) min(abs(xc - x));
min_z = @(z) min(abs(zc - z));
idx_x = @(x) find(xc - x == min_x(x));
idx_z = @(z) find(zc - z == min_z(z));

dist2pt = @(c,p) sqrt((p(:,1)-c(:,1)).^2 + (p(:,2)-c(:,2)).^2);
[~, W1_cell] = min(dist2pt(Gg.cells.centroids, [0.9,0.3]));
[~, W2_cell] = min(dist2pt(Gg.cells.centroids, [1.7,0.7]));

W1_x = ii(W1_cell);
W1_z = jj(W1_cell);

W = addWell([], Gg, rock, W1_cell, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W1, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1], ...
                'status', true); % starts injecting

W2_x = ii(W2_cell);
W2_z = jj(W2_cell);

W = addWell(W, Gg, rock, W2_cell, ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate_W2, ... % volumetric injection rate
                'radius', 0.001, ...
                'comp_i', [0 1], ...
                'status', false); % starts off

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

schedule.control(2:3) = schedule.control(1);
schedule.control(2).W(2).status = true;
schedule.control(3).W(1).status = false;
schedule.control(3).W(2).status = false;

schedule.step.control(t <= t_control2) = 1;
schedule.step.control(t > t_control2 & t <= t_control3) = 2;
schedule.step.control(t > t_control3) = 3;

nearWell = false(Gg.cells.num, numel(W));
openBC = false(Gg.cells.num, 1);

horzWellDistRate = [0.05, 0.10]; % [0.03, 0.08] % ratio of total horizontal length considered "close" to well
vertWellDistRate = [0.08, 0.12]; % [0.06, 0.12]
for i = 1:numel(W)
    c = W(i).cells(1);
    hdist_i = abs(ii - ii(c)); 
    hdist_x = abs(xc - xc(c));
    vdist_j = abs(jj - jj(c)); 
    vdist_z = abs(zc - zc(c));
    % store cells close to well, if fine-scale needed
    nearWell(hdist_i < horzWellDistRate(1)*max(ii) & ... % max(ii)
             hdist_x < horzWellDistRate(2)*max(xc) & ...             
             vdist_j < vertWellDistRate(1)*max(jj) & ...
             vdist_z < vertWellDistRate(2)*max(zc), i) = true; % max(kk)
end
if ~isempty(bc)
    openBC(sum(Gg.faces.neighbors(bc.face, :), 2)) = true;
end

isFine.well = nearWell;
isFine.bc = openBC;

remFineCells = any(isFine.well,2) | isFine.bc | ...
                isFine.faultCells | isFine.bufferCells;
sealingCells = false(Gg.cells.num, 1);
sealingCells(isFine.sealingCells) = true;

isFine.allFineCells = sealingCells | remFineCells;

isFine.bcFaces = false(Gg.faces.num, 1);
bcFacesAll = addConfiningLayersFF(Gg,'cells',isFine.bc);
isFine.bcFaces(bcFacesAll) = true;

isFine.wellFaces = false(Gg.faces.num, numel(W));
for i=1:numel(W)
    wellFacesAll = addConfiningLayersFF(Gg,'cells',isFine.well(:,i));
    isFine.wellFaces(wellFacesAll, i) = true;
end

z = max(Gg.faces.centroids(:,2)) - Gg.cells.centroids(:,2);
state0 = initResSol(Gg, 1.013e5*Pascal + fluid.rhoWS*norm(gravity)*z, [1,0]); % atmospheric pressure

%% VE-to-VE vertical transition between facies
external_faces = find(all(gN == 0, 2));
gN_int = find(all(gN ~= 0, 2)); % interior neighbors
gN1_int = gN(gN_int, 1);
gN2_int = gN(gN_int, 2);

facies_transition = zeros(Gg.faces.num, 1);
facies_transition(gN_int) = Gg.facies(gN1_int) ~= Gg.facies(gN2_int);

faciesFaces = find(facies_transition);

%% Erase VE-connection if entry pressure exceeded
p_entry = zeros(numel(unique(Gg.facies)), 1);
p_entry_cells = zeros(Gg.cells.num, 1);
for fac=1:numel(p_entry)
    p_entry(fac) = fluid.pcWG{fac}(1e-5);
    fac_cells = Gg.facies == fac;
    p_entry_cells(fac_cells) = p_entry(fac);
end

facies_N = Gg.faces.neighbors(faciesFaces, :);
facies_centroid = reshape(Gg.cells.centroids(facies_N, 2), [], 2); % reshaoe to get two columns (one for each neighbor)

for i=1:size(facies_centroid, 1)
    [~, row_max] = max(facies_centroid(i,:));
    facie_upper = facies_N(i, row_max);
    row_min = 3 - row_max;
    facie_lower = facies_N(i, row_min);
    % CHECK INSTEAD IF P_ENTRY IS HIGHER THAN P_ENTRY IN LAYER WHERE CLOSEST INJECTION OCCURS!
    W_cells = [W1_cell, W2_cell];
    W_centroids = Gg.cells.centroids(W_cells, 1); % closest in horizontal distance

    [~, closest_idx] = min(abs(Gg.cells.centroids(facie_upper, 1) - W_centroids));
    closest_well = W_cells(closest_idx);
    above_well = Gg.cells.centroids(facie_upper, 2) > Gg.cells.centroids(closest_well, 2);
    %isVE = isnan(Gg.i(facie_upper)) && isnan(Gg.i(facie_lower));    
    in_middle = Gg.cells.centroids(facie_upper, 2) > min(poly.p31fAt.bnodes(:,2)) && ...
                Gg.cells.centroids(facie_upper, 2) < max(poly.p31fAq.internal_top(:,2)) && ...
                Gg.cells.centroids(facie_upper, 1) > min(poly.p31fAq.internal_east(:,1)) && ...
                Gg.cells.centroids(facie_upper, 1) < max(poly.p25fA.internal_west(:,1));

    if p_entry_cells(facie_upper) <= p_entry_cells(facie_lower) || ...
            (p_entry_cells(facie_upper) <= p_entry_cells(closest_well) && ...
            above_well) || (Gg.facies(facie_upper) == 4 && Gg.facies(facie_lower) == 5)
        faciesFaces(i) = 0;
    end
end

faciesFaces = faciesFaces(faciesFaces ~= 0); % remove unwanted connections
isFine.faciesFaces = false(Gg.faces.num, 1);
isFine.faciesFaces(faciesFaces) = true;

%% Convert to hybrid model
fluid.rhoWS = 500*kilogram/meter^3;
%model = GenericBlackOilModel(Gg, rock, fluid);
model = TwoPhaseFluidFlowerModel(Gg, rock, fluid, 1, 1, 'useCNVConvergence', true);
model.oil = false;
model = validateModel(model);

remFineFaces = find(any(isFine.wellFaces,2) | isFine.bcFaces);
faultFaces = find(isFine.faultFaces | isFine.bufferFaces);
otherFineFaces = unique([remFineFaces; faultFaces]);
[model_hybrid, model_coarse] = convertToHybridModel_FF(model, isFine.allFineCells, ...
                                                        'sealingFaces', faciesFaces, ...
                                                        'otherFineFaces', otherFineFaces, ...                                                    
                                                        'sealingCells', sealingCells, ...
                                                        'sealingBFaces', isFine.sealingBFaces, ...
                                                        'setSubColumnsFine', false, ...
                                                        'multiplier', 1, ...
                                                        'sumTrans', true, ...
                                                        'p_entry', p_entry);
Gh = model_hybrid.G;
Gf = model.G;


%% Plot discretized grid
bc_faces = find(any(Gg.faces.neighbors == 0, 2));
fafa = [find(isFine.faciesFaces); find(isFine.sealingBFaces); ...
        find(isFine.bufferFaces); find(any(isFine.wellFaces,2)); ...
        bc_faces; external_faces];

f4 = figure(4);
clf(4)
disc = model_hybrid.G.cells.discretization;
ph = model_hybrid.G.partition;
%plotGrid(Gg, find(disc(ph) == 1798), 'edgecolor', 'yellow')
plotCellData(Gg, disc(ph),'edgealpha', 0)
plotOutlinedGrid(Gg, W, bc, isFine.sealingBFaces | any(isFine.wellFaces,2) | isFine.faultFaces | isFine.bufferFaces | isFine.bcFaces);
hold on
plotFaces(Gf, fafa, 'edgecolor', 'yellow')
cmap = colorcube(12500);
%cmap= 'jet';
colormap(cmap)
xlabel('[m]')
ylabel('[m]')
xlim([0,2.8])
ylim([0,1.2])

%saveas(f4, strcat(plot_dir, 'hybrid_partition'), 'png');

%% Upscale settings
schedule_hybrid = upscaleSchedule(model_hybrid, schedule);
state0_hybrid = upscaleState(model_hybrid, model, state0);
nls = NonLinearSolver('maxIterations', 70);
nls.useRelaxation = true;
nls.maxTimestepCuts = 8;

%% Simulate hybrid!
%[ws_hybrid, states_hybrid] = simulateScheduleAD(state0_hybrid, model_hybrid, schedule_hybrid, 'NonLinearSolver', nls);
problem_hybrid = packSimulationProblem(state0_hybrid, model_hybrid, schedule_hybrid, ...
                                    strcat('FF/', simtype_dir), 'NonLinearSolver', nls, 'Name', output_filename);
% warn = warning('query', 'last');
% w_id = warn.identifier;
% warning('off', w_id);
[ok, status] = simulatePackedProblem(problem_hybrid);

[ws_hybrid, states_hybrid, report_hybrid] = getPackedSimulatorOutput(problem_hybrid);

%% Simulate fine!
%[ws, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
problem = packSimulationProblem(state0, model, schedule, ...
                                    strcat('FF/', simtype_dir), 'NonLinearSolver', nls, 'Name', output_filename);

[ok, status] = simulatePackedProblem(problem);

[ws, states, report] = getPackedSimulatorOutput(problem);

%% Reconstruct fine state
states_hybrid_fs = convertHybridStates_FF(model_hybrid, model, states_hybrid);

%% Plot settings
c1 = [255, 255, 255]/255;
c2 = [48, 37, 255]/255;
c3 = [0, 255, 0]/255;
cc = interp1([0; 1e-5; 1], [c1; c2; c3], (0:0.001:1)');

t = cumsum(schedule.step.val)/hour();
t1 = find(schedule.step.control == 1);
t1 = t1(end);
t2 = find(schedule.step.control == 2);
t2 = t2(end);
t3 = find(schedule.step.control == 3);
t3 = t3(end);

swr = fluid.krPts.w(3,2);

%% Plot hybrid
fafap = [find(isFine.faciesFaces); find(isFine.sealingBFaces); ...
        find(isFine.bufferFaces)];
W1_faces = Gg.cells.faces(Gg.cells.facePos(W1_cell):Gg.cells.facePos(W1_cell+1), 1);
W2_faces = Gg.cells.faces(Gg.cells.facePos(W2_cell):Gg.cells.facePos(W2_cell+1), 1);

for i = [t1, t2, t3]
    f5 = UtilFunctionsFF.fullsizeFig(5); clf(5)     
    plotCellData(Gh, states_hybrid{i}.s(:,2), 'edgec', 'none');
    plotFaces(Gh, (1:Gh.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on    
    colormap(cc); 
    %caxis([0, 1-swr]); 
    caxis([0, max(states_hybrid{end}.sGmax)])
    colorbar('location','southoutside');
    axis tight off
    title({'Hybrid saturation', sprintf('hour: %.1f',t(i))})       
    %saveas(f5, sprintf(strcat(plot_dir, 'coarse_sat_%d'), i), 'png');       
end

for i = [t1, t2, t3]
    f6 = UtilFunctionsFF.fullsizeFig(6); clf(6)
    %f1 = figure(1); clf
    plotCellData(Gf, states_hybrid_fs{i}.s(:,2), 'edgec', 'none');
    plotFaces(Gf, (1:Gf.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(Gf, fafap, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.4)
    plotFaces(Gf, [W1_faces; W2_faces], 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 2)
    colormap(cc); 
    %caxis([0, 0.1]); 
    caxis([0, max(states_hybrid_fs{end}.sGmax)]);
    colorbar('location','southoutside');
    axis tight off
    title({'Reconstructed fine-scale saturation', sprintf('hour: %.1f',t(i))})   
    %setDaspect(run3D, standard);
    %saveas(f6, sprintf(strcat(plot_dir, 'reconstructed_sat_%d'), i), 'png');       
end

%% Fine-scale solution
for i = [t1, t2, t3] %1:20:numel(states)
    f7 = UtilFunctionsFF.fullsizeFig(7); clf(7)     
    plotCellData(Gf, states{i}.s(:,2), 'edgec', 'none');
    plotFaces(Gf, (1:Gf.faces.num)', 'edgec', 'k', 'facealpha', 0, 'edgealpha', 0.1, 'linewidth', 0.1)
    hold on
    plotFaces(Gf, fafap, 'edgec', 'k', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 0.4)
    plotFaces(Gf, [W1_faces; W2_faces], 'edgec', 'red', 'facealpha', 0, 'edgealpha', 1, 'linewidth', 2)
    colormap(cc); 
    % caxis([0, 1-swr]); 
    caxis([0, max(states{end}.sMax(:,2))])
    colorbar('location','southoutside');
    axis tight off
    title({'Full-dimensional saturation', sprintf('hour: %.1f',t(i))})       
    %saveas(f7, sprintf(strcat(plot_dir_fine, 'fine_sat_%d'), i), 'png');       
end



%% FUNTCION DEFINTIONS
function [poly_obj, nodes_overlap, pts_overlap] = glueToUpperPolygon(all_polys, poly_obj, poly_num, poly_num_upper, ...
                                                                        nodes_overlap, pts_overlap, G_glob, varargin)   
    fac_scale = 1;
    if nargin == 8
        inter_horz = varargin{1};
    elseif nargin > 8
        inter_horz = varargin{1};
        fac_scale = varargin{2};
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
    
    poly = cartesianSubgrid(poly, Lx, Lz, Nx, Nz, num_x_overlap, [], fac_scale);                
    
    poly.G.nodes.coords(poly.top_mask, :) = nodes_overlap.(p_idx); % glue to overlapping nodes on top side
    
    [poly, closest_bottom] = correct_dx_poly(poly, G_glob, ...
                                                poly.bottom_mask, 'bottom'); % manually define nodes on bottom side
       
    poly = interpolateZ_remaining(poly, closest_bottom, ...
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
    dist_pinch = @(node) sqrt((node(:,1)-x_pinch).^2 + (node(:,2)-z_pinch).^2);
    top2pinch = dist_pinch(top_nodes);
    bottom2pinch = dist_pinch(bottom_nodes);
    net_dist = top2pinch + bottom2pinch;
    [~, pin_idx] = min(net_dist);
  
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

    polyC = logicalIndicesWestOverlap(polyC, polyA);           

    p_idxC = strcat(p_idx, 'C');
    polyC.p_idx = p_idxC;
    poly_obj.(p_idxC) = polyC;
    nodes_overlap.(p_idxC) = top_nodesC; % choose top nodes since polygon is overlapping with upper polygon
    pts_overlap.(p_idxC) = top_ptsC;
end