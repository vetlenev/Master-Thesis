function [G, cellw] = getMesh_bo(mesh, fig_mesh)
%
%
%

switch mesh.type
%     case 'coarse'
%         meshpath = 'mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/nodes_coordinates_full.dat';
%         %cellw = [1201 1289];
%         cellw = []; %5722; %355;  %5716    
%     case 'coarseUpper'  
%         meshpath = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/nodes_coordinates_full.dat');
%         ucids = load(fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids_2D_coarse.mat'));
%         cellw = []; 
%     case 'fine'
%         meshpath = '';
%         cellw = [];
%     case 'ref300'
%         meshpath = 'mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref300/nodes_coordinates_full.dat';
%         cellw = [];
%     case 'ref200'
%         meshpath = 'mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref200/nodes_coordinates_full.dat';
%         cellw = [];
%     case 'ref100'
%         meshpath = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref100/nodes_coordinates_full.dat';
%         cellw = [];
%     case 'ref50'
%         meshpath = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref50/nodes_coordinates_full.dat';
%         cellw = [];
%     case 'ref25'
%         meshpath = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref25/nodes_coordinates_full.dat';
%         ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/ucids_ref25_2D.mat');
%         cellw = [];  
%     case 'ref12.5'
%         meshpath = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref12.5/nodes_coordinates_full.dat';
%         %ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref12.5/ucids_ref12.5_2D.mat');
%         cellw = [];
%         
%     %  second scenario
%     case 'sc21'
%         meshpath = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/initial/nodes_coordinates_full.dat';
%         ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/init/ucids_sc21_2D.mat');
%         cellw = [];
        
    case 'sc2'
        meshpath = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/throwDiv/nodes_coordinates.dat');
        ucids = load(fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/ucids_sc2_2D.mat'));
        cellw = [];     
end

vertices = dlmread(meshpath);
x = vertices(:,1); y = vertices(:, 2).*-1;

% Delaunay triangulation
%thick = [linspace(80,80,50)];
if isfield(mesh, 'thick')
    thick = mesh.thick;
else
    warning('Check extruded dimension. Correct number of layers?')
    thick = [logspace(3.4, 2, 25) 71 logspace(2, 3.4, 25)];
end
G = triangleGrid([x, y]);

if mesh.reduce == 1
    if strcmp(mesh.type, 'sc2')
        G = removeCells(G, [ucids.unit_cell_ids{56:58}]); % old, resLM1, marga
    elseif strcmp(mesh.type, 'sc21')
        G = removeCells(G, [ucids.unit_cell_ids{51:53}]); % old, resLM1, marga
    else
        G = removeCells(G, [ucids.unit_cell_ids{1:3}]);   % old, resLM1, marga
    end
elseif mesh.reduce == 2
    G = removeCells(G, [ucids.unit_cell_ids{[1:3 6:7 11:13]}]); % old, resLM1, marga, MMUM, you and associated fault domains
end

% Extrude grid
G = makeLayeredGrid(G, thick);

% Permute axes
G.nodes.coords = G.nodes.coords(:,[3 1 2]);
%G.cells.centroids = G.cells.centroids(:,[3 1 2]);

% Compute grid geometry
G = computeGeometry(G);

% Gravity
gravity reset on; %g = norm(gravity);

% plot Grid
if fig_mesh == 1
    %axisarg = {'FontSize', 11, 'TickLabelInterpreter','latex'};
    %latx={'Interpreter','latex'};
    %fontsz = {'fontsize',12};
    %tit_arg = {'fontsize',13,'Interpreter','latex'};
    
    %h = figure(1);
    %plotGrid(G, 'FaceColor', 'k', 'FaceColor', [255 255 204]./255, 'FaceAlpha', 0.2, 'EdgeColor', [0.5 0.5 0.5]);
    plotGrid(G)
    %title('\textbf{Grid}', tit_arg{:});
    %set(gca, axisarg{:}); xlabel('x [m]', fontsz{:}, latx{:}); ylabel('y [m]', fontsz{:}, latx{:});
    %zlabel('z [m]', fontsz{:}, latx{:})
    %set(gca,'Zdir','reverse') % z coord (here depth) positive
    %set(gca,'Xdir','reverse')
    set(gca,'Ydir','reverse')
    axis equal
    %xlim([0 sum(thick)]); ylim([0 45000]); zlim([0 max(G.faces.centroids(:,G.griddim))]);
    %set(h, 'Position', [50, 50, 800, 200])
    view([-65,10])
end 

return
end