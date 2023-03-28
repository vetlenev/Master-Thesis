L_global = 100 * meter;
H1 = 10 * meter;
H2 = 8 * meter;
Xres_global = 60;
Yres1 = 10;
Yres2 = 8;

fac = 0.42; % relative length of G2 compared to G2

L2 = fac * L_global;
Xres2 = ceil(Xres_global * fac); % ceil ensures horizontal cell size for G2 is smaller than for G1


%% First attempt at making a glued grid

G1 = cartGrid([Xres_global, Yres1], [L_global, H1]);

G2 = cartGrid([Xres2, Yres2], [L2, H2]);

% shifting grid G2 upwards (on top of G1)
G2.nodes.coords(:,2) = G2.nodes.coords(:,2) + max(G1.nodes.coords(:,2));

% The following glued grid works, but it is not ideal, due to the poor
% matching of node positions along the glued interface.
G3 = glue2DGrid(G1, G2);


%% Second attempt at making a glued grid

dx_global = L_global / Xres_global; % cell size of G1
L2_tmp = Xres2 * dx_global; % grid size conforming with G1

G2_better = cartGrid([Xres2, Yres2], [L2_tmp, H2]);

% shifting grid G2 downwards
G2_better.nodes.coords(:,2) = G2_better.nodes.coords(:,2) + max(G1.nodes.coords(:,2));

% shifting right-side nodes to comply with the original length requirement
g2_node_xpos = reshape(G2_better.nodes.coords(:,1), Xres2+1, []);
g2_node_xpos(end, :) = L2;
G2_better.nodes.coords(:,1) = g2_node_xpos(:);

% The following glued grid has a perfect match among faces wherever possible
G3_better = glue2DGrid(G1, G2_better);

%% Show the two glued grids next to each other
figure()
subplot(1,2,1); plotGrid(G3); title('Nonconforming faces.');
subplot(1,2,2); plotGrid(G3_better); title('conforming faces.');

%set(gcf, 'position', [118 2409 1113 396]);