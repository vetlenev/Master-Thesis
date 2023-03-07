%% Upscale a problem with a diagonal trend
% Based on permeabilityExample1.m, permeabilityExample2.m,
% upscalePermeabilityPeriodic.m, 
% In this example we will upscale a problem with diagonal trend and compare
% the results obtained with a flow-based method with pressure-drop and
% periodic boundary conditions (incompTPFA and Mimetic solvers).

clear, clc, close all

mrstModule add incomp upscaling mimetic ls-utils;

permType = 'anisotropic';

%% Set up model
[Lx,Ly] = deal(5,100);
[nx,ny] = deal(50,50);
G = computeGeometry(cartGrid([nx ny], [Lx Ly]));
L = max(G.faces.centroids)-min(G.faces.centroids);
ids = false(ny);
ids = full(spdiags(true(ny, 20), 10:29, ids));
ids = full(spdiags(true(ny, 25),  -25:-1, ids));
ids = reshape(transpose(flipud(ids)), G.cells.num, 1);
switch permType
    case 'isotropic'
        rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]); 
        rock.perm(ids) = 0.1*milli*darcy;
    case 'anisotropic'
        rock.perm = repmat([30 -2.5 100]*milli*darcy, G.cells.num, 1);
        rock.perm(ids, :) = repmat([1.1e-4 -3.5e-4 1e-2]*milli*darcy, sum(ids), 1);
end
rock.poro = repmat(.3, [G.cells.num, 1]);
hT = computeTrans(G, rock);

% Coarse grid
upscaled = [1 5];
G_ups = computeGeometry(cartGrid(upscaled, [Lx Ly]));
p = partitionCartGrid(G.cartDims, upscaled);
CG = generateCoarseGrid(G, p);

% Fluid
fluid = initSingleFluid('mu' , 1, 'rho', 1);

% Plot
figure(1); subplot(1,2,1)
plotCellData(G,log10(rock.perm(:,1)/(milli*darcy)),'EdgeColor','none'); colorbar
axis equal; xlim([0 Lx]); ylim([0 Ly]);


%% Pressure-drop boundary conditions
% TPFA
Kpd = convertTo(upscalePerm(G, CG, rock, 'method', 'tpfa'), milli*darcy);       

p2 = partitionCartGrid(G.cartDims, [1 1]);
CG2 = generateCoarseGrid(G, p2); 
KpdSingle = convertTo(upscalePerm(G, CG2, rock, 'method', 'tpfa'), milli*darcy);

%% Periodic boundary conditions
% Structures with boundary conditions
d = G.griddim;
[bcl,bcr, Dp]=deal(cell(d,1));
bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
for j = 1:d
   bcl{j} = pside([], G, bcsides{j, 1}, 0);
   bcr{j} = pside([], G, bcsides{j, 2}, 0);
   Dp{j}  = 0;
end
Dp{1} = 4*barsa;

% Periodic grid
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, Dp);

% Mimetic solver with function upscalePermeability.m
Sp2 = computeMimeticIPGp(G, Gp, rock); 
psolver = @(state0, Gp, fluid, bcp) incompMimetic(state0, Gp, Sp2, fluid, 'bcp', bcp);
Kmm = convertTo(myupscalePermeabilityPeriodic(Gp, bcp, Dp{1}, psolver, fluid, L), milli*darcy);
[Vmm,Dmm] = eig(Kmm);
assert(all(eig(Kmm)>0))

%% Compare outflow
% Define coarse grid permeability
% PERMEABILITIES ALWAYS IN m^2 to simulate!
crock{1}.perm = repmat(KpdSingle*milli*darcy, CG.cells.num, 1);               % TPFA p drop
crock{2}.perm = repmat([Kmm(1) Kmm(2) Kmm(4)]*milli*darcy, CG.cells.num, 1);  % Mimetic periodic

% Define fluid
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3); % water

% Define initial conditions
watCol = 50*meter;
topy = 1000*meter;
g = 9.861;
p_r    = (watCol + topy + min(G.cells.centroids(:,2)))*fluid.rhoWS*g;
[z_0, z_max] = deal(min(G.cells.centroids(:,2)), max(G.cells.centroids(:,2)));
equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_r);
p0 = flipud(reshape(deval(equil, G.cells.centroids(:,2)), [], 1));  clear equil
state0 = initResSol(G, p0, 1);   

figure(1); subplot(1,2,2)
plotCellData(G, p0,'EdgeColor','none'); c = colorbar;
axis equal, xlabel('x [m]'), ylabel('y [m]');
xlim([0 Lx]); ylim([0 Ly]);

p_s    = (watCol + topy + min(G_ups.cells.centroids(:,2)))*fluid.rhoWS*g;
[z_0u, z_maxu] = deal(min(G_ups.cells.centroids(:,2)), max(G_ups.cells.centroids(:,2)));
if z_0u ~= z_maxu
    equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_s);
    p0u = flipud(reshape(deval(equil, G_ups.cells.centroids(:,2)), [], 1));  clear equil
elseif G_ups.cells.num == 1
    p0u = (watCol + topy + G_ups.cells.centroids(:,2))*fluid.rhoWS*g;
end
state0c = initResSol(G_ups, p0u, 1);

% 1. Pressure drop with no-flow in other boundaries
% Fine-scale problem
bc        = pside([], G, 'North', 0.1*min(p0));
faces     = bc.face;
bc        = pside(bc, G, 'South',  10*max(p0));
x         = incompTPFA(state0, G, hT, fluid, 'bc', bc);

S         = computeMimeticIP(G, rock);
xm        = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale problem
bc_ups    = pside([], G_ups, 'North', 0.1*min(p0));
faces_ups = bc_ups.face;
bc_ups    = pside(bc_ups, G_ups, 'South', 10*max(p0));

T_ups     = computeTrans(G_ups, crock{1});
x_ups{1}  = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
T_ups2    = computeTrans(G_ups, crock{2});
x_ups{2}  = incompTPFA(state0c, G_ups, T_ups2, fluid, 'bc', bc_ups);
Sc        = computeMimeticIP(G_ups, crock{2});
x_upsm    = incompMimetic(state0c, G_ups, Sc, fluid, 'bc', bc_ups);

flux1 = sum(x.flux(faces));
flux2 = sum(xm.flux(faces));
flux3 = sum(x_ups{1}.flux(faces_ups));
flux4 = sum(x_ups{2}.flux(faces_ups));
flux5 = sum(x_upsm.flux(faces_ups));
disp(' Pressure drop with sealed boundaries: ')
disp(['Sum outflux on fine scale (TPFA)  : ', num2str(flux1)]);
disp(['Sum outflux on fine scale (Mimetic)  : ', num2str(flux2)]);
disp(['Sum outflux on coarse scale (TPFA with p drop, TPFA k) : ', num2str(flux3)]);
disp(['Sum outflux on coarse scale (TPFA with periodic, mimetic k) : ', num2str(flux4)]);
disp(['Sum outflux on coarse scale (Mimetic with periodic, mimetic k): ', num2str(flux5)]);
disp('___________________________________________________________________')
clear x xm x_ups x_upsm bc bc_ups T_ups T_ups2 Sc


% 2. P drop, open boundaries at constant p
% Fine scale
pTop = 0.1*min(p0);
pBot = 10*max(p0);
bfac = boundaryFaces(G);
fSide = bfac(G.faces.centroids(bfac, 1) < .001);
fSideCentr = G.faces.centroids(fSide, :);

p_s    = (watCol + topy + min(fSideCentr(:,2)))*fluid.rhoWS*g;
[z_0, z_max] = deal(min(fSideCentr(:,2)), max(fSideCentr(:,2)));
equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_s);
pSideVals = flipud(reshape(deval(equil, fSideCentr(:,2)), [], 1));  clear equil

bc = pside([], G, 'North', pTop);
bc = pside(bc, G, 'South', pBot);
bc = pside(bc, G, 'East', pSideVals);
bc = pside(bc, G, 'West', pSideVals);

x  = incompTPFA(state0, G, hT, fluid, 'bc', bc);
xm = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale
bfacu = boundaryFaces(G_ups);
fSideu = bfacu(G_ups.faces.centroids(bfacu, 1) < .001);
fSideCentru = G_ups.faces.centroids(fSideu, :);

p_su    = (watCol + topy + min(fSideCentru(:,2)))*fluid.rhoWS*g;
[z_0u, z_maxu] = deal(min(fSideCentru(:,2)), max(fSideCentru(:,2)));
if z_0u ~= z_maxu
    equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_su);
    pSideValsu = flipud(reshape(deval(equil, fSideCentru(:,2)), [], 1));  clear equil
elseif G_ups.cells.num == 1
    pSideValsu = (watCol + topy + fSideCentru(:,2))*fluid.rhoWS*g;
end

bc_ups = pside([], G_ups, 'North', pTop);
bc_ups = pside(bc_ups, G_ups, 'South', pBot);
bc_ups = pside(bc_ups, G_ups, 'East', pSideValsu);
bc_ups = pside(bc_ups, G_ups, 'West', pSideValsu);

T_ups    = computeTrans(G_ups, crock{1});
x_ups{1} = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
T_ups2   = computeTrans(G_ups, crock{2});
x_ups{2} = incompTPFA(state0c, G_ups, T_ups2, fluid, 'bc', bc_ups);
S_ups    = computeMimeticIP(G_ups, crock{2});
x_upsm   = incompMimetic(state0c, G_ups, S_ups, fluid, 'bc', bc_ups);

flux1 = sum(x.flux(faces));
flux2 = sum(xm.flux(faces));
flux3 = sum(x_ups{1}.flux(faces_ups));
flux4 = sum(x_ups{2}.flux(faces_ups));
flux5 = sum(x_upsm.flux(faces_ups));
disp('Pressure drop with open boundaries: ')
disp(['Sum outflux on fine scale (TPFA): ', num2str(flux1)]);
disp(['Sum outflux on fine scale (Mimetic): ', num2str(flux2)]);
disp(['Sum outflux on coarse scale (TPFA with pdrop (diagonal) k): ', num2str(flux3)]);
disp(['Sum outflux on coarse scale (TPFA with periodic mimetic k): ', num2str(flux4)]);
disp(['Sum outflux on coarse scale (Mimetic with periodic mimetic k): ', num2str(flux5)]);