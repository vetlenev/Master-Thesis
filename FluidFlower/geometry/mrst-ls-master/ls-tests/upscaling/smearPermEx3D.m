%% Upscale a problem with a diagonal trend
% Based on permeabilityExample1.m, permeabilityExample2.m,
% upscalePermeabilityPeriodic.m, 
% In this example we will upscale a problem with diagonal trend and compare
% the results obtained with a flow-based method with pressure-drop and
% periodic boundary conditions (incompTPFA and Mimetic solvers).

clear, clc, close all

mrstModule add incomp upscaling mimetic ls-utils;

permType = 'anisotropic';
topDepth = 1000*meter;

%% Set up model
[Lx,Ly,Lz] = deal(5,5,50);
[nx,ny,nz] = deal(50,10,50);
G = computeGeometry(cartGrid([nx ny nz], [Lx Ly Lz]));
G.cells.centroids(:,3) = G.cells.centroids(:,3)+topDepth;
G.faces.centroids(:,3) = G.faces.centroids(:,3)+topDepth;
G.nodes.coords(:,3) = G.nodes.coords(:,3)+topDepth;
L = max(G.faces.centroids)-min(G.faces.centroids);

% cells with low and high perm
ids = false(nx);
ids = full(spdiags(true(nx, 20), 30:49, ids));
ids = full(spdiags(true(nx, 25),  -25:-1, ids));
ids = reshape(transpose(ids), G.cartDims(1)*G.cartDims(3), 1);
xyNum = G.cartDims(1)*G.cartDims(2);
id3 = zeros(G.cells.num,1);
for n=1:G.cartDims(3)
    ii = (n-1)*xyNum+1:n*xyNum;
    idsXZ = (n-1)*G.cartDims(1)+1:n*G.cartDims(1);
    id3(ii) = repmat(ids(idsXZ),G.cartDims(2),1);
end
id3 = logical(id3);

switch permType
    case 'isotropic'
        rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]); 
        rock.perm(id3) = 0.1*milli*darcy;
    case 'anisotropic'
        rock.perm = repmat([30 -2.5 100]*milli*darcy, G.cells.num, 1);
        %rock.perm = repmat([1.1e-4 -3.5e-4 1e-2]*milli*darcy, G.cells.num, 1);
        rock.perm(id3, :) = repmat([1.1e-4 -3.5e-4 1e-2]*milli*darcy, sum(id3), 1);
end
rock.poro = repmat(.3, [G.cells.num, 1]);
hT = computeTrans(G, rock);

% Coarse grid
upscaled = [1 1 5];
G_ups = computeGeometry(cartGrid(upscaled, [Lx Ly Lz]));
G_ups.cells.centroids(:,3) = G_ups.cells.centroids(:,3)+topDepth;
G_ups.faces.centroids(:,3) = G_ups.faces.centroids(:,3)+topDepth;
G_ups.nodes.coords(:,3) = G_ups.nodes.coords(:,3)+topDepth;
p = partitionCartGrid(G.cartDims, upscaled);
CG = generateCoarseGrid(G, p);

% Fluid
fluid = initSingleFluid('mu' , 1, 'rho', 1);

% Plot
figure(1), subplot(1,2,1)
plotCellData(G,log10(rock.perm(:,1)/(milli*darcy)),'EdgeColor','none'); colorbar
xlim([0 Lx]), ylim([0 Ly]), zlim([topDepth topDepth+Lz])
axis equal, view([-30, 15]), xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')


%% Pressure-drop boundary conditions
% TPFA
Kpd = convertTo(upscalePerm(G, CG, rock, 'method', 'tpfa'), milli*darcy);          

% Mimetic
Kpdm = convertTo(upscalePerm(G, CG, rock, 'method', 'mimetic'), milli*darcy);                                 


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

% % TPFA
psolver = @(state0, Gp, fluid, bcp) incompTPFA(state0, Gp, hT, fluid, 'bcp', bcp);
K = convertTo(myupscalePermeabilityPeriodic(Gp, bcp, Dp{1}, psolver, fluid, L), milli*darcy);
[V,D] = eig(K);
assert(all(eig(K)>0))

% Mimetic solver with function upscalePermeability.m
Sp2 = computeMimeticIPGp(G, Gp, rock); 
psolver = @(state0, Gp, fluid, bcp) incompMimetic(state0, Gp, Sp2, fluid, 'bcp', bcp);
Kmm = convertTo(myupscalePermeabilityPeriodic(Gp, bcp, Dp{1}, psolver, fluid, L), milli*darcy);
[Vmm,Dmm] = eig(Kmm);
assert(all(eig(Kmm)>0))               % perm must be positive definitexlim([0 Lx])

%% Compare outflow
% Define coarse grid permeability
% PERMEABILITIES ALWAYS IN m^2 to simulate!
crock{1}.perm = Kpd*milli*darcy;                                                         % TPFA p drop
crock{2}.perm = Kpdm*milli*darcy;                                                        % Mimetic p drop
%crock{3}.perm = repmat([K(1) K(2) K(3) K(5) K(6) K(9)]*milli*darcy, CG.cells.num, 1);    % TPFA periodic
crock{4}.perm = repmat([Kmm(1) Kmm(2) Kmm(3) Kmm(5) Kmm(6) Kmm(9)]*milli*darcy, CG.cells.num, 1);            % Mimetic periodic

% Define fluid
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);   % water

% Define initial conditions (gravity ON downwards from z!)
gravity reset on
watCol = 50*meter;
g = norm(gravity);
p_r    = (watCol + min(G.cells.centroids(:,3)))*fluid.rhoWS*g;
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
state0 = initResSol(G, p0, 1);   

figure(1); subplot(1,2,2)
plotCellData(G, p0); colorbar
axis equal, view([-45, 15]), xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')


p_s    = (watCol + min(G_ups.cells.centroids(:,3)))*fluid.rhoWS*g;
[z_0u, z_maxu] = deal(min(G_ups.cells.centroids(:,3)), max(G_ups.cells.centroids(:,3)));
if z_0u ~= z_maxu
    equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_s);
    p0u = reshape(deval(equil, G_ups.cells.centroids(:,3)), [], 1);  clear equil
elseif G_ups.cells.num == 1
    p0u = (watCol + G_ups.cells.centroids(:,3))*fluid.rhoWS*g;
end
state0c = initResSol(G_ups, p0u, 1);
gravity reset

% 1. Pressure drop with no-flow in other boundaries
% Fine-scale problem

bc        = pside([], G, 'Top', 0.1*min(p0));
faces     = bc.face;
bc        = pside(bc, G, 'Bottom',  10*max(p0));
x         = incompTPFA(state0, G, hT, fluid, 'bc', bc);

S         = computeMimeticIP(G, rock);
xm        = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale problem
bc_ups    = pside([], G_ups, 'Top', 0.1*min(p0));
faces_ups = bc_ups.face;
bc_ups    = pside(bc_ups, G_ups, 'Bottom', 10*max(p0));
T_ups     = computeTrans(G_ups, crock{1});
x_ups     = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);

Sc        = computeMimeticIP(G_ups, crock{2});
x_upsm    = incompMimetic(state0c, G_ups, Sc, fluid, 'bc', bc_ups);


% Note that for this case, since there are very low permeability layers in
% diagonal, a p drop in the vertical direction will generate most of the
% outflow through the lateral boundaries. 
% Negative flux is outflow, i.e. out of the domain.
flux1 = sum(x.flux(faces));
flux2 = sum(xm.flux(faces));
flux3 = sum(x_ups.flux(faces_ups));
flux4 = sum(x_upsm.flux(faces_ups));
disp(' Pressure drop with sealed boundaries: ')
disp(['Sum outflux on fine scale (TPFA)  : ', num2str(flux1)]);
disp(['Sum outflux on fine scale (Mimetic)  : ', num2str(flux2)]);
disp(['Sum outflux on coarse scale (TPFA) : ', num2str(flux3)]);
disp(['Sum outflux on coarse scale (Mimetic): ', num2str(flux4)]);
disp('___________________________________________________________________')
clear x xm x_ups x_upsm bc bc_ups

% 2. P bottom boundary, open lateral boundary
pTop = 0.1*min(p0);
pBot = 10*max(p0);
bfac = boundaryFaces(G);
fSide = bfac(G.faces.centroids(bfac, 1) < .001);
fSideCentr = G.faces.centroids(fSide, :);

p_s    = (watCol + min(fSideCentr(:,3)))*fluid.rhoWS*g;
[z_0, z_max] = deal(min(fSideCentr(:,3)), max(fSideCentr(:,3)));
equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_s);
pSideVals = reshape(deval(equil, fSideCentr(:,3)), [], 1);  clear equil

bc = pside([], G, 'Top', pTop);
bc = pside(bc, G, 'Bottom', pBot);
bc = pside(bc, G, 'East', pSideVals);
bc = pside(bc, G, 'West', pSideVals);
%bc = pside(bc, G, 'North', pSideVals);
%bc = pside(bc, G, 'South', pSideVals);

x  = incompTPFA(state0, G, hT, fluid, 'bc', bc);
xm = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale
bfacu = boundaryFaces(G_ups);
fSideu = bfacu(G_ups.faces.centroids(bfacu, 1) < .001);
fSideCentru = G_ups.faces.centroids(fSideu, :);

p_su    = (watCol + min(fSideCentru(:,3)))*fluid.rhoWS*g;
[z_0u, z_maxu] = deal(min(fSideCentru(:,3)), max(fSideCentru(:,3)));
if z_0u ~= z_maxu
    equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_su);
    pSideValsu = reshape(deval(equil, fSideCentru(:,3)), [], 1);  clear equil
elseif G_ups.cells.num == 1
    pSideValsu = (watCol + fSideCentru(:,3))*fluid.rhoWS*g;
end

bc_ups = pside([], G_ups, 'Top', pTop);
bc_ups = pside(bc_ups, G_ups, 'Bottom', pBot);
bc_ups = pside(bc_ups, G_ups, 'East', pSideValsu);
bc_ups = pside(bc_ups, G_ups, 'West', pSideValsu);

x_ups{1} = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
%T_ups2   = computeTrans(G_ups, crock{3});
%x_ups{2} = incompTPFA(state0c, G_ups, T_ups2, fluid, 'bc', bc_ups);
T_ups3   = computeTrans(G_ups, crock{4});
x_ups{3} = incompTPFA(state0c, G_ups, T_ups3, fluid, 'bc', bc_ups);
S_ups    = computeMimeticIP(G_ups, crock{4});
x_upsm   = incompMimetic(state0c, G_ups, S_ups, fluid, 'bc', bc_ups);

flux1 = sum(x.flux(faces));
flux2 = sum(xm.flux(faces));
flux3 = sum(x_ups{1}.flux(faces_ups));
%flux4 = sum(x_ups{2}.flux(faces_ups));
flux5 = sum(x_ups{3}.flux(faces_ups));
flux6 = sum(x_upsm.flux(faces_ups));
disp('Pressure drop with open boundaries: ')
disp(['Sum outflux on fine scale (TPFA): ', num2str(flux1)]);
disp(['Sum outflux on fine scale (Mimetic): ', num2str(flux2)]);
disp(['Sum outflux on coarse scale (TPFA with pdrop (diagonal) k): ', num2str(flux3)]);
%disp(['Sum outflux on coarse scale (TPFA with periodic TPFA k): ', num2str(flux4)]);
disp(['Sum outflux on coarse scale (TPFA with periodic mimetic k): ', num2str(flux5)]);
disp(['Sum outflux on coarse scale (Mimetic with periodic mimetic k): ', num2str(flux6)]);
