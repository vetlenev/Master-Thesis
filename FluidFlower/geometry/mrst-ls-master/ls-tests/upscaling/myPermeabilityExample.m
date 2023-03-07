%% Upscale a problem with a diagonal trend
% Based on permeabilityExample1.m, permeabilityExample2.m,
% upscalePermeabilityPeriodic.m, 
% In this example we will upscale a problem with diagonal trend and compare
% the results obtained with a flow-based method with pressure-drop and
% periodic boundary conditions (incompTPFA and Mimetic solvers).

clear, clc, close all

mrstModule add incomp upscaling mimetic ls-utils;

permType = 'isotropic';

%% Set up model
[Lx,Ly] = deal(10);
[nx,ny] = deal(20);
G = computeGeometry(cartGrid([nx ny], [Lx Ly]));
[i, j] = gridLogicalIndices(G);
%ids = sin(2*sum(G.cells.centroids(:,:),2)*pi/Lx)>0;
%ids=false(ny);
%ids = full(spdiags(true(ny, 13), -3:3, ids));
%ids(6:20, :) = true; %ids(:, 8:13) = true;
ids = rand(ny)>0.75;
ids = reshape(transpose(ids), G.cells.num, 1);
switch permType
    case 'isotropic'
        rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]); 
        rock.perm(ids) = 0.01*milli*darcy;
        ikp = 1;
    case 'anisotropic'
        rock.perm = repmat([30 -2.5 100]*milli*darcy, G.cells.num, 1);
        rock.perm(ids, :) = repmat([1.1e-4 -3.5e-4 1e-2]*milli*darcy, sum(ids), 1);
        ikp = [1 2];
end
rock.poro = repmat(.3, [G.cells.num, 1]);
hT = computeTrans(G, rock);

% Coarse grid
upscaled = [1 1];
G_ups = computeGeometry(cartGrid(upscaled, [Lx Ly]));
p = partitionCartGrid(G.cartDims, upscaled);
CG = generateCoarseGrid(G, p);

% Fluid
fluid = initSingleFluid('mu' , 1, 'rho', 1);

% Plot
clf, plotCellData(G,rock.perm(:,1)/(milli*darcy),'EdgeColor','k'); colorbar

% Structures with boundary conditions
d = G.griddim;
[bcl,bcr, Dp]=deal(cell(d,1));
bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
for j = 1:d;
   bcl{j} = pside([], G, bcsides{j, 1}, 0);
   bcr{j} = pside([], G, bcsides{j, 2}, 0);
   Dp{j}  = 0;
end
Dp{1} = 4*barsa;
L  = max(G.faces.centroids)-min(G.faces.centroids);


%% Pressure-drop boundary conditions
[v,dp] = deal(zeros(d, 1));
for i=1:d
   bc = addBC([], bcl{i}.face, 'pressure', Dp{1});
   bc = addBC(bc, bcr{i}.face, 'pressure', Dp{2});

   xr = initResSol(G, 1*barsa, 1);
   xr = incompTPFA(xr, G, hT, fluid, 'bc', bc);

   v(i)  = sum(xr.flux(bcr{i}.face)) / ...
           sum(G.faces.areas(bcr{i}.face));
   dp(i) = Dp{1}/L(i);
end
Kpd = convertTo(v./dp, milli*darcy);      


% Mimetic
Kpdm = convertTo(upscalePerm(G, CG, rock, 'method', 'mimetic'), milli*darcy);                                 


%% Linear boundary conditions
bfac = boundaryFaces(G);
fSide = bfac(G.faces.centroids(bfac, 1) < .1);
fSideCentr = G.faces.centroids(fSide, 2);
pSideVals = Dp{1} - fSideCentr.*(Dp{1}-Dp{2})./Ly;

[v,dp] = deal(zeros(d, 1));
idxs = 1:d;
for i=idxs
   bc = addBC([], bcl{i}.face, 'pressure', Dp{1});
   bc = addBC(bc, bcr{i}.face, 'pressure', Dp{2});
   ids = idxs(idxs~=i);
   bc = addBC(bc, bcl{ids}.face, 'pressure', pSideVals);
   bc = addBC(bc, bcr{ids}.face, 'pressure', pSideVals);
   
   xr = initResSol(G, 1*barsa, 1);
   xr = incompTPFA(xr, G, hT, fluid, 'bc', bc);
   
   v(i,i) = sum(xr.flux(bcr{i}.face)) / sum(G.faces.areas(bcr{i}.face));
   v(i,ids) =  sum(xr.flux([bcr{ids}.face; bcl{ids}.face])) ./ ...
               sum(G.faces.areas([bcr{ids}.face; bcl{ids}.face]));
   %v(i)  = sum(xr.flux(bcr{i}.face)) / ...
   %        sum(G.faces.areas(bcr{i}.face));
   dp(i) = Dp{1}/L(i);
end
Kl = convertTo(v./dp, milli*darcy);
assert(all(eig(Kl)>0));


%% Periodic boundary conditions
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, Dp);
ofaces = cell(d,1);
for j=1:d, ofaces{j} = bcp.face(bcp.tags==j); end
v  = nan(d);
dp = Dp{1}*eye(d);
nbcp = bcp;

% TPFA
for i=1:d
   for j=1:d, nbcp.value(bcp.tags==j) = dp(j,i); end

   xr = initResSol(Gp, 100*barsa, 1);
   xr = incompTPFA(xr, Gp, hT, fluid, 'bcp', nbcp);

   for j=1:d
      v(j,i) = sum(xr.flux(ofaces{j})) / ... *bcp.sign(bcp.tags==j)) / ...
               sum(Gp.faces.areas(ofaces{j}));
   end
end
dpf = bsxfun(@rdivide, dp, L);
K = convertTo(v/dpf, milli*darcy);                                         
[V,D] = eig(K);                                                            

% Mimetic solver
Sp = computeMimeticIPGp(G, Gp, rock, 'InnerProduct', 'ip_tpf'); % to match TPFA
vm = nan(d);
bcp_new = bcp;
for i=1:d
    for j=1:d
        bcp_new.value(bcp.tags==j)=dp(j,i);%/L(j);
    end
    xrm = initResSol(Gp, 100*barsa, 1);
    xrm = incompMimetic(xrm, Gp, Sp, fluid, 'bcp', bcp_new);
    for j=1:d
        flux=sum(xrm.flux(ofaces{j}).*bcp.sign(bcp.tags==j));
        vm(j,i)=flux/sum(Gp.faces.areas(ofaces{j}));
    end
end
Km = convertTo(-vm/dpf, milli*darcy);                                         
[Vm,Dm] = eig(Km);

% Mimetic solver with function upscalePermeability.m
Sp2 = computeMimeticIPGp(G, Gp, rock); 
psolver = @(state0, Gp, fluid, bcp) incompMimetic(state0, Gp, Sp2, fluid, 'bcp', bcp);
Kmm = convertTo(myupscalePermeabilityPeriodic(Gp, bcp, Dp{1}, psolver, fluid, L), milli*darcy);
[Vmm,Dmm] = eig(Kmm);


%% Compare outflow
% Define coarse grid permeability
% PERMEABILITIES ALWAYS IN m^2 to simulate!
crock{1}.perm = repmat(Kpd(ikp)'*milli*darcy, CG.cells.num, 1);               % TPFA p drop
crock{2}.perm = repmat(Kpdm(ikp)*milli*darcy, CG.cells.num, 1);               % Mimetic p drop
crock{3}.perm = repmat([K(1) K(2) K(4)]*milli*darcy, CG.cells.num, 1);        % TPFA periodic
crock{4}.perm = repmat([Kmm(1) Kmm(2) Kmm(4)]*milli*darcy, CG.cells.num, 1);  % Mimetic periodic

% Define fluid
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3); % water

% Define initial conditions
state0 = initResSol(G, 0.5*barsa, 1);            % note random p, here just x,y coordinates
state0c = initResSol(G_ups, 0.5*barsa, 1);

% 1. Pressure drop with no-flow in other boundaries
% Fine-scale problem
bc        = pside([], G, 'North', 0);
faces     = bc.face;
bc        = pside(bc, G, 'South',  1*barsa());
x         = incompTPFA(state0, G, hT, fluid, 'bc', bc);

S         = computeMimeticIP(G, rock);
xm        = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale problem
bc_ups    = pside([], G_ups, 'North', 0);
faces_ups = bc_ups.face;
bc_ups    = pside(bc_ups, G_ups, 'South', 1*barsa());
T_ups     = computeTrans(G_ups, crock{1});
x_ups     = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);

Sc        = computeMimeticIP(G_ups, crock{2});
x_upsm    = incompMimetic(state0c, G_ups, Sc, fluid, 'bc', bc_ups);

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

% 2. P drop, open boundaries at constant p

% Fine-scale
pTop = 0.1*barsa();
pBot = 2*barsa();
pSideAvg = pTop + (pBot-pTop)/2;
bfac = boundaryFaces(G);
fSide = bfac(G.faces.centroids(bfac, 1) < .1);
fSideCentr = G.faces.centroids(fSide, 2);
pSideVals = flipud(pTop  + fSideCentr.*(pBot-pTop)./Ly);

bc = pside([], G, 'North', pTop);
bc = pside(bc, G, 'East', pSideVals);
bc = pside(bc, G, 'West', pSideVals);
bc = pside(bc, G, 'South', pBot);

x  = incompTPFA(state0, G, hT, fluid, 'bc', bc);
xm = incompMimetic(state0, G, S, fluid, 'bc', bc);

% Coarse-scale
bc_ups = pside([], G_ups, 'North', pTop);
bc_ups = pside(bc_ups, G_ups, 'East', pSideAvg);
bc_ups = pside(bc_ups, G_ups, 'West', pSideAvg);
bc_ups = pside(bc_ups, G_ups, 'South', pBot);

x_ups{1} = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
T_ups2   = computeTrans(G_ups, crock{3});
x_ups{2} = incompTPFA(state0c, G_ups, T_ups2, fluid, 'bc', bc_ups);
S_ups    = computeMimeticIP(G_ups, crock{4});
x_upsm   = incompMimetic(state0c, G_ups, S_ups, fluid, 'bc', bc_ups);

flux1 = sum(x.flux(faces));
flux2 = sum(xm.flux(faces));
flux3 = sum(x_ups{1}.flux(faces_ups));
flux4 = sum(x_ups{2}.flux(faces_ups));
flux5 = sum(x_upsm.flux(faces_ups));
disp('Constant flux with open boundaries: ')
disp(['Sum outflux on fine scale (TPFA): ', num2str(flux1)]);
disp(['Sum outflux on fine scale (Mimetic): ', num2str(flux2)]);
disp(['Sum outflux on coarse scale (TPFA with pdrop (diagonal) k): ', num2str(flux3)]);
disp(['Sum outflux on coarse scale (TPFA with periodic TPFA k): ', num2str(flux4)]);
disp(['Sum outflux on coarse scale (Mimetic with periodic mimetic k): ', num2str(flux5)]);