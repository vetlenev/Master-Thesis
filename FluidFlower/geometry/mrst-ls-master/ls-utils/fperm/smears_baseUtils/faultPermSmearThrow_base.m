%% Compute permeability of a fault with clay smears
%
% In the configuration used here, FW on the left and HW is on the right.
%
% The smears are consistent with (1) the number of sand and clay layers at
% each side of the fault; (2) the thickness of these layers; (3) fault 
% displacement and thickness; (4) the residual friction angle of the sands
% and clays; (5) the SSFc; and (6) the maximum length of a smear segment in
% each subdomain (to control whether smears are broken or not).
% Based on these inputs (1) to (6), the following are calculated: (I.) the
% probability of encountering a cell with smear in each smear subdomain;
% (II.) the length of each smear segment within each subdomain; and (III.)
% the thickness of each smear (and accordingly, the thickness of each
% subdomain).
%
%__________________________________________________________________________
clc, clear, close all
mrstModule add ls-proj ls-utils coarsegrid upscaling incomp mimetic mpfa ...
           streamlines mrst-gui

% Same result in the stochastic smear placement simulation
% for n=1:100
% rng('default');
% rng(n);
%rng(13)

%% 1. Initialize
% Input independent variables
% Note: fw and hw layers' thickness (T) is the true layer thickness
%       to agree with Egholm et a. (2008). Relation to the apparent thickness
%       is Tap = T/(cosBl*sinB) where Bl is layer dip angle and B is
%       fault dip angle.
% Note: Collapse T and isclay in both fw and hw so that adjacent layers of
%       the same lithology are just one (eg [0 1] and not [0 1 1 1])
f    = struct('T', 2, 'dip', 60);           % fThickness not from grid
fw   = struct('T', [10 30 20 10 20 10], 'isclay', repmat([true true], 1, 3), 'Bl', 0);
hw   = struct('T', [5 20 25 25 10 15], 'isclay', repmat([true true], 1, 3), 'Bl', 0);
sPerms = repmat([20, 20], 1, 1);          % 1x2 or N(input layer) x 2
cPerms = [1e-4 1e-3; 1e-5 1e-4; 1e-2 1e-1; 1e-5 1e-4; 1e-4 1e-3; 1e-3 1e-2; 1e-4 1e-3; 1e-5 1e-4; 1e-2 1e-1; 1e-5 1e-4; 1e-4 1e-3; 1e-3 1e-2];      % "
sand = struct('perm', sPerms*milli*darcy, 'poro', 0.3, 'phi', 30);
clay = struct('perm', cPerms*milli*darcy, 'poro', 0.1, 'phi', [10 5 8 12 3 16 10 5 8 4 3 5], 'SSFc', [4 3 6 2 10 8 4 3 6 5 10 8]);
g    = struct('resL', 1, 'resT', 0.1);          % grid resolution [m]
smear.inL  = {'Lsmear', 20, 'Lsmear', 10, 'Lsmear', 'Lsmear', 'Lsmear', 20, 'Lsmear', 10, 'Lsmear', 'Lsmear'};                        % MUST BE CELL. Smear segment length [m] 
smear.modeledAs = 'objects';                    % algorithm to place smears
smear.tolerance = 0.001;                        % tolerance in P matching  
displayFigures = [0 1 1];                       % [grid, smear location, perm] 

%% 2. Create and populate grid
% Cartesian grid in MRST has indexing that starts at (0, 0) coordinate
% and runs faster on the X axis. Mapping grid and values for properties 
% are stored in structure M.

tic
% 2.1 Compute initial values
[f, sand, clay, fw, hw, smear] = smearOrganizeIn(f, sand, clay, fw, hw, smear);

% 2.2 Create and plot grid
nels = max([round(f.T/g.resT), round(f.D/g.resL)]);
G = computeGeometry(cartGrid([nels, nels], [f.T f.D]));
G.xzFaceDim = [f.T/nels, f.D/nels];

% f1 = figure(1);
% plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1); 
% axis equal; colorbar; xlim([0 f.T]); ylim([0 f.D]); 
% xlabel('x [m]'), ylabel('y [m]')
% title(['Number of cells = ' num2str(G.cells.num)])
% set(f1, 'position', [400, 0, 350, 900]); view([-30, 75]) 

% 2.3 Mapping matrix structure
M = smearMap(G, fw, hw, f, sand, clay, smear);

% 2.4 Smears (final 1s)
%     Continuity along possible diagonals 
[M, smear] = smearContinuity(fw, f, M, sand, clay, smear);

%     Smear distribution
if any(M.Psmear < 1) % place objects (segments) of Lsmear, and match P(smear)
    assert(strcmp(smear.modeledAs, 'objects'));
    M = smearObjects(M, smear, f, G);
elseif isempty(M.Psmear)
    disp('No smear: P(smear) = 0')
else
    disp('Continuous smears: P(smear) = 1')
end

% assign porosity and permeability to each grid cell
[poroG, permG, k] = poroPermToGrid(G, M, f, sand, clay);

%% 3. Compute equivalent/upscaled permeability and porostiy
% Porosity (additive)
% Here, for 1 cell in upscaled grid and regular cells we could just do 
% Poro = (sum(sum(M.vals))*clay.poro +  ...
%        (G.cartDims(1)^2 - sum(sum(M.vals)))*sand.poro) / G.cartDims(1)^2;
p = partitionCartGrid(G.cartDims, [1 1]);
CG = generateCoarseGrid(G, p);
Poro = accumarray(p, poroG)./accumarray(p,1);

% Upscaled 2D Permeability tensor (Perm is not an additive property)
gcellsCheck      = 0;               % check if Perm obtained with grid with AR of only 2 gives same Perm
kbc              = 'sealed';        % 'sealed' or 'open' (periodic)
kmethod          = 'mpfa';          % tpfa, mpfa or mimetic
kacc             = 1;               % use mex acceleration in mpfa
outfl            = 0;               % compare outflux of fine and coarse models
outUps           = [1 10];          % [nx, ny] in coarse grid to compare outflow
plotFigs         = 0;  
[K] = smearPermUpscaled(G, permG, f, plotFigs, gcellsCheck, kbc, kmethod, kacc, outfl, outUps);

% Compute perm in third dimension and save the 3 independent components
smearCells = sum(sum(M.vals));
nct = G.cells.num;
kyy = (smearCells/nct)*clay.perm(2) + ((nct - smearCells)/nct)*sand.perm(2);     
Perm = [K(1), kyy, K(4)];
% if Perm(1)/(milli*darcy) < 9e-2 && Perm(1)/(milli*darcy) > 1e-3
%    disp(num2str(n))
%    break
% end
toc

% Plots
smearFigs(displayFigures, G, M, f, fw, hw, poroG, permG)

%end
%pause
%close all force
%save('resL7.5_10beds_fT2fdip60_phi1030_SSFc6_LsSmear.mat')
