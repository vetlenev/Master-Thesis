%% Black Oil Test: fluid properties
% Based on: singlePhaseAD.m + black oil examples
%
% The purpose of this test model is to illustrate how to add saturation
% regions and capillary pressures either from .DATA input files (tabulated)
% or from a more complex function defined in the program.

clc, clear, close all force
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui linearsolvers ...
               ls-tests


%% Option for capillary pressure
pcToUse = 'function';    % 'input' (the one in .DATA file) or 'function' (customized)

           
%% Set up model geometry
meshpath = 'mrst-2019a/myprojects/2.5Dtestmesh/nodes_coordinates.dat';
vertices = dlmread(meshpath);
x = vertices(:,1); y = vertices(:, 2).*-1;
% delaunay triangulation and extrude grid
%thick = [50 linspace(50,300,15)];
thick = [linspace(200,20,10) 15 linspace(20,200,10)];
G = makeLayeredGrid(triangleGrid([x, y]), thick);
% Permute axes
G.nodes.coords = G.nodes.coords(:,[3 1 2]);
% Compute grid geometry
G = computeGeometry(G);
% Gravity
gravity reset on
% Plot Grid
plotGrid(G, 'FaceColor', 'k', 'FaceColor', [255 255 204]./255, 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); set(gca,'Ydir','reverse')
view([-125, 25]); axis equal tight; 

%% Define rock & fluid model
% Define regions
% ------------------------ IMPORTANT NOTE: -------------------------------
% If you use an indexing vector that just contains the cell
% indices of the fault cells, make sure it is sorted from smaller to larger 
% before using assignFaultPc.
% -------------------------------------------------------------------------
id_seal = G.cells.centroids(:,3) < 300; % i.e. top "seal" layer (= caprock)

% Poro and perm
seal_satRegNum = 2;
rock = makeRock(G, 100*milli*darcy, 0.3); % homogeneous perm
rock.perm(id_seal, 1) = 0.1*milli*darcy;  % if permeability is very low, s = 0 and pc = 0
rock.poro(id_seal) = 0.17;
rock.regions.saturation = ones(G.cells.num,1); 
rock.regions.saturation(id_seal) = seal_satRegNum;

%rock.regions.rocknum = ones(G.cells.num,1); rock.regions.rocknum(id_seal) = 2; % only for keyword ROCK, not required for ROCKTAB
%rock.perm(idc)=1*milli*darcy; 
% note that rock.regions implemented fields are just 'saturation',
% 'pvt', 'imbibition', and 'surfactant'. For rock compressibility regions
% that differ from pvt regions, assignROCK is modified similar to
% assignROCKTAB so that regions are correctly picked when evaluating
% fluid.pvMultR(p) directly.

% For the deck, we create an ECLIPSE deck file which contains PROPS to use. 
% The remaining keywords are just to avoid errors in reading the deck, but 
% are not used here.
fn = 'ls-tests/2.5D_testMeshAndProps/regAndPcForJosimar/ftest_co2brine.DATA';
deck = convertDeckUnits(readEclipseDeck(fn));

% Indicate regions for both fluid and rock
%rockn = ones(G.cells.num,1); rockn(id_seal)=2; % it can be different from satn. TEST
%deck.REGIONS.ROCKNUM = rockn;

% Initialize fluid structure
fluid = initDeckADIFluid(deck);

% Do we need to modify/add capillary pressure?
if strcmp(pcToUse, 'function')
    % This could be any function you want. Here:
    % assignFaultPc requires a reference Pc value for a reference poro and
    % perm, and then the scaling is performed for each cell in the region
    % accordingly. In order to be conservative, you would take the
    % reference Pc for the caprock (i.e. pure clay), for which you can use
    % the values in the .DATA input file (second table under SGOF). 
    refPc.val  = fluid.pcOG{seal_satRegNum};       % you would use the id for fault cells
    refPc.poro = 0.15;  
    refPc.perm = 0.01*milli*darcy; 
    
    fault.poro = rock.poro(id_seal);                % must be passed in ascending order of global fault cell id
    fault.perm = rock.perm(id_seal);                % "
    fault.satRegNum = 2;                       
    
    fluid = assignFaultPcSimple(fluid, fault, refPc);   
end


%% Wells
simTime = 20*year;

numwell = 1;                                  % 4 horizontal injectors
injmass = 10^9*(simTime/year);                % 1 Mt/year of gas
injrate = injmass/(fluid.rhoGS(1)*simTime*numwell);
wellLoc = [sum(thick)/2, 2360, 500]; 
dirw = 'z';
distx = pdist2(G.cells.centroids(:,1), wellLoc(:,1));
disty = pdist2(G.cells.centroids(:,2), wellLoc(:,2));
distz = pdist2(G.cells.centroids(:,3), wellLoc(:,3));
dist = sqrt(distx.^2 + disty.^2 + distz.^2);
[~, wellInx] = min(dist,[],1);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', dirw, ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1]); % order is always 'WOG' ('OG' in GenericBlackOilModel)

% 10y injection = simTime
timesteps = [linspace(0.01, 0.01, 10) linspace(0.1, 0.1, 9) ...
             linspace(1, 1, 9) linspace(2, 2, 5)]*year;
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

%% Initialize state
% Need to initialize 'p', 's', 'rs' and 'rv'
% p: impose vertical equilibrium (fully saturated in oil)
g = norm(gravity);
rho_wr = 1020*kilogram/meter^3;
p_r = g*rho_wr*50;                         % p(z0) below 50 m of water col.
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);            % no dissolved gas at the beginning
rv0 = 0;                                % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);


%% Plot well and initial pressure
clf
show = true(G.cells.num,1);
show(wellInx) = false;
plotCellData(G,p0/barsa, show,'EdgeColor','k', 'FaceAlpha', 0.8);
c = colorbar; %xlim([0 5025]); ylim([0 5000]); zlim([0 1000]);
c.Label.String = 'p [bar]';
plotWell(G,W, 'height',10);
set(gca,'Xdir','reverse')
%set(gca,'Ydir','reverse')
view([55,25]), camproj perspective 
axis equal tight;


%% Correct Pc and Assign model (Two phase live oil - dry gas)
% Pore volume
%fluid = assignPvMult(fluid, deck.PROPS.ROCK(:, 2), rock.regions.rocknum);

% Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);

% Nonlinear solver and threads
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);    %MEX is giving me an error in mexFaceAverageDiagonalJac
model = model.validateModel();
numComponents = 2;
solver  = AMGCL_CPRSolverAD('tolerance', 1e-6, 'block_size', numComponents, ...
                            'useSYMRCMOrdering', true, ...
                            'coarsening', 'aggregation', 'relaxation', 'ilu0');
nls = NonLinearSolver('LinearSolver', solver);

% Modifications to FlowPropertyFunctions, if any, must be done AFTER setting the autodiff
% backend (in setAcceleration_bo) since a copy of the backend will be made
% when setting the flow property functions.
%model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));

% This gives an error, so we get the outputs at the end (before figures)
%model.OutputStateFunctions = {'ComponentTotalMass', 'RelativePermeability', 'CapillaryPressure'};


%% Boundary conditions
cells = 1:G.cells.num;
cells = cells(~id_seal);
idx = [G.cells.centroids(cells, 1) < 110 ...
       G.cells.centroids(cells, 1) > 2110];
idy = [G.cells.centroids(cells, 2) < 100 ...
       G.cells.centroids(cells, 2) > 4930];
ids = any([any(idx, 2) any(idy, 2)], 2);
cellsext = cells(ids);
model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e3;


% Plot BC cells
%clf, plotGrid(G,'FaceColor', 'none'); set(gca,'Xdir','reverse'); 
%view([55,25]), camproj perspective; axis equal tight;
%plotCellData(G, p0, cellsext); colorbar; xlabel('x [m]')
%pause

%% Simulate schedule
% Add schdule
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', []);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % field 2 empty well (after injection)
schedule.control(2).W.val = 0;
schedule.step.control(cumsum(schedule.step.val) > 5*year) = 2;              % inject 5 years only

% Simulate schedule
N = 4;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N; 
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

for n=1:numel(states)
    states{n}.dp = states{n}.pressure - p0;  
    states{n}.FlowProps.CapillaryPressure = model.getProp(states{n}, 'CapillaryPressure');
    states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
end
states{1}.reg = rock.regions.saturation;

cnum = 1:G.cells.num;
figure(10)
plotToolbar(G, states, cnum(G.cells.centroids(:,1)>wellLoc(1))); set(gca,'Xdir','reverse')
view(55, 25); camproj perspective; axis equal tight 
plotWell(G,W); c = colorbar; caxis([0 1])

% figure(11)
% plotCellData(G,states{end}.rs, 'faceAlpha', 0.7);
% set(gca,'Xdir','reverse')
% view(55, 25); camproj perspective; axis equal tight 
% c = colorbar; %C = flipud(load('colormap.txt')); colormap(C/255);
%caxis([0 90]); cmocean('thermal');