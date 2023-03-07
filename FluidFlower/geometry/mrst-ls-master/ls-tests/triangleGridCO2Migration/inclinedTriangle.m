%% Migration of CO2 at top of formation
% Free-phase CO2 in most "triangles" at the top of the reservoir, even 
% after 200y (i.e. the CO2 does not fully migrate updip and concentrate at
% the highest point of the reservoir (as expected). Rather, it remains in
% all top cells where it arrived at some point during the simulation).

clear, close all

% Modules
mrstModule add ls-tests ls-utils ad-props deckformat mrst-gui ad-core ...
           ad-blackoil linearsolvers
      
% Options
gridType = 'tri';
krHyst   = true;
cellSize = 10;

% Location
folderName = [gridType '_krHyst_' num2str(krHyst) '_cellSize_' num2str(cellSize)];
topDir     = 'C:\Users\lsalo\matlab\sim_data\mrst\triangleGridCO2Migration\';
       
%% DATA
% Gravity (simulate 15 degrees of inclination)
theta = 5;
R = makehgtform('zrotate', -pi*theta/180);
gravity reset on;
gravity y; gravity(gravity()*-1);
gravity(R(1:3,1:3)*gravity().' );

% Mesh
if strcmp(gridType, 'tri') || strcmp(gridType, 'vor')
    meshpath = 'mrst-dev/mrst-ls/ls-tests/triangleGridCO2Migration/nodeCoord.dat';
    triSize = 10;
    vertices = dlmread(meshpath);
    x = vertices(:,1) + abs(min(vertices(:, 1))); 
    y = vertices(:, 2)*-1;
    if strcmp(gridType, 'tri')
        G = computeGeometry(triangleGrid([x, y]));
        valminBC = 100;
        zlim = -2500;
    else
        G = computeGeometry(pebi(triangleGrid([x, y])));
        valminBC = 50;
        zlim = -2450;
    end
elseif strcmp(gridType, 'quad')
    Lx = 5000;
    Ly = 500;
    quadSize = cellSize;
    nx = Lx/quadSize; ny = Ly/quadSize;
    G = cartGrid([nx, ny], [Lx, Ly]);
    G.nodes.coords(:, 2) = (G.nodes.coords(:, 2) + 2400)*-1;
    G = computeGeometry(G);
    valminBC = 100;
    zlim = -2500;
else
    error('Unknown grid type')
end
%plotGrid(G); axis equal;

if strcmp(gridType, 'tri')
    % Modify centroid locations of lower tri in top tri layer, so that their
    % depth is midway between upper (inverted) tri at left and right

    % Find upper (inverted) tri (top boundary)
    ztop = max(G.faces.centroids(:, 2));
    externalFaces = any([G.faces.neighbors(:, 1) == 0, ...
                         G.faces.neighbors(:, 2) == 0], 2);
    externalCells = G.faces.neighbors(externalFaces, :);
    externalCells = externalCells(externalCells > 0);
    topCells = externalCells(G.cells.centroids(externalCells, 2) > (ztop - triSize));
    %plotGrid(G); hold on; plotGrid(G, topCells, 'facecolor', 'b')

    % Find triangles sharing a face with these (the ones to modify)
    f2cn = gridCellNo(G);
    idFacesTopCells = find(ismember(f2cn, topCells));
    facesTopCells = G.cells.faces(idFacesTopCells);
    %plotGrid(G); hold on; plotFaces(G, facesTopCells, 'edgecolor', 'r')
    cellsTheseFaces = G.faces.neighbors(facesTopCells, :);
    lowerTopCells = cellsTheseFaces(~ismember(cellsTheseFaces, topCells));
    lowerTopCells = unique(lowerTopCells(lowerTopCells~=0));
    %plotGrid(G); hold on; plotGrid(G, lowerTopCells, 'facecolor', 'r')

    % Loop over each one of these, get centroid of shallow neighbors (inverted)
    % and assign midway.
    for n=1:numel(lowerTopCells)
       cid = lowerTopCells(n);
       fid = G.cells.faces(ismember(f2cn, cid));
       neighbors = G.faces.neighbors(fid, :); % all, including cid and 0
       neighbors(neighbors == 0) = [];
       neighbors = unique(reshape(neighbors, numel(neighbors), 1));
       neighbors = neighbors(G.cells.centroids(neighbors, 2) > ...
                             G.cells.centroids(cid, 2));
       %plotGrid(G); hold on; plotGrid(G, neighbors, 'facecolor', 'r')
       G.cells.centroids(cid, 2) = mean(G.cells.centroids(neighbors, 2));
    end
    
    % Plot Centroid locations
    plotGrid(G, 'edgeColor', 'none'); hold on
    cellsTheseFaces(cellsTheseFaces == 0) = [];
    text(G.cells.centroids(cellsTheseFaces, 1), ...
        G.cells.centroids(cellsTheseFaces, 2), 'o')
                        
end

% Rock
rock.poro = 0.3*ones(G.cells.num, 1);
rock.perm = 200*milli*darcy*ones(G.cells.num, 1);
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.rocknum    = ones(G.cells.num, 1);

% Fluid 
if krHyst == true 
    fn = 'simple_hyst.DATA';
else
    fn = 'josimarGoM.DATA';
end
deck  = convertDeckUnits(readEclipseDeck(fn));
cP    = deck.PROPS.ROCK(:, 2);
fluid = initDeckADIFluid(deck);

if isfield(fluid, 'krHyst') && fluid.krHyst == 1
   assert(strcmp(fn, 'simple_hyst.DATA')) % otherwise change below
   numreg = max(rock.regions.saturation);
   fluid.krHyst = numreg + 1;            
   rock.regions.imbibition = rock.regions.saturation + numreg;
   fluid = addScanKr(fluid, rock.regions.imbibition);
end

% Initialize
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
d  = gravity() ./ g;
z_rot = G.cells.centroids * d(1:2).';
p_r = g*rho_wr*min(abs(G.faces.centroids(:,2)));  % add ref water column

[z_0, z_max] = deal(min(z_rot), max(z_rot));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, z_rot), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

figure(1); subplot(1,2,1)
plotCellData(G, state0.pressure, 'LineWidth', 0.1, 'edgeColor', 'none')
view([-5, 90]); axis equal off; colorbar; 

% Wells & timesteps
dist = pdist2(G.cells.centroids, [500 -2750]);
[~, perforatedCell] = min(dist,[],1);

injTime = 30*year;
mrate = 1e-3*10^9/year;             % Inject 1 Mt/y per km in third dim
injmass = mrate*injTime;
rhoInj = fluid.rhoGS;
injrate = mrate/rhoInj;                                      % [Sm^3/s]
W = addWell([ ], G, rock, perforatedCell, 'Name', 'I1', 'Dir', 'y', ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', z_rot(perforatedCell), 'radius', 0.3); 

simTime     = 300*year;        
reportTimes = [[1, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
               240, 300, 365, 456.25, 547.5, 638.75, 730]*day, ...
               [2.2:.2:32, 32.5:.5:35 36:150 155:5:300]*year];
timesteps   = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

%% MODEL SETUP and BC
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
model.minimumPressure = min(state0.pressure);

% Mex and Solver
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'deferredAssembly', true);
model = model.validateModel();
nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, 'targetChangeAbs', 0.25);
nls.timeStepSelector = stepSel;
nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');

% Changes to model. Note that this must be done after setting the autodiff
% backend (in setAcceleration_bo) since a copy of the backend will be made
% when setting the flow property functions.
model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
model.PVTPropertyFunctions  = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));

% Boundary Conditions (Grid Dependent)
bc = [];
lim = [valminBC max(G.faces.centroids(:, 1))-valminBC];
east = G.cells.centroids(:,1) < lim(1);
west = G.cells.centroids(:,1) > lim(2);
vert = G.cells.centroids(:,2) < zlim;
cellsExt = any([east, all([west, vert], 2)], 2);
model.operators.pv(cellsExt) = model.operators.pv(cellsExt)*1e4;

figure(1); subplot(1,2,2)
hold on; plotGrid(G, 'edgeColor', 'none'); 
plotGrid(G, find(cellsExt), 'facecolor', 'r', 'edgeColor', 'none')
view([-5, 90]); axis equal off;

% Schedule
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules 
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % nr of wells must be the same for the entire simulation
schedule.control(2).W.val = 0;                                              % field 2 rate 0 (after injection)
schedule.step.control(cumsum(schedule.step.val) > injTime) = 2; 


%% SIMULATION
% we run a parallel simulation with N threads
%pause
N = 4;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                   % Specify threads manually
%schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);                      % Simple schedule, this injects for the total simTime
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                                'NonLinearSolver', nls, ...
%                                                'Verbose', true); 
problem = packSimulationProblem(state0, model, schedule, folderName, ...
                                'Name', folderName, 'Directory', topDir, ...
                                'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem); 


%% RESULTS PLOTS
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
end

% Save data
% save([folderName '.mat'])

fh = figure(5);
plotToolbar(G, states); axis equal off; view([-5 90])
set(fh, 'position', [500, 200, 600, 450]);

% pressure and saturation eoi, eos
figure(6)
plotCellData(G,states{79}.dp, 'EdgeColor', 'none')
c = colorbar; caxis([0 600000])
axis equal off; view([-5 90])

figure(7)
plotCellData(G,states{end}.dp, 'EdgeColor', 'none')
caxis([0 600000]); axis equal off; view([-5 90])

figure(8)
plotCellData(G,states{79}.s(:,2), 'EdgeColor', 'none')
c = colorbar; caxis([0 0.7]);
axis equal off; view([-5 90])

figure(9)
plotCellData(G,states{end}.s(:,2), 'EdgeColor', 'none')
caxis([0 0.7]); axis equal off; view([-5 90])
