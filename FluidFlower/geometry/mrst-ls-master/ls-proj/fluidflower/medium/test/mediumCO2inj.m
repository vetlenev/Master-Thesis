%% Perform CO2 injection (blackoil module)
%--------------------------------------------------------------------------
% 3 saturation regions, kr hysteresis, fault Pc scaled
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
mrstVerbose on


%% Options
hyster = 'on';
vapoil = 'off';
disgas = true;
mesh = {500, 'composite', 'mediumG_composite_nc500.mat'};      % 150 (coarse), 250 (medium), 500 (fine)
wellno = 1;                                                    % n injectors
rate = 1;                                                      % mL/min should be ok up to a few mL/min
topboundary = 'straight';
bctype = 'cp';                                                % constant pressure (cp) or pore volume multiplier (pvm)
permtype = 'logn';
sigmap = 'varied';


%% Save data
folderName = ['mediumMesh', num2str(mesh{1}), '_', mesh{2}, '_bc_', bctype, ...
              '_perm_', permtype, '_sigmap_', num2str(sigmap), '_hyster_', hyster, ...
              '_vapoil_', vapoil, '_rate_mm_', num2str(rate), ...
              '_Ninj_', num2str(wellno)];
topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\tests\';

%% Plots
fig3D = @() figure('Position', [0, 0, 1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );


%% Mesh
% Load simpleExtrudedMesh
pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/test/');
G_dat  = load(fullfile(pth, mesh{3}));
sealID = [3 6 8 11 13];
uGr = {[1 7], [3 8], [4 9], [5 10], [6 11], 12, 13, 14, 2};
G = G_dat.G;
faultID = 2;

% Plot
% fig3D();
% %%colormap(jet); plotCellData(G, G_dat.p); view([90, 0]); 
% Ncompart = max(G_dat.p);
% cmap = copper(Ncompart);
% colr = [10, 2, 6, 10, 2, 10, 2, 8, 6];
% for n=1:numel(uGr)
%     plotCellData(G, colr(n)*ones(sum(ismember(G_dat.p, uGr{n})), 1), ...
%                  ismember(G_dat.p, uGr{n}), 'edgealpha', 0.2)
% end
% outlineCoarseGrid(G, G_dat.compartID,'EdgeColor','w','LineWidth',2);
% setAxProps(gca), %camlight();
% colormap(copper); %c = colorbar; set(c, 'YTick', sort(colr));
% axis equal off, view([90, 0])
% ylim([0.43 0.46]), zlim([0.43 0.47])


%% Rock
% Layers follow uGr (bot to top, left to right for faulted, fault at the end)
rock.poro = nan(G.cells.num, 1);
rock.perm = nan(G.cells.num, 1);    % isotropic
poroMean = [0.45, 0.38, 0.4, 0.45, 0.38, 0.45, 0.38, 0.45, 0.4];
permMean = [1500, 30, 200, 600, 30, 600, 30, 500, 100];   % [D]
sigmas = [40, 20, 10, 40, 20, 40, 20, 20, 20];
Ly = 0.89;
Lz = 0.47;
Ny = 1000;
Nz = 500;
rng(5931)
for n=1:numel(uGr)
    % Cell ids of corresponding unit
    cid = ismember(G_dat.compartID, uGr{n});
    nc = sum(cid);
    
    % Porosity
    r = 0.03; % half range
    rock.poro(cid) = (poroMean(n) - r) + rand(nc, 1)*(2*r);
    
    % Permeability
    if strcmp(permtype, 'randm')    % random
        logr = 0.2;
        rock.perm(cid) = 10.^(rand(nc, 1)*(2*logr) + ...
                              (log10(permMean(n)) - logr));
    elseif strcmp(permtype, 'logn') % lognormal
        stdev = sigmap*permMean(n);
        var_lnk = log(sigmas(n)*2);
        corr_leny= 0.1*Ly;
        corr_lenz= 0.1*Lz;
        [perm, var_lnk_actual]= random_perm(var_lnk, corr_leny, corr_lenz, ...
                                            Ny, Nz, Ly, Lz);
        perm = perm + (permMean(n) - mean(perm(:)));
        perm = reshape(perm, 1, Ny, Nz);
        rock.perm(cid) = sampleFromBox(G, perm, find(cid));
%         N = [2, 200, 100];
%         if ismember(faultID, uGr{n})
%            N = [200, 100, 2]; 
%         end
%         K = reshape(logNormLayers(N, [permMean(n) permMean(n)], ...
%                     'sigma', sigma*1.6, 'std', 1000), N);
% %         if sigma == 5
% %             q = permMean(n)/6;
% %         elseif sigma == 3
% %             q = permMean(n)/7;
% %         elseif sigma < 3
% %             q = permMean(n)/7.3;
% %         end
% %         rock.perm(cid) = sampleFromBox(G, q.*K, find(cid));
%         rock.perm(cid) = sampleFromBox(G, K, find(cid));
    end
end
rock.perm = rock.perm*darcy; % SI units

% Plot permeability
% cid = ismember(G_dat.compartID, uGr{1});
% fig3D(); plotCellData(G, log10(rock.perm/(darcy)), 'edgealpha', 0.2)
% setAxProps(gca), colormap(turbo), c = colorbar; %caxis([3 4])
% axis equal off; view([90 0]), %ylim([0.3 0.6]); zlim([0.42 0.47])


%% Fluids
% 3 Saturation regions
fault.satRegNum = 3;
rock.regions.saturation = ones(G.cells.num, 1);                              % reservoir units
rock.regions.saturation(ismember(G_dat.compartID, sealID)) = 2;
rock.regions.saturation(ismember(G_dat.compartID, faultID)) = fault.satRegNum;

% 1 Regions for rock compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% load deck and initialize fluid
pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/');
fn  = fullfile(pth, 'krHyst_fPcFault_co2wat_test_PVDG_05PVTO.DATA');
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);

% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Scaled Pc in fault
refPc.val  = fluid.pcOG{fault.satRegNum};
refPc.poro = 0.37; %0.07;                             % Table 3.3, sample 5 in Treviño & Meckel, 2017 (Pc curve is the green one in Fig. 3.4a)
refPc.perm = 200.8*milli*darcy; %0.00208*milli*darcy; % "
fault.poro = rock.poro(G_dat.compartID==faultID);
fault.k.vals = rock.perm(G_dat.compartID==faultID);
fluid = assignFaultPc(fluid, fault, refPc);

% Hysteresis
if isfield(fluid, 'krHyst') && fluid.krHyst == 1
   numReg = max(rock.regions.saturation);
   fluid.krHyst = 1 + numReg;                % porous rock is first imb region (4), fault is third (6)
   rock.regions.imbibition = rock.regions.saturation + numReg;
   fluid = addScanKr(fluid, rock.regions.imbibition);
end

% Plot density, viscosity
% f = fluid;
% p = linspace(1, 5, 50)*barsa;
% rhoG = f.bG(p).*f.rhoGS;
% muG  = f.muG(p);
% rs   = f.rsSat(p);
% rhoB = f.bO(p, rs, true(numel(p), 1)).*(rs.*f.rhoGS + f.rhoOS);
% rhoB_noCO2 = f.bO(p', zeros(numel(p), 1), false(numel(p), 1)).*f.rhoOS;
% muB  = f.muO(p, rs, true(numel(p), 1));
% plot(p, rhoB, '-b'); hold on; plot(p, rhoB_noCO2, '--b')
% 
% p2 = 4*barsa*ones(50,1);
% rs = linspace(0,f.rsSat(p2(1)),50)';
% rhoB2 = [f.bO(p2(1:end-1),rs(1:end-1),false(49,1)); ...
%          f.bO(p2(1),rs(end),true)].*(rs.*f.rhoGS + f.rhoOS);
% plot(rs, rhoB2, '-k'); xlabel('rs'), ylabel('\rho [kg/m^3]')

%% Initialize
gravity reset on
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
water_column = 0.3; % m
p_r = 1*barsa + g*rho_wr*water_column; % p at shallowest z
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
if ~disgas   
    equil = ode23(@(z,p) g .* fluid.bO(p)*fluid.rhoOS, [z_0, z_max], p_r);
else
    equil = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
end
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = zeros(G.cells.num, 1);             % no vaporized water at the beginning
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

%plot
% fig3D(); plotCellData(G, p0/barsa, 'edgealpha', 0.2)
% setAxProps(gca), colormap(jet), c = colorbar;
% axis equal off


%% Wells
%t = [5*minute [0.5, 14]*day];
t = [6*minute  [1, 48]*hour];     % rampup, injtime, simtime                          
injrate = rate*(milli*litre)/(minute*wellno);                               % [Sm^3/s]
injvol = injrate*t(2);
rhoInj = fluid.rhoGS;
mrate = injvol*rhoInj/(t(2)*wellno);                                        
wellInx = find(G_dat.wellNo==6);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', G.cells.centroids(wellInx, G.griddim), ...
            'Radius', 0.9*1e-3);                                % refDepth must be included for 'bhp' wells.
% reportTimes = [(1:t(1)/minute)*minute, ...                      % rampup inj rate
%                [10:10:120 150:30:720]*minute, ...               % injection
%                (722:1:740)*minute, ...                          % rampdown inj rate
%                [12.5:.5:48 49:72 76:2:96 100:4:120]*hour, ...   % rest
%                (5.5:.25:t(end)/day)*day];                       % rest
% reportTimes = [(1:t(1)/minute)*minute, ...                        % rampup inj rate
%                [10:5:120 130:10:720]*minute, ...                  % injection
%                (722:1:740)*minute, ...                            % rampdown inj rate
%                (750:10:7200)*minute];                             % rest
% reportTimes = [(0.5:0.5:t(1)/minute)*minute, ...                  % rampup inj rate
%                [6:1:360]*minute, ...                              % injection
%                (360.5:.5:370)*minute, ...                         % rampdown inj rate
%                (371:1:2880)*minute];                              % rest
reportTimes = [1:240 245:5:2880]*minute;                          % rest
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==t(end), 'sum of timesteps must equal simTime')


%% Model
% If you get the badoil error, delete a few timesteps before and decrease
% the toleranceCNV. Another option could be to decrease the maxIterations
% below the itnum at which you get the badoil error.
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'vapoil', false, 'water', false);
                         
% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true);
model.toleranceCNV = 1e-3;
%model.toleranceMB = 1e-5;
model = model.validateModel();

% nls = getNonLinearSolver(model, 'TimestepStrategy', 'ds', ...
%                          'useCPR', true);
nls = getNonLinearSolver(model, 'TimestepStrategy', 'none', ...
                         'useCPR', true);
stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, ...
                                     'targetChangeAbs', 0.25);
nls.timeStepSelector = stepSel;
if strcmp(bctype, 'pvm')
    nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'bicgstab');
end
nls.LinearSolver.maxIterations = 50;
%nls.LinearSolver = AMGCL_CPRSolverBlockAD('tolerance', 1e-4, 'Solver', 'gmres', ...
%                                          'preconditioner', 'amg', ...
%                                          'coarsening', 'smoothed_aggregation');
nls.useRelaxation = 1;
%nls.relaxationType = 'sor';
%nls.enforceResidualDecrease = 1;
%nls.useLinesearch = 1;
nls.maxIterations = 6;
nls.maxTimestepCuts = 12;
nls.acceptanceFactor = 5;
%nls.alwaysUseStabilization=1;
%nls.convergenceIssues=true;
%nls = [];

% Model changes (after acceleration!)
model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
if strcmp(hyster, 'on')
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
    %model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));
end
model.minimumPressure = min(state0.pressure);

%% BCs
if strcmp(topboundary, 'straight')
    f = false(G.faces.num, 1);
    f(boundaryFaces(G)) = true;
    f = all([G.faces.centroids(:, 3) < 1e-4, f], 2);
    if strcmp(bctype, 'cp')
        bc = addBC([], find(f), 'pressure', p_r, 'sat', [1 0]);
    elseif strcmp(bctype, 'pvm')    
        cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
        cellsext(cellsext==0) = [];
        %     fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
        %     setAxProps(gca), axis equal off, view([90 0]);
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e5;
        bc = [];
    end
else
    % Find top faces
    f = false(G.faces.num, 1);
    f(boundaryFaces(G)) = true;
    zmean = mean(G.cells.centroids(:, 3));
    f = all([G.faces.centroids(:, 3) < zmean, f], 2);
    f = all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 2))], 2);
    f = find(all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 1))], 2));

%     fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r')
%     setAxProps(gca), colormap(jet), c = colorbar;
%     axis equal off

    if strcmp(bctype, 'cp')
        [z_0, z_max] = deal(min(G.faces.centroids(f,3)), max(G.faces.centroids(f,3)));
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
        fp_val = reshape(deval(equil, G.faces.centroids(f,3)), [], 1);  clear equil
        bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0]);
    elseif strcmp(bctype, 'pvm')
        cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
        cellsext(cellsext==0) = [];
        %fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
        %setAxProps(gca), axis equal off, view([90 0])
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e4; 
        bc = [];
    end
end

%% Schedule
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
v = injrate;
% injrates = [0.005*v 0.05*v 0.1*v 0.25*v 0.75*v v 0.75*v 0.5*v 0.25*v ...
%             0.1*v 0.05*v 0];
injrates = [0.05*v 0.1*v 0.25*v 0.5*v 0.75*v v 0.75*v 0.5*v 0.25*v ...
            0.1*v 0.05*v 0];
tmp = cell(numel(injrates), 1);                                             % create 2 schedules
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
for n=1:numel(injrates)
    schedule.control(n).W = W;                                              % field 1 used during injection
    schedule.control(n).W.val = injrates(n);                                % nr of wells must be the same for the entire simulation
    schedule.control(n).bc = bc;
end
% Smaller than rampup
idStep = find(cumsum(schedule.step.val) < t(1));
schedule.step.control(1:5) = 1:5;
schedule.step.control(6:idStep(end)) = 6;
schedule.step.control(idStep(end)+1:end) = 6;

% Larger than injtime
idStep2 = find(cumsum(schedule.step.val) > t(2), 1);
schedule.step.control(idStep2:idStep2+numel(injrates)-2-5) = ...
    max(schedule.step.control)+1:numel(injrates); % set timesteps after injTime to use Well field 2
schedule.step.control(find(schedule.step.control > max(idStep)+1, 1, 'last')+1:end) = ...
                                                           numel(injrates); 

%% Simulation
N = 4;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                           'NonLinearSolver', nls, ...
%                                           'Verbose', true);
 
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  


%% Results
% compute quantitites
model = model.validateModel();
bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;
    rho = model.getProp(states{n}, 'ComponentPhaseDensity');
    states{n}.FlowProps.ComponentPhaseDensity = rho{:,1}; % 2 components in brine phase.
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    bhp(n) = wellSols{n}.bhp;
end
states{1}.reg = model.rock.regions.saturation;


% basic overview
%cmap = flipud(cmocean('tempo')); 
cmap = hot; 
seal1 = find(ismember(G_dat.compartID, uGr{2}));
seal2 = find(ismember(G_dat.compartID, uGr{5}));
fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); 
plotGrid(G, seal1, 'edgecolor', [0.5 0.5 0.5]);
plotGrid(G, seal2, 'edgecolor', [0.5 0.5 0.5]);
setAxProps(gca), colormap(cmap), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), ylim([0.4 0.5]), zlim([0.40 0.47])

% Find cell and check kr (krHyst model ok?)
% tid = 61;
% dist = pdist2(G.cells.centroids(:,2:3), [0.8 0.411]);
% [~, inx] = min(dist,[],1);
% %inx = wellInx;
% krout = states{tid}.FlowProps.RelativePermeability{2}(inx);
% 
% % fig3D();
% % plotGrid(G, 'facecolor', 'none', 'edgecolor', [0.5 0.5 0.5], 'edgealpha', 0.2);
% % hold on, plotGrid(G, inx, 'facecolor', 'r');
% % setAxProps(gca), axis equal, view([90 0])
% 
% sg = states{tid}.s(inx, 2);
% sgmax = states{tid}.sMax(inx, 2);
% kr = fluid.krGi{1}(sg, sgmax);
% disp(['sG, sGmax: ' num2str(sg), ', ' num2str(sgmax)]);
% disp(['Output kr: ' num2str(krout)]);
% disp(['Model kr: ' num2str(kr)]);