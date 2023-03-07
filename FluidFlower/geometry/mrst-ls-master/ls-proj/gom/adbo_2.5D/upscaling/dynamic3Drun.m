function [states, G2, rock2, fluid, state0, rate] = dynamic3Drun(G, rock , fluid, ...
                                                                 trep, fault3D, opt, states_plot)
%
%
%

% 2. Prepare fluid structure
if opt.dyn_ispc == 0
    fluid.pcOGu = fluid.pcOG;
    nri = numel(fluid.krG);
    fluid.pcOG = cell(1, nri);
    for k=1:nri
        fluid.pcOG{k} = @(p) zeros(numel(p), 1);
    end
elseif opt.dyn_ispc == 1
    % We need to have same nreg for kr as for pc
    krg = fluid.krG;
    fluid.krG = cell(1, numel(fluid.isclay));
    fluid.krG(fluid.isclay) = krg(2);
    fluid.krG(~fluid.isclay) = krg(1);
    
    % No empty kr, pc indices (otherwise mrst sim won't work).
    pcog = fluid.pcOG;
    idu = find(~cellfun(@isempty, pcog));
    if any(diff(idu) > 1)
        fluid.pcOG = cell(1, numel(idu));
        fluid.pcOGu = pcog;
        fluid.pcOG(1:numel(idu)) = pcog(idu);
        
        % And adjust rock.regions.saturation
        idv = 1:numel(idu);
        for k=1:numel(idu)
            if idu(k) ~= idv(k)
                rock.regions.saturation(rock.regions.saturation == idu(k)) = idv(k);
            end
        end
    end
end

% 3. Get new grid adding a boundary layer of sand-perm cells, and
%    modifiy depth for correct p initialization and varying density and
%    viscosity of fluids.
physDim = max(G.faces.centroids);
physDim(end) = physDim(end) + G.cellDim(end);
cellDim = G.cartDims;
cellDim(end) = cellDim(end) + 1;
G2 = cartGrid(cellDim, physDim);
%    layer center from mesh - half layer thickness, approx.
zmin = opt.zmax{1}(end) - round(opt.thick{1}(end)/2);
G2.nodes.coords(:,end) = G2.nodes.coords(:,end) + zmin;
G2 = computeGeometry(G2);

% 4. Get rock, with boundary layer of sand-perm cells.
idG = reshape(1:G.cells.num, G.cartDims(1)*G.cartDims(2), G.cartDims(3));
idG = reshape(fliplr(idG), G.cells.num, 1);
rock2.poro = nan(G2.cells.num, 1);
rock2.poro(1:G.cells.num) = rock.poro(idG);
rock2.poro(G.cells.num+1:end) = mean(rock.poro(~fault3D.Grid.isSmear));

rock2.perm = nan(G2.cells.num, 6);
rock2.perm(1:G.cells.num, :) = rock.perm(idG, :);
rock2.perm(G.cells.num+1:end, :) = repmat(mean(rock.perm(~fault3D.Grid.isSmear, :)), ...
                                                G.cartDims(1)*G.cartDims(2), 1);

rock2.regions.rocknum = ones(G2.cells.num,1);
if ~opt.dyn_incomp_run
    rock2.regions.rocknum(1:G.cells.num) = rock.regions.rocknum(idG);
end

% Only 2 regions for relperm, sand and clay smear.
regsat = ones(G.cells.num, 1);
regsat(fault3D.Grid.isSmear) = 2;
rock2.regions.saturation = ones(G2.cells.num, 1);
rock2.regions.saturation(1:G.cells.num) = regsat(idG);
%rock2.regions.saturation = nan(G2.cells.num, 1);
%rock2.regions.saturation(1:G.cells.num) = rock.regions.saturation(idG);
%uCellsBot = min(fault3D.Grid.units(~fault3D.Grid.isSmear));
%rock2.regions.saturation(G.cells.num+1:end) = uCellsBot;

switch opt.dyn_perm_case
    case 'sand'
        rock2.perm = mean(rock2.perm(rock2.regions.saturation==1,end)) ...
                     *ones(G2.cells.num,1);
        %rock2.perm = 1e-20*ones(G2.cells.num, 1);
        rock2.poro = mean(rock2.poro(rock2.regions.saturation==1)) ...
                    *ones(G2.cells.num,1);
    case 'geomean'
        rock2.perm = geomean(rock2.perm(:,end))*ones(G2.cells.num,1);
        rock2.poro = mean(rock2.poro)*ones(G2.cells.num,1);
end

% check plot
%     figure(randi(1000, 1)); subplot(1,2,1)
%     plotCellData(G, log10(rock.perm(:,1)/(milli*darcy)), 'edgecolor', 'none');
%     view([30 20]); ax=gca; ax.ZDir = 'normal'; ax.DataAspectRatio = [0.05 1 1];
%     colormap(copper); colorbar;
%     subplot(1,2,2)
%     plotCellData(G2, log10(rock2.perm(:,1)/(milli*darcy)), 'edgecolor', 'none');
%     view([30 20]); ax=gca; ax.DataAspectRatio = [0.05 1 1];
%     colormap(copper); colorbar;
% 
%     figure(randi(1000, 1)); subplot(1,2,1)
%     plotCellData(G, rock.regions.saturation, 'edgecolor', 'none');
%     view([30 20]); ax=gca; ax.ZDir = 'normal'; ax.DataAspectRatio = [0.05 1 1];
%     colormap(flipud(copper)); colorbar;
%     subplot(1,2,2)
%     plotCellData(G2, rock2.regions.saturation, 'edgecolor', 'none');
%     view([30 20]); ax=gca; ax.DataAspectRatio = [0.05 1 1];
%     colormap(flipud(copper)); colorbar;

% 5. Initialize
g = 9.807;         % m/s^2
water_column = 50; % m
if opt.zmax{1}(1) > 1000 && opt.zmax{1}(1) < 3000
    h_wat = water_column + min(G2.faces.centroids(:,3));
    pO = fluid.rhoOS*g*h_wat;
    p_r = g*fluid.bO(pO)*fluid.rhoOS*h_wat;
end
[z_0, z_max] = deal(min(G2.cells.centroids(:,3)), max(G2.cells.centroids(:,3)));
equil  = ode23(@(z,p) g*fluid.bO(p)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G2.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G2.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = 0;                                  % dead oil (brine)
rv0 = 0;                                  % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% Incompressible?
if opt.dyn_incomp_run   % we leave a bit of compressibility
    %fluid.pvMultR{1} = @(p) ones(numel(p), 1);
    fluid.pvMultR = fluid.pvMultR{1};
    %fluid.bG = @(x) fluid.bG(mean(p0))*ones(numel(x),1);
    fluid.muG = @(x) fluid.muG(mean(p0))*ones(numel(x),1);
    %fluid.bO = @(x) fluid.bO(mean(p0))*ones(numel(x),1);
    fluid.muO = @(x) fluid.muO(mean(p0))*ones(numel(x),1);
end

% plot
%     h = figure(40);
%     colormap(turbo);
%     plotToolbar(G2, p0/barsa), set(gca,'Xdir','reverse');
%     view(30,20); ax = gca; ax.DataAspectRatio = [0.05 1 1];
%     title('Initial pressure')
%     c = colorbar;
%     c.Label.String = '$p_0$ [bar]';  c.Label.FontSize = 20;
%     c.Label.Interpreter = 'latex';   c.FontSize = 16;
%     set(h, 'Position', [100, 100, 600, 600])

% 6. Model and Acceleration
model = GenericBlackOilModel(G2, rock2, fluid, 'disgas', false, ...
                             'vapoil', false, 'water', false);

% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
    'deferredAssembly', true);
model.toleranceCNV = 1e-3;
model.toleranceMB = 2e-7;
model = model.validateModel();
nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
    'useCPR', true);
nls.LinearSolver.maxIterations = 50;
nls.useRelaxation = false;
nls.useLinesearch = true;
nls.maxIterations = 10;
nls.maxTimestepCuts = 14;

% Model changes (after acceleration!)
model.operators.p0 = state0.pressure;                            % Needed for MyPvMult
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
model.minimumPressure = min(state0.pressure);

% 7. BCs
pvt = sum(model.operators.pv);
if opt.zmax{1}(1) > 1000 && opt.zmax{1}(1) < 3000
    rate = opt.dyn_mrate*pvt/day();
end
bc = fluxside([],  G2, 'Bottom', rate, 'sat',  [0 1]);
bc = pside(bc, G2, 'Top', min(p0), 'sat', [1 0]);

% 8. Schedule
tsim = opt.dyn_tsim_year*year;
%     trep = [60*minute (4:4:24)*hour [5:5:30 60:30:360]*day 1*year ...
%             (1.5:.5:10)*year [11:1:100 102:2:500 505:5:1000]*year ...
%             [1100:100:10000]*year];
%     trep = [60*minute (4:4:24)*hour [5:5:30 60:30:360]*day ...
%             ([1:.2:5 6:100])*year];
timesteps = [trep(1) diff(trep)];
assert(sum(timesteps)== tsim, 'sum of timesteps must equal simTime')
schedule = simpleSchedule(timesteps, 'W', [], 'bc', bc);

% 9. Simulation
N = 2;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', ...
                                 nls, 'Verbose', true);
%states = [];
%ds = norm(states{end-1}.s(:,1) - states{end}.s(:,1), inf);
%assert(ds < 1e-3)

% Check plot
model = model.validateModel();
for k=1:numel(states)
    states{k}.dp = states{k}.pressure - state0.pressure;
    states{k}.FlowProps.RelativePermeability = model.getProp(states{k}, 'RelativePermeability');
end
states{1}.reg = model.rock.regions.saturation;

if states_plot
    plotToolbar(G2, states, 'edgealpha', 0.2); hold on
    cmap = turbo;
    colormap(cmap), c = colorbar; %clim([0 40000])
    %axis off
    view([30 20]), hold off
    set(gca, 'ColorScale', 'linear')
    caxis([0 1])
    ax = gca;
    ax.DataAspectRatio = [0.1 1 1];
end
    
        
end