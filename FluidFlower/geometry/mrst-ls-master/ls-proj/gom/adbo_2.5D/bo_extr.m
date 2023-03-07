%% Base case model (2D extruded) for GoM CO2 fault migration
%
%____________________________________________________________________
clc, clear, close all force;
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils
       
mrstVerbose on

%% CHOOSE SIMULATION OPTIONS
%fname     = 'sc1_center_sc2Mesh63x23523_SGRauto_Spe02_80C_krHyst_50y500y';
%fname     = 'sc1_NextToF_FW_noPcFault_sc2Mesh63x23523_SGRauto_Spe02_80C_krHyst_50y500y';
%fname     = 'sc1_NextToF_sc2Mesh63x23523_SGRauto_Spe02_80C_krHyst_50y500y';
%topDir    = 'C:\Users\lsalo\matlab\sim_data\mrst\gomBasinFaultLeak\sc1\nextToF';
%topDir     = 'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc1\';
%fname = 'sc2_baseCase_nextToF_sc2Mesh63x23523_SGRauto_Spe02_krHyst_pvtBase_50y500y';
fname = 'sc2_baseCase_nextToF_noPcFault_Fk175mD_sc2Mesh63x23523_SGRauto_Spe02_krHyst_pvtBase_50y500y';
topDir     = 'C:\Users\lsalo\matlab\sim_data\mrst\gomBasinFaultLeak\sc2\baseCase_nextToF';
%fname = 'sc2_predict_nextToF_sc2Mesh63x23523_highPerm_noHyst_pvtBase_30y200y';
%topDir     = 'C:\Users\lsalo\matlab\sim_data\mrst\gomBasinFaultLeak\sc2';
%topDir     = 'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc2\baseCase_nextToF';
figs      = [0, 0, 0, 0];                
scenario  = 2;
%studyCase = 3;  % for sc 1; 2 = nextToF (HW), 3 = nextToF (FW), 4 = domain center
studyCase = 'sc2_baseCase';
%studyCase = 'sc2_predict_theta30_nohyst';
[mesh, opt, wells, t, resPlots] = simOpts(scenario, studyCase);
opt.permCase = 'high'; % for sc2 with predict only


%% GET SIMULATION VARIABLES
% Domain grid & set gravity
[G, cellw] = getMesh_bo(mesh, figs(1));

% Simulation grid, Rock, unit cell ids, fault props
[G, rock, ucids, fault] = getRockParams_bo(mesh, G, opt, figs(2));

% Fluid & Regions
[fluid, rock, cP, fault, refPc] = getFluid_bo(mesh, G, rock, opt, ...
                                              ucids.unit_cell_ids, fault, figs(3));

% Initialize
state0 = init_bo(G, fluid, opt.watCol, ucids.unit_cell_ids, figs(4));

% Wells & timesteps
[W, timesteps, wellInx] = getWell_bo(G, rock, fluid, wells, t);


%% MODEL SETUP and BC
% modify pvMult to take into account:
%   - In rock compressibility: initial p at each cell
% modify fluid.Pc to take into account:
%   - In capillary pressure: pc is a fcn of poro and perm in the FZ (selected cases only)
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);
if strcmp(opt.faultPc, 'scaled')
    [fluid] = assignFaultPc(fluid, fault, refPc);
end

% TEST ONLY to force fault migration (ucids 61 mesh dependent!)
disp('****************************************************************')
disp('---------------  ATTENTION: PC FAULT SET TO 0!  ----------------')
fluid.pcOG{3} = @(sg) zeros(numel(sg), 1);                     % pc=0, only for test!
rock.perm(ucids.unit_cell_ids{61},[1 4 6]) = 175*milli*darcy;   % kxx, kyy, kzz = 10 mD
rock.perm(ucids.unit_cell_ids{61},[2 3 5]) = 0;                % kxy, kxz, kyz = 0
disp('****************************************************************')

% Define base model class
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, ...
                             'water', false);

% Mex and Solver
[model, nls] = setAcceleration_bo(model, opt);

% Changes to model. Note that this must be done after setting the autodiff
% backend (in setAcceleration_bo) since a copy of the backend will be made
% when setting the flow property functions.
model.operators.p0          = state0.pressure;                            % Needed for MyPvMult
if isfield(fluid, 'krHyst')
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
end
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Density', BlackOilDensity(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Mobility', Mobility(model));
%model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhasePressures', MyPhasePressures(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PressureReductionFactors', BlackOilPressureReductionFactors(model));
% model.FlowPropertyFunctions = model.FlowP ropertyFunctions.setStateFunction('ShrinkageFactors', BlackOilShrinkageFactors(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('Viscosity', BlackOilViscosity(model));
% model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RsMax', RsMax(model));
%model.FlowPropertyFunctions = MyFlowPropertyFunctions(model);
model.minimumPressure = min(state0.pressure);

% Boundary Conditions
[bc, model] = getBC_bo(G, fluid, model, opt, mesh, ucids);
schedule = getSchedule_bo(timesteps, W, bc, t);

%% SIMULATION
% Run a parallel simulation with N threads
N = 8;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
%nls.LinearSolver.tolerance = 1e-3;
%nls.LinearSolver.setSolver('bicgstab')
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                            'NonLinearSolver', nls, ...
%                                            'Verbose', true);
 
problem = packSimulationProblem(state0, model, schedule, fname, 'Name', fname, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);                                            

%% RESULTS, PLOTS AND SAVE DATA
model = model.validateModel();
bhp = zeros(1, numel(states));
c_refr =  ucids.unit_cell_ids{35}(~ismember(ucids.unit_cell_ids{35}, ...
                                            ucids.unit_cell_ids{36})); %grid sc2!
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    pc = model.getProp(states{n}, 'CapillaryPressure');
    states{n}.FlowProps.CapillaryPressure = pc{1,2};
    states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    bhp(n) = wellSols{n}.bhp;
    dp_wellc(n) = states{n}.dp(W.cells);
    dp_resr(n) = mean(states{n}.dp(c_refr));
    dp_resl(n) = mean(states{n}.dp(ucids.unit_cell_ids{36})); % grid sc2!
end
states{1}.reg = model.rock.regions.saturation;

% Save data
if ~isempty(fname)
    disp(['ATTENTION: Problem and model data will be saved in: ' pwd])
    %disp('Press Ctrl+C to CANCEL or any other key to SAVE data.')
    
    save(['proMod_' fname '.mat'], 'problem', 'model', 'W', 'ucids', ...
         'state0', 'ncont', 'mesh', 'timesteps', 'fault', '-v7.3')
    %save([fname '_states.mat'], 'states', '-v7.3') % larger than 2GB
end

% Mesh Plots
if opt.sc == 1 && studyCase == 2
    injCase = 'nextToF';
elseif opt.sc == 1 && studyCase == 3
    injCase = 'nextToF_FW';
else
    injCase = 'center';
end
%ncont = load('ls-proj/gom/adbo_2.5D/data_files/ncont_29lyr.mat');
%ncont = ncont.ncont
plotMeshRes_bo(model, ucids, W, states, state0, ncont, mesh, resPlots, opt, injCase) 

% Mesh video
close all
%latx = {'interpreter','latex'};
%dir = 'C:/Users/lsalo/matlab/vids/gomFaultLeakage/sc1/';
dir = 'C:/Users/lsalo/matlab/vids/gomFaultLeakage/sc2/baseCase_S/';
idp = numel(timesteps);
generateFigs = 1;
v = VideoWriter([dir 'sc2S_FW_5fps.mp4'], 'MPEG-4');
v.Quality = 100;
v.FrameRate = 5;
open(v);
disp('Generating video...')
s_only = 1;
for k=1:idp
    if generateFigs == 1
        if opt.sc == 2, sealId = 2:2:22;
        else,           sealId = [2 22];
        end
        if strcmp(injCase, 'nextToF'),          lims = [1.2 1.53];
        elseif strcmp(injCase, 'nextToF_FW'),   lims = [1.1 1.4];
        end
        
        % Make tiled plot
        hf = figure(37);
        c_sh = [0.5 0.5 0.5];
        id = false(G.cells.num, 1);
        id([ucids.unit_cell_ids{1:22} ucids.unit_cell_ids{25:32} ucids.unit_cell_ids{58:59}]) = true;
        resId = 1;
        wellLoc = G.cells.centroids(W.cells,:);
        id(G.cells.centroids(:,1) < wellLoc(1)) = false;
        id(G.cells.centroids(:,1) > wellLoc(1)) = false;
        id(G.cells.centroids(:,2) < 9000) = false; id(G.cells.centroids(:,2) > 1.7*10^4) = false;
        %id(G.cells.centroids(:,3) < 10^3) = false; id(G.cells.centroids(:,3) > 4*10^3) = false;
        ids = false(G.cells.num, 1);  ids([ucids.unit_cell_ids{sealId}]) = id([ucids.unit_cell_ids{sealId}]);
        idf = false(G.cells.num, 1);  idf(ucids.unit_cell_ids{end}) = id(ucids.unit_cell_ids{end});
        if s_only == 1
            tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
            ax2 = nexttile; % -------- 1
            plotToolbar(G, states{k}.s(:,2), ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
                'lineWidth', 0.5, 'EdgeAlpha', 1) % seal edges
            hold on
            plotToolbar(G, states{k}.s(:,2), idf, 'FaceAlpha', 0, 'EdgeColor', 'none') % fault edges
            colormap(flipud(cmocean('deep')))
            plotToolbar(G, states{k}.s(:,2), id, 'FaceAlpha', 0.85);                   % full states plot
            set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
            plotWell(G, W, 'color', 'k', 'height', -1000);
            text(22000,14000, 1400, 'TS', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            text(22000,14000, 1800, 'SR', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            text(22000,11800, 1800, ['t = ' num2str(round(sum(timesteps(1:k))/year, 1)) ...
                ' y'], 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            c = colorbar; caxis([0 0.7]); c.Label.String = '$S_g$ [-]';
            c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
            c.Ticks = 0:.1:1;
            ylim(lims*10^4)
            zlim([1300 2500])
            axis off
            
            ax4 = nexttile; % --------- 4
            id2 = false(G.cells.num, 1);
            id2(ucids.unit_cell_ids{resId}) = true;
            plotToolbar(G, states{k}.s(:,2), id2, 'EdgeAlpha', 0.02, ...
                'lineWidth', 0.1, 'EdgeColor', [0.5 0.5 0.5]);
            set(gca,'Xdir','reverse'), view(55, 25),
            camproj perspective; axis equal tight off
            plotWell(G, W, 'color', 'w', 'height', -1500); 
            %c = colorbar; 
            caxis([0 0.7])
            %c.Ticks = 0:.1:1;
            colormap(ax4, flipud(cmocean('tempo')));
            colormap(ax2, flipud(cmocean('tempo')));
            set(hf, 'units', 'pixels', 'position', [200, 20, 1500, 600]);
            set(hf,'color','w');
        else
            tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
            ax = nexttile; % ------- 1
            plotToolbar(G, states{k}.dp/barsa, ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
                'lineWidth', 0.5, 'EdgeAlpha', 1); % seal edges
            hold on
            plotToolbar(G, states{k}.dp/barsa, idf, 'FaceAlpha', 0, ...
                'EdgeColor', 'none'); % fault edges
            plotToolbar(G, states{k}.dp/barsa, id, 'FaceAlpha', 0.85);                   % full states plot
            set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
            plotWell(G, W, 'color', 'k', 'height', -1000);
            text(22000,14000, 1400, 'TS', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            text(22000,14000, 1800, 'SR', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            text(22000,11800, 1800, ['t = ' num2str(round(sum(timesteps(1:k))/year, 1)) ...
                ' y'], 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            c = colorbar; caxis([0 10]); c.Label.String = '$\Delta p$ [bar]';
            c.Ticks = 0:1:10;
            c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
            ylim(lims*10^4)
            zlim([1300 2500])

            ax2 = nexttile; % -------- 2
            plotToolbar(G, states{k}.s(:,2), ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
                'lineWidth', 0.5, 'EdgeAlpha', 1) % seal edges
            hold on
            plotToolbar(G, states{k}.s(:,2), idf, 'FaceAlpha', 0, 'EdgeColor', 'none') % fault edges
            colormap(flipud(cmocean('deep')))
            plotToolbar(G, states{k}.s(:,2), id, 'FaceAlpha', 0.85);                   % full states plot
            set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
            plotWell(G, W, 'color', 'k', 'height', -1000);
            %text(22000,14000, 1400, 'TS', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            %text(22000,14000, 1800, 'SR', 'fontSize', 12, 'color', 'w', 'fontweight', 'bold')
            c = colorbar; caxis([0 0.7]); c.Label.String = '$S_g$ [-]';
            c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
            c.Ticks = 0:.1:1;
            ylim(lims*10^4)
            zlim([1300 2500])
            axis off
            colormap(ax2, flipud(cmocean('deep')));

            ax3 = nexttile; % ----------- 3
            id2 = false(G.cells.num, 1);
            id2(ucids.unit_cell_ids{resId}) = true;
            %id2(G.cells.centroids(:,2) < 10^4) = false;
            %id2(G.cells.centroids(:,2) > 1.7*10^4) = false;
            colormap(cmocean('balance'))
            plotToolbar(G, states{k}.dp/barsa, id2, 'EdgeAlpha', 0.02, ...
                'lineWidth', 0.1, 'EdgeColor', [0.5 0.5 0.5]);
            set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
            plotWell(G, W, 'color', 'w', 'height', -1500); c = colorbar; 
            caxis([0 10]); c.Ticks = 0:1:10;

            ax4 = nexttile; % --------- 4
            id2 = false(G.cells.num, 1);
            id2(ucids.unit_cell_ids{resId}) = true;
            plotToolbar(G, states{k}.s(:,2), id2, 'EdgeAlpha', 0.02, ...
                'lineWidth', 0.1, 'EdgeColor', [0.5 0.5 0.5]);
            set(gca,'Xdir','reverse'), view(55, 25),
            camproj perspective; axis equal tight off
            plotWell(G, W, 'color', 'w', 'height', -1500); c = colorbar; 
            caxis([0 0.7])
            c.Ticks = 0:.1:1;
            colormap(ax4, flipud(cmocean('tempo')));
            colormap(ax, cmocean('balance'));
            colormap(ax2, flipud(cmocean('tempo')));
            colormap(ax3, cmocean('balance'));
            set(hf, 'units', 'pixels', 'position', [200, 20, 1000, 500]);
            set(hf,'color','w');
        end

        % Export graphics
        %exportgraphics(hf,[dir num2str(k) '.jpg'],'ContentType','image',...
        %               'Resolution', 300, 'BackgroundColor','w')
        print(hf, [dir num2str(k)], '-djpeg', '-r300')
        close(hf)
    end
    % Read graphics and save frame
    I = imread([dir num2str(k) '.jpg']);
    
    % Add frame to video
    %frame = getframe(gcf);
    frame = im2frame(I);
    writeVideo(v,frame);
    close(gcf)
    
    if mod(k, 5) == 0
        disp([num2str(k) '/' num2str(idp) ' frames completed.'])
    end
end
close(v)
disp(['Done. Video file saved to ' dir]);  

% Fault property plots
zc = [1355 1871; 1450 1974];
clrs = copper(8);
latx = {'interpreter', 'latex'};
[zcord, ids] = sort(G.cells.centroids(fault.fcells, end));
hf = figure(27);
tiledlayout(2, 2, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile
hold on
fill([0 0.6 0.6 0 0], [zc(2,1) zc(2,1) zc(1,2) zc(1,2) zc(2,1)], clrs(1,:), ...
     'edgecolor', 'none', 'facealpha', 0.5)
fill([0 0.6 0.6 0 0], [zc(2,1) zc(2,1) zc(1,1) zc(1,1) zc(2,1)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 0.6 0.6 0 0], [zc(2,2) zc(2,2) zc(1,2) zc(1,2) zc(2,2)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 0.6 0.6 0 0], [zc(1,1) zc(1,1) 0 0 zc(1,1)], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
fill([0 0.6 0.6 0 0], [3000 3000 zc(2,2) zc(2,2) 3000], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
%plot([0 0.5], [zc(1) zc(1)], '--', 'color', [0.6 0.6 0.6])
%plot([0 0.5], [zc(1,2) zc(1,2)], '--', 'color', [0.6 0.6 0.6])
%plot([0 0.5], [zc(2,1) zc(2,1)], '-', 'color', [0.6 0.6 0.6])
%plot([0 0.5], [zc(2,2) zc(2,2)], '-', 'color', [0.6 0.6 0.6])
plot(fault.kPred.SGR(ids), zcord, '-k', 'linewidth', 1)
text(0.025, 2200, 'Storage', 'fontsize', 11)
text(0.025, 2450, 'reservoir (SR)', 'fontsize', 11)
text(0.025, 1650, 'Top seal (TS)', 'fontsize', 11)
hold off
set(gca,'YDir','reverse')
xlim([0 0.6]); xticks(0:.1:.6)
ylim([0 3000]);
grid on
xlabel('SGR', 'fontsize', 12)
ylabel('$z$ [m]', latx{:}, 'fontsize', 14)
nexttile % ------------
hold on
fill([-3.5 2 2 -3.5 -3.5], [zc(2,1) zc(2,1) zc(1,2) zc(1,2) zc(2,1)], clrs(1,:), ...
     'edgecolor', 'none', 'facealpha', 0.5)
fill([-3.5 2 2 -3.5 -3.5], [zc(2,1) zc(2,1) zc(1,1) zc(1,1) zc(2,1)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([-3.5 2 2 -3.5 -3.5], [zc(2,2) zc(2,2) zc(1,2) zc(1,2) zc(2,2)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([-3.5 2 2 -3.5 -3.5], [zc(1,1) zc(1,1) 0 0 zc(1,1)], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
fill([-3.5 2 2 -3.5 -3.5], [3000 3000 zc(2,2) zc(2,2) 3000], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
plot(log10(fault.k.vals(ids)/(milli*darcy)), zcord, '-k', 'linewidth', 1)
text(1.3, 2200, 'SR', 'fontsize', 10)
text(1.3, 1650, 'TS', 'fontsize', 10)
hold off
set(gca,'YDir','reverse')
xlim([-3.5 2]); xticks(-3:1:2)
ylim([0 3000]);
grid on
xlabel('$\log_{10}(k$ [mD])', latx{:},'fontsize', 12)
%ylabel('$z$ [m]', latx{:}, 'fontsize', 14)
nexttile % ------------
hold on
fill([0 0.5 0.5 0 0], [zc(2,1) zc(2,1) zc(1,2) zc(1,2) zc(2,1)], clrs(1,:), ...
     'edgecolor', 'none', 'facealpha', 0.5)
fill([0 0.5 0.5 0 0], [zc(2,1) zc(2,1) zc(1,1) zc(1,1) zc(2,1)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 0.5 0.5 0 0], [zc(2,2) zc(2,2) zc(1,2) zc(1,2) zc(2,2)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 0.5 0.5 0 0], [zc(1,1) zc(1,1) 0 0 zc(1,1)], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
fill([0 0.5 0.5 0 0], [3000 3000 zc(2,2) zc(2,2) 3000], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
plot(fault.poro(ids), zcord, '-k', 'linewidth', 1)
text(0.33, 2200, 'SR', 'fontsize', 10)
text(0.33, 1650, 'TS', 'fontsize', 10)
hold off
set(gca,'YDir','reverse')
xlim([0 0.4]); xticks(0:.1:.4)
ylim([0 3000]);
grid on
xlabel('$\phi$ [-]', latx{:},'fontsize', 12)
ylabel('$z$ [m]', latx{:}, 'fontsize', 14)
nexttile % ------------
hold on
fill([0 12 12 0 0], [zc(2,1) zc(2,1) zc(1,2) zc(1,2) zc(2,1)], clrs(1,:), ...
     'edgecolor', 'none', 'facealpha', 0.5)
fill([0 12 12 0 0], [zc(2,1) zc(2,1) zc(1,1) zc(1,1) zc(2,1)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 12 12 0 0], [zc(2,2) zc(2,2) zc(1,2) zc(1,2) zc(2,2)], clrs(3,:), ...
     'edgecolor', 'none', 'facealpha', 0.35)
fill([0 12 12 0 0], [zc(1,1) zc(1,1) 0 0 zc(1,1)], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
fill([0 12 12 0 0], [3000 3000 zc(2,2) zc(2,2) 3000], clrs(8,:), ...
     'edgecolor', 'none', 'facealpha', 0.2)
pcv =fluid.pcOG{3}(1e-3*ones(numel(fault.fcells), 1))/barsa;
plot(pcv(ids), zcord, '-k', 'linewidth', 1)
text(0.85, 2200, 'SR', 'fontsize', 10)
text(0.85, 1650, 'TS', 'fontsize', 10)
hold off
set(gca,'YDir','reverse')
xlim([0 1]); xticks(0:0.2:1)
ylim([0 3000]);
grid on
xlabel('$P_\mathrm{e}$ [bar]', latx{:},'fontsize', 12)
%ylabel('$z$ [m]', latx{:}, 'fontsize', 14)
set(hf, 'position', [100, 100, 550, 500])

% Performance plots
plotPerf_bo(report);

% Line / point plots
if strcmp(injCase, 'nextToF')
cellsClosestTo = [22500, 13583, 1900;   % top reservoir cell
                  22500, 13583, 1880;   % caprock, above top res cell
                  22500, 12920, 1980];  % fault, contact SR (HW) with TS
elseif strcmp(injCase, 'nextToF_FW')
cellsClosestTo = [22500, 12755, 1880;   % top reservoir cell, approx above inj
                  22500, 12755, 1865;   % caprock, above top res cell
                  22500, 12780, 1875];  % fault, contact SR (FW) with TS    
end
monitoringCells = plotLine_bo(G, timesteps, W, states, state0, bhp, ...
                              cellsClosestTo, ucids, injCase);